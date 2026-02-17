#!/usr/bin/env python3
"""
18_mapbiomas_validation.py
===========================
Cross-validates the study's Olofsson forest area estimates against
MapBiomas Colombia Collection 1 as an independent corroboration source.

Purpose:
    Verify that the direction and approximate magnitude of forest decline
    detected by the study's Random-Forest / Olofsson pipeline are consistent
    with an independent, continent-wide land-cover product.

GEE Asset:
    projects/mapbiomas-colombia/assets/collection1/mapbiomas_colombia_collection1_integration_v1

MapBiomas years used:
    2013 (T1), 2016 (T2), 2020 (T3), 2022 (closest proxy for T4=2024,
    since MapBiomas Colombia Collection 1 may not cover 2024).

Reclassification to 5-class scheme:
    Forest   : MapBiomas codes 1-6, 49
    Pastures : MapBiomas codes 12, 15, 18, 19, 20, 21, 36, 39, 41, 46, 47, 48
    Urban    : MapBiomas code  24
    Water    : MapBiomas codes 26, 33
    Other    : all remaining codes

Inputs:
    - gee_config.py          -> STUDY_AREA_BBOX (ee.Geometry.Rectangle)
    - outputs/phase3_stats/olofsson_area_estimates.json

Outputs:
    - outputs/phase3_stats/mapbiomas_validation.json

Requirements:
    - Google Earth Engine Python API (`earthengine-api`)
    - Authenticated GEE session (run `earthengine authenticate` first)
    - GEE project ID set in .env as GEE_PROJECT_ID

Usage:
    python scripts/18_mapbiomas_validation.py

References:
    MapBiomas Colombia. (2024). Collection 1 of annual land use and land
        cover maps. https://colombia.mapbiomas.org
    Olofsson, P. et al. (2014). Good practices for estimating area and
        assessing accuracy of land change. RSE, 148, 42-57.
"""

import os
import sys
import json
import time
import math
from datetime import datetime

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'outputs', 'phase3_stats')
OLOFSSON_PATH = os.path.join(OUTPUT_DIR, 'olofsson_area_estimates.json')
OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'mapbiomas_validation.json')

# ---------------------------------------------------------------------------
# GEE initialisation (must happen before importing gee_config which also
# calls ee.Initialize)
# ---------------------------------------------------------------------------

sys.path.insert(0, BASE_DIR)

import ee  # noqa: E402

# gee_config.py initialises EE and exports STUDY_AREA_BBOX
from gee_config import STUDY_AREA_BBOX  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAPBIOMAS_ASSET = (
    'projects/mapbiomas-public/assets/colombia/collection1/'
    'mapbiomas_colombia_collection1_integration_v1'
)

# MapBiomas band names follow the pattern 'classification_YYYY'
MAPBIOMAS_YEARS = [2013, 2016, 2020, 2022]

# Study period mapping (Olofsson keys -> MapBiomas year)
PERIOD_TO_MB_YEAR = {
    'pre_acuerdo':    2013,
    'transicion':     2016,
    'post_acuerdo_1': 2020,
    'post_acuerdo_2': 2022,   # proxy: 2022 for T4 (2024 may not be available)
}

# Pixel resolution for area calculation (MapBiomas Colombia = 30 m)
MAPBIOMAS_SCALE = 30  # metres

# Bounding box (must match gee_config.py)
BBOX = [-75.0, 6.0, -73.5, 8.0]

# ---------------------------------------------------------------------------
# MapBiomas Colombia legend reclassification
# ---------------------------------------------------------------------------

# Target classes
CLASS_FOREST   = 1
CLASS_PASTURES = 2
CLASS_URBAN    = 3
CLASS_WATER    = 4
CLASS_OTHER    = 5

CLASS_NAMES = {
    CLASS_FOREST:   'Forest',
    CLASS_PASTURES: 'Pastures',
    CLASS_URBAN:    'Urban',
    CLASS_WATER:    'Water',
    CLASS_OTHER:    'Other',
}

# MapBiomas source codes -> target class
FOREST_CODES   = [1, 2, 3, 4, 5, 6, 49]
PASTURE_CODES  = [12, 15, 18, 19, 20, 21, 36, 39, 41, 46, 47, 48]
URBAN_CODES    = [24]
WATER_CODES    = [26, 33]

# Build a lookup dict: mapbiomas_code -> target_class
_RECLASS_MAP = {}
for code in FOREST_CODES:
    _RECLASS_MAP[code] = CLASS_FOREST
for code in PASTURE_CODES:
    _RECLASS_MAP[code] = CLASS_PASTURES
for code in URBAN_CODES:
    _RECLASS_MAP[code] = CLASS_URBAN
for code in WATER_CODES:
    _RECLASS_MAP[code] = CLASS_WATER


def build_reclassify_image(source_image):
    """
    Reclassify a single-band MapBiomas classification image into the
    5-class scheme using ee.Image.remap().

    Any code not in the explicit lists is mapped to CLASS_OTHER.
    """
    from_values = []
    to_values = []
    for code, target in _RECLASS_MAP.items():
        from_values.append(code)
        to_values.append(target)

    reclassified = source_image.remap(
        from_values, to_values, defaultValue=CLASS_OTHER
    ).rename('class').toInt8()

    return reclassified


# ---------------------------------------------------------------------------
# Zonal area extraction
# ---------------------------------------------------------------------------

def extract_class_areas_ha(classified_image, region, scale):
    """
    Compute per-class area (ha) within *region* by counting pixels.

    Uses ee.Image.reduceRegion with a grouped frequency histogram.
    Band order: area (band 0) then class (band 1) — the grouped reducer
    sums band 0 and groups by band 1.

    Returns:
        dict  {class_int: area_ha}  (server-side values fetched to client)
    """
    # Pixel area image in hectares
    pixel_area_ha = ee.Image.pixelArea().divide(10000).rename('area')

    # Band order matters: area first (summed), class second (grouped)
    combined = pixel_area_ha.addBands(classified_image)

    # Sum pixel area per class
    result = combined.reduceRegion(
        reducer=ee.Reducer.sum().group(
            groupField=1,
            groupName='class'
        ),
        geometry=region,
        scale=scale,
        maxPixels=1e10,
        bestEffort=True,
    )

    # Parse grouped output
    groups = ee.List(result.get('groups')).getInfo()
    area_dict = {}
    for g in groups:
        cls = int(g['class'])
        area = float(g['sum'])
        area_dict[cls] = round(area, 1)

    return area_dict


# ---------------------------------------------------------------------------
# Load Olofsson estimates
# ---------------------------------------------------------------------------

def load_olofsson_forest_areas(path):
    """
    Read per-period total forest area (Bosque denso + Bosque secundario)
    from the Olofsson JSON, along with 95 % CIs.

    Returns:
        dict  {period_key: {year, forest_ha, forest_lower_ha, forest_upper_ha,
                            forest_se_ha}}
    """
    with open(path, 'r') as f:
        data = json.load(f)

    out = {}
    for period_key, pdata in data['periods'].items():
        cls1 = pdata['per_class']['1']   # Bosque denso
        cls2 = pdata['per_class']['2']   # Bosque secundario

        forest_ha = cls1['adjusted_area_ha'] + cls2['adjusted_area_ha']
        forest_se = math.sqrt(cls1['area_se_ha']**2 + cls2['area_se_ha']**2)
        forest_ci95 = 1.96 * forest_se

        out[period_key] = {
            'year': pdata['year'],
            'label': pdata['label'],
            'forest_ha': round(forest_ha, 0),
            'forest_se_ha': round(forest_se, 0),
            'forest_ci95_ha': round(forest_ci95, 0),
            'forest_lower_ha': round(forest_ha - forest_ci95, 0),
            'forest_upper_ha': round(forest_ha + forest_ci95, 0),
            'bosque_denso_ha': cls1['adjusted_area_ha'],
            'bosque_secundario_ha': cls2['adjusted_area_ha'],
        }

    return out


# ---------------------------------------------------------------------------
# Concordance metrics
# ---------------------------------------------------------------------------

def compute_concordance(mapbiomas_forest, olofsson_forest):
    """
    Compare forest area time series between MapBiomas and Olofsson.

    Parameters:
        mapbiomas_forest : list of (year, forest_ha)  -- chronologically sorted
        olofsson_forest  : list of (year, forest_ha, forest_lower, forest_upper)

    Returns:
        dict with concordance statistics
    """
    # --- Direction agreement ---
    # For each consecutive pair, check if both sources agree on direction
    n_intervals = len(mapbiomas_forest) - 1
    direction_agree = 0

    mb_directions = []
    ol_directions = []
    for i in range(n_intervals):
        mb_delta = mapbiomas_forest[i + 1][1] - mapbiomas_forest[i][1]
        ol_delta = olofsson_forest[i + 1][1] - olofsson_forest[i][1]

        mb_dir = 'increase' if mb_delta > 0 else ('decrease' if mb_delta < 0 else 'stable')
        ol_dir = 'increase' if ol_delta > 0 else ('decrease' if ol_delta < 0 else 'stable')

        mb_directions.append(mb_dir)
        ol_directions.append(ol_dir)

        if mb_dir == ol_dir:
            direction_agree += 1

    direction_agreement_pct = (direction_agree / n_intervals * 100) if n_intervals > 0 else None

    # --- Overall trend agreement ---
    mb_overall_change = mapbiomas_forest[-1][1] - mapbiomas_forest[0][1]
    ol_overall_change = olofsson_forest[-1][1] - olofsson_forest[0][1]
    overall_direction_agree = (
        (mb_overall_change > 0 and ol_overall_change > 0) or
        (mb_overall_change < 0 and ol_overall_change < 0) or
        (mb_overall_change == 0 and ol_overall_change == 0)
    )

    # --- Magnitude comparison ---
    mb_total_forest_first = mapbiomas_forest[0][1]
    ol_total_forest_first = olofsson_forest[0][1]
    mb_pct_change = (mb_overall_change / mb_total_forest_first * 100) if mb_total_forest_first > 0 else 0
    ol_pct_change = (ol_overall_change / ol_total_forest_first * 100) if ol_total_forest_first > 0 else 0

    # --- Does MapBiomas forest fall within Olofsson 95% CI? ---
    ci_overlap = []
    for i, (mb_year, mb_ha) in enumerate(mapbiomas_forest):
        ol_lower = olofsson_forest[i][2]
        ol_upper = olofsson_forest[i][3]
        within_ci = ol_lower <= mb_ha <= ol_upper
        ci_overlap.append({
            'mb_year': mb_year,
            'ol_year': olofsson_forest[i][0],
            'mb_forest_ha': round(mb_ha, 0),
            'ol_forest_ha': round(olofsson_forest[i][1], 0),
            'ol_lower_ha': round(ol_lower, 0),
            'ol_upper_ha': round(ol_upper, 0),
            'within_olofsson_ci95': within_ci,
        })

    n_within = sum(1 for c in ci_overlap if c['within_olofsson_ci95'])

    return {
        'n_intervals': n_intervals,
        'direction_agreement_count': direction_agree,
        'direction_agreement_pct': round(direction_agreement_pct, 1) if direction_agreement_pct is not None else None,
        'interval_directions': [
            {
                'from_year': mapbiomas_forest[i][0],
                'to_year': mapbiomas_forest[i + 1][0],
                'mapbiomas_direction': mb_directions[i],
                'olofsson_direction': ol_directions[i],
                'agree': mb_directions[i] == ol_directions[i],
            }
            for i in range(n_intervals)
        ],
        'overall_trend': {
            'mapbiomas_change_ha': round(mb_overall_change, 0),
            'olofsson_change_ha': round(ol_overall_change, 0),
            'mapbiomas_pct_change': round(mb_pct_change, 2),
            'olofsson_pct_change': round(ol_pct_change, 2),
            'direction_agree': overall_direction_agree,
        },
        'ci_overlap': ci_overlap,
        'n_within_olofsson_ci95': n_within,
        'pct_within_olofsson_ci95': round(n_within / len(ci_overlap) * 100, 1) if ci_overlap else None,
    }


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    t0 = time.time()
    print('=' * 70)
    print('MAPBIOMAS COLOMBIA CROSS-VALIDATION')
    print(f'Date : {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'Asset: {MAPBIOMAS_ASSET}')
    print(f'BBOX : {BBOX}')
    print('=' * 70)

    # ------------------------------------------------------------------
    # 1. Load MapBiomas image
    # ------------------------------------------------------------------
    print('\n[1/5] Loading MapBiomas Colombia Collection 1 ...')

    try:
        mapbiomas = ee.Image(MAPBIOMAS_ASSET)
        band_names = mapbiomas.bandNames().getInfo()
        print(f'  Bands available: {len(band_names)}')
        print(f'  First 5: {band_names[:5]}')
        print(f'  Last 5 : {band_names[-5:]}')
    except Exception as e:
        print(f'  ERROR loading MapBiomas asset: {e}')
        print('  Make sure you have access to the MapBiomas Colombia asset.')
        print('  Asset ID: ' + MAPBIOMAS_ASSET)
        sys.exit(1)

    region = STUDY_AREA_BBOX

    # ------------------------------------------------------------------
    # 2. Extract and reclassify per year
    # ------------------------------------------------------------------
    print('\n[2/5] Extracting and reclassifying MapBiomas per year ...')

    mapbiomas_results = {}

    for year in MAPBIOMAS_YEARS:
        band_name = f'classification_{year}'
        print(f'\n  --- {year} (band: {band_name}) ---')

        if band_name not in band_names:
            print(f'  WARNING: Band {band_name} not found in asset. Skipping.')
            mapbiomas_results[year] = None
            continue

        try:
            classification = mapbiomas.select(band_name)
            reclassified = build_reclassify_image(classification)

            print('  Computing zonal areas ...')
            areas = extract_class_areas_ha(reclassified, region, MAPBIOMAS_SCALE)

            # Ensure all classes present
            for cls_id in CLASS_NAMES:
                if cls_id not in areas:
                    areas[cls_id] = 0.0

            total = sum(areas.values())
            print(f'  Total area: {total:,.0f} ha')
            for cls_id in sorted(areas.keys()):
                pct = areas[cls_id] / total * 100 if total > 0 else 0
                print(f'    {CLASS_NAMES.get(cls_id, f"Class {cls_id}"):<12}: '
                      f'{areas[cls_id]:>12,.0f} ha  ({pct:5.1f}%)')

            mapbiomas_results[year] = {
                'year': year,
                'band': band_name,
                'area_ha': {CLASS_NAMES[k]: v for k, v in sorted(areas.items())},
                'total_ha': round(total, 1),
                'forest_ha': round(areas.get(CLASS_FOREST, 0), 0),
                'pastures_ha': round(areas.get(CLASS_PASTURES, 0), 0),
                'urban_ha': round(areas.get(CLASS_URBAN, 0), 0),
                'water_ha': round(areas.get(CLASS_WATER, 0), 0),
                'other_ha': round(areas.get(CLASS_OTHER, 0), 0),
            }

        except Exception as e:
            print(f'  ERROR processing year {year}: {e}')
            mapbiomas_results[year] = None

    # ------------------------------------------------------------------
    # 3. Load Olofsson estimates
    # ------------------------------------------------------------------
    print('\n[3/5] Loading Olofsson area estimates ...')

    if not os.path.isfile(OLOFSSON_PATH):
        print(f'  ERROR: Olofsson file not found: {OLOFSSON_PATH}')
        sys.exit(1)

    olofsson = load_olofsson_forest_areas(OLOFSSON_PATH)
    for pk, od in olofsson.items():
        print(f'  {od["label"]}: {od["forest_ha"]:,.0f} ha '
              f'[{od["forest_lower_ha"]:,.0f} - {od["forest_upper_ha"]:,.0f}]')

    # ------------------------------------------------------------------
    # 4. Concordance analysis
    # ------------------------------------------------------------------
    print('\n[4/5] Computing concordance metrics ...')

    # Build parallel time series
    period_keys_ordered = ['pre_acuerdo', 'transicion', 'post_acuerdo_1', 'post_acuerdo_2']

    mapbiomas_series = []
    olofsson_series = []
    valid_periods = []

    for pk in period_keys_ordered:
        mb_year = PERIOD_TO_MB_YEAR[pk]
        mb_data = mapbiomas_results.get(mb_year)
        ol_data = olofsson.get(pk)

        if mb_data is None or ol_data is None:
            print(f'  Skipping {pk}: missing data')
            continue

        mapbiomas_series.append((mb_year, mb_data['forest_ha']))
        olofsson_series.append((
            ol_data['year'],
            ol_data['forest_ha'],
            ol_data['forest_lower_ha'],
            ol_data['forest_upper_ha'],
        ))
        valid_periods.append(pk)

    if len(mapbiomas_series) < 2:
        print('  ERROR: Not enough valid periods for concordance analysis.')
        concordance = {'error': 'Insufficient data (fewer than 2 valid periods)'}
    else:
        concordance = compute_concordance(mapbiomas_series, olofsson_series)

        print(f'\n  Direction agreement: {concordance["direction_agreement_count"]}/'
              f'{concordance["n_intervals"]} intervals '
              f'({concordance["direction_agreement_pct"]}%)')
        print(f'  Overall trend agrees: {concordance["overall_trend"]["direction_agree"]}')
        print(f'    MapBiomas : {concordance["overall_trend"]["mapbiomas_change_ha"]:+,.0f} ha '
              f'({concordance["overall_trend"]["mapbiomas_pct_change"]:+.2f}%)')
        print(f'    Olofsson  : {concordance["overall_trend"]["olofsson_change_ha"]:+,.0f} ha '
              f'({concordance["overall_trend"]["olofsson_pct_change"]:+.2f}%)')
        print(f'  MapBiomas within Olofsson 95% CI: '
              f'{concordance["n_within_olofsson_ci95"]}/{len(concordance["ci_overlap"])} periods '
              f'({concordance["pct_within_olofsson_ci95"]}%)')

    # ------------------------------------------------------------------
    # 5. Write output JSON
    # ------------------------------------------------------------------
    print('\n[5/5] Writing output ...')

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    output = {
        'metadata': {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'purpose': (
                'Independent cross-validation of study forest area estimates '
                'against MapBiomas Colombia Collection 1'
            ),
            'mapbiomas_asset': MAPBIOMAS_ASSET,
            'study_bbox': BBOX,
            'mapbiomas_scale_m': MAPBIOMAS_SCALE,
            'reclassification': {
                'Forest': f'MapBiomas codes {FOREST_CODES}',
                'Pastures': f'MapBiomas codes {PASTURE_CODES}',
                'Urban': f'MapBiomas codes {URBAN_CODES}',
                'Water': f'MapBiomas codes {WATER_CODES}',
                'Other': 'All remaining codes',
            },
            'note_t4_proxy': (
                'MapBiomas year 2022 used as proxy for study T4 (2024) '
                'because Collection 1 may not include 2024.'
            ),
            'olofsson_source': OLOFSSON_PATH,
        },
        'mapbiomas_areas': {
            str(year): data
            for year, data in mapbiomas_results.items()
            if data is not None
        },
        'olofsson_forest': {
            pk: olofsson[pk]
            for pk in period_keys_ordered
            if pk in olofsson
        },
        'concordance': concordance,
    }

    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f'  Written: {OUTPUT_PATH}')

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)
    print(f'  MapBiomas years processed : {[y for y in MAPBIOMAS_YEARS if mapbiomas_results.get(y) is not None]}')
    print(f'  Valid comparison periods   : {len(valid_periods)}')
    if isinstance(concordance, dict) and 'overall_trend' in concordance:
        print(f'  Overall direction agree   : {concordance["overall_trend"]["direction_agree"]}')
        print(f'  Interval direction agree  : {concordance["direction_agreement_pct"]}%')
        print(f'  Within Olofsson 95% CI    : {concordance["pct_within_olofsson_ci95"]}%')
    print(f'  Output file               : {OUTPUT_PATH}')
    print(f'  Elapsed                   : {elapsed:.1f} s')
    print('=' * 70)


if __name__ == '__main__':
    main()
