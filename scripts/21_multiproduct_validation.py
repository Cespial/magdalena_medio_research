#!/usr/bin/env python3
"""
21_multiproduct_validation.py
==============================
Cross-validates the study's Olofsson forest area estimates against multiple
independent global/regional land-cover products in Google Earth Engine.

Purpose:
    Establish multi-source concordance for the direction and magnitude of
    forest change detected by the study's Random-Forest / Olofsson pipeline.
    Directional agreement across 5+ independent products strengthens the
    conclusion that observed trends are robust, even when absolute areas differ
    due to resolution, classification scheme, and temporal coverage.

Products validated:
    1. MapBiomas Colombia Collection 1 (loaded from prior script output)
    2. ESA WorldCover v200 (10 m, 2021 only)
    3. Google Dynamic World v1 (10 m, 2015-present)
    4. MODIS MCD12Q1 v061 (500 m, 2013-2023)
    5. Hansen Global Forest Change v1.12 (30 m, 2000-2024)

Inputs:
    - gee_config.py                -> STUDY_AREA_BBOX, PERIODS
    - outputs/phase3_stats/olofsson_area_estimates.json
    - outputs/phase3_stats/mapbiomas_validation.json

Outputs:
    - outputs/phase3_stats/multiproduct_validation.json

Requirements:
    - Google Earth Engine Python API (`earthengine-api`)
    - Authenticated GEE session (run `earthengine authenticate` first)
    - GEE project ID set in .env as GEE_PROJECT_ID

Usage:
    python scripts/21_multiproduct_validation.py

References:
    Zanaga, D., et al. (2022). ESA WorldCover 10 m 2021 v200.
    Brown, C. F., et al. (2022). Dynamic World. Nature Scientific Data.
    Friedl, M. A., et al. (2022). MODIS/Terra+Aqua Land Cover Type
        Yearly L3 Global 500 m, MCD12Q1 v061.
    Hansen, M. C., et al. (2013). High-resolution global maps of 21st-century
        forest cover change. Science, 342(6160), 850-853.
    MapBiomas Colombia. (2024). Collection 1 of annual land use and land
        cover maps. https://colombia.mapbiomas.org
    Olofsson, P., et al. (2014). Good practices for estimating area and
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
MAPBIOMAS_PATH = os.path.join(OUTPUT_DIR, 'mapbiomas_validation.json')
OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'multiproduct_validation.json')

# ---------------------------------------------------------------------------
# GEE initialisation
# ---------------------------------------------------------------------------

sys.path.insert(0, BASE_DIR)

import ee  # noqa: E402

from gee_config import STUDY_AREA_BBOX, PERIODS  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BBOX = [-75.0, 6.0, -73.5, 8.0]

# Ordered period keys and their target years
PERIOD_KEYS = ['pre_acuerdo', 'transicion', 'post_acuerdo_1', 'post_acuerdo_2']
TARGET_YEARS = {pk: PERIODS[pk]['map_year'] for pk in PERIOD_KEYS}
# -> {'pre_acuerdo': 2013, 'transicion': 2016, 'post_acuerdo_1': 2020, 'post_acuerdo_2': 2024}

# Interval labels for directional agreement
INTERVALS = [
    ('T1_T2', 'pre_acuerdo', 'transicion'),
    ('T2_T3', 'transicion', 'post_acuerdo_1'),
    ('T3_T4', 'post_acuerdo_1', 'post_acuerdo_2'),
]

# ---------------------------------------------------------------------------
# GEE product identifiers
# ---------------------------------------------------------------------------

ESA_WORLDCOVER_ID = 'ESA/WorldCover/v200'
DYNAMIC_WORLD_ID = 'GOOGLE/DYNAMICWORLD/V1'
MODIS_LC_ID = 'MODIS/061/MCD12Q1'
HANSEN_ID = 'UMD/hansen/global_forest_change_2024_v1_12'


# ===========================================================================
# Helper: pixel-area-based forest area computation
# ===========================================================================

def compute_forest_area_ha(forest_mask, region, scale):
    """
    Compute total forest area (ha) for a binary forest mask (1 = forest).

    Uses ee.Image.pixelArea() for geodetically accurate area rather than
    simple pixel counting.

    Parameters:
        forest_mask : ee.Image  -- single band, 1 = forest, 0 = non-forest
        region      : ee.Geometry
        scale       : int       -- nominal resolution in metres

    Returns:
        float  -- forest area in hectares
    """
    area_image = forest_mask.multiply(ee.Image.pixelArea()).divide(10000)

    result = area_image.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=region,
        scale=scale,
        maxPixels=1e10,
        bestEffort=True,
    )

    area_ha = result.values().get(0).getInfo()
    return round(float(area_ha), 0)


# ===========================================================================
# 1. Load Olofsson (this study) forest areas
# ===========================================================================

def load_this_study():
    """
    Read per-period total forest area (Bosque denso + Bosque secundario)
    from the Olofsson JSON.

    Returns:
        dict  {period_key: forest_ha}   keyed by period
    """
    with open(OLOFSSON_PATH, 'r') as f:
        data = json.load(f)

    out = {}
    for pk in PERIOD_KEYS:
        pdata = data['periods'][pk]
        cls1 = pdata['per_class']['1']   # Bosque denso
        cls2 = pdata['per_class']['2']   # Bosque secundario
        forest_ha = cls1['adjusted_area_ha'] + cls2['adjusted_area_ha']
        out[pk] = round(forest_ha, 0)

    return out


# ===========================================================================
# 2. Load MapBiomas (from prior script 18 output)
# ===========================================================================

def load_mapbiomas():
    """
    Load MapBiomas forest area from the mapbiomas_validation.json produced
    by 18_mapbiomas_validation.py.

    MapBiomas has data for 2013, 2016, 2020, 2022 (proxy for T4=2024).

    Returns:
        dict  {period_key: forest_ha}
    """
    with open(MAPBIOMAS_PATH, 'r') as f:
        data = json.load(f)

    mb_areas = data['mapbiomas_areas']

    # Mapping from period key to the MapBiomas year available
    mb_year_map = {
        'pre_acuerdo': '2013',
        'transicion': '2016',
        'post_acuerdo_1': '2020',
        'post_acuerdo_2': '2022',   # proxy for 2024
    }

    out = {}
    for pk, mb_yr in mb_year_map.items():
        if mb_yr in mb_areas and mb_areas[mb_yr] is not None:
            out[pk] = mb_areas[mb_yr]['forest_ha']
        else:
            out[pk] = None

    return out


# ===========================================================================
# 3. ESA WorldCover v200 (2021 only)
# ===========================================================================

def compute_esa_worldcover(region):
    """
    ESA WorldCover v200: 10 m global land cover for 2021.
    Class 10 = Tree cover -> Forest.

    Since only 2021 is available, this serves as a single reference point
    closest to T3 (2020).  We store the value under each period key but
    mark it as the same 2021 snapshot.

    Returns:
        dict  {period_key: forest_ha or None}
    """
    print('\n  [ESA WorldCover] Loading ESA/WorldCover/v200 ...')

    wc = ee.ImageCollection(ESA_WORLDCOVER_ID).first()

    # Class 10 = Tree cover
    forest_mask = wc.select('Map').eq(10)

    forest_ha = compute_forest_area_ha(forest_mask, region, scale=10)
    print(f'  [ESA WorldCover] Forest area (2021): {forest_ha:,.0f} ha')

    # Assign the 2021 value to the closest period (T3 = 2020), None for others
    out = {
        'pre_acuerdo': None,
        'transicion': None,
        'post_acuerdo_1': forest_ha,   # 2021 as proxy for T3=2020
        'post_acuerdo_2': None,
    }

    return out


# ===========================================================================
# 4. Google Dynamic World v1
# ===========================================================================

def compute_dynamic_world(region):
    """
    Google Dynamic World v1: 10 m near-real-time land cover (2015-present).
    Band 'label' contains per-pixel most-likely class.
    Class 1 = trees -> Forest.

    For each study period, create a temporal composite by taking the mode
    of the 'label' band over the period window defined in gee_config.PERIODS.

    Dynamic World starts mid-2015, so T1 (2013) is not available.

    Returns:
        dict  {period_key: forest_ha or None}
    """
    print('\n  [Dynamic World] Loading GOOGLE/DYNAMICWORLD/V1 ...')

    out = {}

    for pk in PERIOD_KEYS:
        year = TARGET_YEARS[pk]
        start = PERIODS[pk]['start']
        end = PERIODS[pk]['end']

        print(f'  [Dynamic World] Period {pk} ({start} to {end}) ...')

        # Dynamic World available from ~2015-06-23 onward
        if year < 2015:
            print(f'    Skipping: Dynamic World not available before 2015')
            out[pk] = None
            continue

        try:
            dw_col = (ee.ImageCollection(DYNAMIC_WORLD_ID)
                      .filterBounds(region)
                      .filterDate(start, end)
                      .select('label'))

            n_images = dw_col.size().getInfo()
            print(f'    Images in window: {n_images}')

            if n_images == 0:
                print(f'    WARNING: No Dynamic World images found for this period')
                out[pk] = None
                continue

            # Mode composite: most frequent class per pixel
            mode_composite = dw_col.mode()

            # Class 1 = trees
            forest_mask = mode_composite.eq(1)

            forest_ha = compute_forest_area_ha(forest_mask, region, scale=10)
            print(f'    Forest area: {forest_ha:,.0f} ha')
            out[pk] = forest_ha

        except Exception as e:
            print(f'    ERROR: {e}')
            out[pk] = None

    return out


# ===========================================================================
# 5. MODIS MCD12Q1 v061
# ===========================================================================

def compute_modis_lc(region):
    """
    MODIS MCD12Q1 v061: 500 m annual land cover (2001-2023).
    Uses LC_Type1 (IGBP classification).

    Forest classes in IGBP:
        1 = Evergreen Needleleaf Forests
        2 = Evergreen Broadleaf Forests
        3 = Deciduous Needleleaf Forests
        4 = Deciduous Broadleaf Forests
        5 = Mixed Forests
        8 = Woody Savannas
        9 = Savannas

    Annual products are keyed by year.  MODIS LC available through 2023,
    so T4 (2024) uses 2023 as proxy.

    Returns:
        dict  {period_key: forest_ha or None}
    """
    print('\n  [MODIS LC] Loading MODIS/061/MCD12Q1 ...')

    # IGBP forest classes
    forest_classes = [1, 2, 3, 4, 5, 8, 9]

    # Mapping from period key to the closest available MODIS year
    modis_year_map = {
        'pre_acuerdo': 2013,
        'transicion': 2016,
        'post_acuerdo_1': 2020,
        'post_acuerdo_2': 2023,   # proxy for 2024 (MODIS LC ends at 2023)
    }

    out = {}

    for pk in PERIOD_KEYS:
        modis_year = modis_year_map[pk]
        date_start = f'{modis_year}-01-01'
        date_end = f'{modis_year}-12-31'

        print(f'  [MODIS LC] Period {pk} (year {modis_year}) ...')

        try:
            modis_img = (ee.ImageCollection(MODIS_LC_ID)
                         .filterDate(date_start, date_end)
                         .first()
                         .select('LC_Type1'))

            # Build forest mask: pixel is forest if LC_Type1 in forest_classes
            # Use remap to set forest classes -> 1, everything else -> 0
            from_vals = list(range(0, 18))  # IGBP classes 0-17
            to_vals = [1 if v in forest_classes else 0 for v in from_vals]

            forest_mask = modis_img.remap(from_vals, to_vals, defaultValue=0)

            forest_ha = compute_forest_area_ha(forest_mask, region, scale=500)
            print(f'    Forest area: {forest_ha:,.0f} ha')
            out[pk] = forest_ha

        except Exception as e:
            print(f'    ERROR: {e}')
            out[pk] = None

    return out


# ===========================================================================
# 6. Hansen Global Forest Change v1.12
# ===========================================================================

def compute_hansen_gfc(region):
    """
    Hansen GFC v1.12: 30 m global forest change (2000-2024).

    Forest definition:
        Forest in year Y = treecover2000 > 25% AND no loss recorded
                           in any year from 2001 up to (Y - 2000).

    The 'lossyear' band encodes the year of loss as an integer offset
    from 2000 (e.g., 13 = 2013, 20 = 2020).

    Note: Hansen only tracks loss, not regrowth.  So forest area is
    monotonically non-increasing by definition.

    Returns:
        dict  {period_key: forest_ha or None}
    """
    print('\n  [Hansen GFC] Loading UMD/hansen/global_forest_change_2024_v1_12 ...')

    try:
        hansen = ee.Image(HANSEN_ID)

        treecover2000 = hansen.select('treecover2000')
        lossyear = hansen.select('lossyear')

        # Baseline: pixels with > 25% canopy cover in 2000
        baseline_forest = treecover2000.gt(25)

    except Exception as e:
        print(f'  ERROR loading Hansen image: {e}')
        return {pk: None for pk in PERIOD_KEYS}

    out = {}

    for pk in PERIOD_KEYS:
        year = TARGET_YEARS[pk]
        loss_offset = year - 2000  # e.g., 2013 -> 13

        print(f'  [Hansen GFC] Period {pk} (year {year}, loss offset <= {loss_offset}) ...')

        try:
            # Loss occurred before or during target year
            # lossyear > 0 ensures we only consider pixels where loss occurred
            # lossyear <= offset ensures loss was before or in the target year
            loss_before_year = lossyear.gt(0).And(lossyear.lte(loss_offset))

            # Forest = baseline AND no loss before target year
            forest_mask = baseline_forest.And(loss_before_year.Not())

            forest_ha = compute_forest_area_ha(forest_mask, region, scale=30)
            print(f'    Forest area: {forest_ha:,.0f} ha')
            out[pk] = forest_ha

        except Exception as e:
            print(f'    ERROR: {e}')
            out[pk] = None

    return out


# ===========================================================================
# Directional agreement analysis
# ===========================================================================

def classify_direction(area_start, area_end):
    """
    Classify the direction of change between two periods.

    Returns:
        str  'decline', 'increase', or 'stable'
    """
    if area_start is None or area_end is None:
        return None

    delta = area_end - area_start
    # Use a 0.5% relative threshold to avoid noise-driven misclassification
    relative_change = abs(delta) / area_start if area_start > 0 else 0

    if relative_change < 0.005:
        return 'stable'
    elif delta < 0:
        return 'decline'
    else:
        return 'increase'


def compute_directional_agreement(products):
    """
    For each interval (T1->T2, T2->T3, T3->T4), count how many products
    agree on decline, increase, or stable.

    Parameters:
        products : dict  {product_name: {period_key: forest_ha or None}}

    Returns:
        dict  {interval_label: {decline: N, increase: N, stable: N,
                                 product_directions: {name: direction}}}
    """
    agreement = {}

    for interval_label, pk_from, pk_to in INTERVALS:
        directions = {}

        for product_name, areas in products.items():
            area_from = areas.get(pk_from)
            area_to = areas.get(pk_to)
            direction = classify_direction(area_from, area_to)
            if direction is not None:
                directions[product_name] = direction

        # Count
        decline_count = sum(1 for d in directions.values() if d == 'decline')
        increase_count = sum(1 for d in directions.values() if d == 'increase')
        stable_count = sum(1 for d in directions.values() if d == 'stable')

        year_from = TARGET_YEARS[pk_from]
        year_to = TARGET_YEARS[pk_to]

        agreement[interval_label] = {
            'from_year': year_from,
            'to_year': year_to,
            'decline': decline_count,
            'increase': increase_count,
            'stable': stable_count,
            'n_products_with_data': len(directions),
            'majority_direction': max(
                ['decline', 'increase', 'stable'],
                key=lambda d: sum(1 for v in directions.values() if v == d)
            ) if directions else None,
            'product_directions': directions,
        }

    return agreement


# ===========================================================================
# Concordance table
# ===========================================================================

def build_concordance_table(products):
    """
    Build a flat concordance table with one row per (product, period).

    Returns:
        list of dicts
    """
    table = []

    for product_name, areas in products.items():
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            area = areas.get(pk)

            table.append({
                'product': product_name,
                'period': pk,
                'year': year,
                'forest_ha': area,
            })

    return table


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    t0 = time.time()
    print('=' * 70)
    print('MULTI-PRODUCT FOREST AREA CROSS-VALIDATION')
    print(f'Date : {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'BBOX : {BBOX}')
    print(f'Years: {[TARGET_YEARS[pk] for pk in PERIOD_KEYS]}')
    print('=' * 70)

    region = STUDY_AREA_BBOX

    # Container for all product results: {product_name: {period_key: forest_ha}}
    all_products = {}
    product_metadata = {}

    # ------------------------------------------------------------------
    # 1. This study (Olofsson estimates)
    # ------------------------------------------------------------------
    print('\n[1/6] Loading this study (Olofsson estimates) ...')

    try:
        this_study = load_this_study()
        all_products['this_study'] = this_study
        product_metadata['this_study'] = {
            'source': 'Olofsson et al. (2014) stratified area estimators',
            'resolution': '30 m (Landsat)',
            'forest_definition': 'Bosque denso + Bosque secundario (classes 1+2)',
            'temporal_coverage': 'All 4 periods',
        }
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            print(f'  T={year} ({pk}): {this_study[pk]:,.0f} ha')
    except Exception as e:
        print(f'  ERROR loading Olofsson estimates: {e}')
        all_products['this_study'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # 2. MapBiomas Colombia (from prior output)
    # ------------------------------------------------------------------
    print('\n[2/6] Loading MapBiomas Colombia (from 18_mapbiomas_validation.py) ...')

    try:
        if not os.path.isfile(MAPBIOMAS_PATH):
            print(f'  WARNING: MapBiomas validation file not found: {MAPBIOMAS_PATH}')
            print('  Run 18_mapbiomas_validation.py first.')
            all_products['mapbiomas'] = {pk: None for pk in PERIOD_KEYS}
        else:
            mapbiomas = load_mapbiomas()
            all_products['mapbiomas'] = mapbiomas
            product_metadata['mapbiomas'] = {
                'source': 'MapBiomas Colombia Collection 1',
                'resolution': '30 m',
                'forest_definition': 'MapBiomas codes 1-6, 49',
                'temporal_coverage': '2013, 2016, 2020, 2022 (proxy for T4)',
                'note': 'Year 2022 used as proxy for T4=2024',
            }
            for pk in PERIOD_KEYS:
                year = TARGET_YEARS[pk]
                val = mapbiomas.get(pk)
                val_str = f'{val:,.0f} ha' if val is not None else 'N/A'
                print(f'  T={year} ({pk}): {val_str}')
    except Exception as e:
        print(f'  ERROR loading MapBiomas: {e}')
        all_products['mapbiomas'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # 3. ESA WorldCover v200
    # ------------------------------------------------------------------
    print('\n[3/6] Computing ESA WorldCover v200 ...')

    try:
        esa_wc = compute_esa_worldcover(region)
        all_products['esa_worldcover'] = esa_wc
        product_metadata['esa_worldcover'] = {
            'source': 'ESA WorldCover v200',
            'gee_asset': ESA_WORLDCOVER_ID,
            'resolution': '10 m',
            'forest_definition': 'Class 10 (Tree cover)',
            'temporal_coverage': '2021 only (used as proxy for T3=2020)',
            'note': 'Single-year product; only T3 value available',
        }
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            val = esa_wc.get(pk)
            val_str = f'{val:,.0f} ha' if val is not None else 'N/A'
            print(f'  T={year} ({pk}): {val_str}')
    except Exception as e:
        print(f'  ERROR with ESA WorldCover: {e}')
        all_products['esa_worldcover'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # 4. Google Dynamic World v1
    # ------------------------------------------------------------------
    print('\n[4/6] Computing Google Dynamic World v1 ...')

    try:
        dw = compute_dynamic_world(region)
        all_products['dynamic_world'] = dw
        product_metadata['dynamic_world'] = {
            'source': 'Google Dynamic World v1',
            'gee_asset': DYNAMIC_WORLD_ID,
            'resolution': '10 m',
            'forest_definition': 'Class 1 (trees), mode composite per period window',
            'temporal_coverage': '2015-present (T1=2013 not available)',
        }
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            val = dw.get(pk)
            val_str = f'{val:,.0f} ha' if val is not None else 'N/A'
            print(f'  T={year} ({pk}): {val_str}')
    except Exception as e:
        print(f'  ERROR with Dynamic World: {e}')
        all_products['dynamic_world'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # 5. MODIS MCD12Q1 v061
    # ------------------------------------------------------------------
    print('\n[5/6] Computing MODIS MCD12Q1 v061 ...')

    try:
        modis = compute_modis_lc(region)
        all_products['modis_lc'] = modis
        product_metadata['modis_lc'] = {
            'source': 'MODIS MCD12Q1 v061 (IGBP classification)',
            'gee_asset': MODIS_LC_ID,
            'resolution': '500 m',
            'forest_definition': 'IGBP classes 1-5 (forests) + 8-9 (woody/savannas)',
            'temporal_coverage': '2013, 2016, 2020, 2023 (proxy for T4)',
            'note': 'Year 2023 used as proxy for T4=2024',
        }
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            val = modis.get(pk)
            val_str = f'{val:,.0f} ha' if val is not None else 'N/A'
            print(f'  T={year} ({pk}): {val_str}')
    except Exception as e:
        print(f'  ERROR with MODIS LC: {e}')
        all_products['modis_lc'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # 6. Hansen Global Forest Change v1.12
    # ------------------------------------------------------------------
    print('\n[6/6] Computing Hansen GFC v1.12 ...')

    try:
        hansen = compute_hansen_gfc(region)
        all_products['hansen_gfc'] = hansen
        product_metadata['hansen_gfc'] = {
            'source': 'Hansen Global Forest Change v1.12',
            'gee_asset': HANSEN_ID,
            'resolution': '30 m',
            'forest_definition': 'treecover2000 > 25% AND no loss before target year',
            'temporal_coverage': '2000-2024 (loss-only; monotonic decline by design)',
            'note': 'Hansen tracks loss only, not regrowth; area is monotonically non-increasing',
        }
        for pk in PERIOD_KEYS:
            year = TARGET_YEARS[pk]
            val = hansen.get(pk)
            val_str = f'{val:,.0f} ha' if val is not None else 'N/A'
            print(f'  T={year} ({pk}): {val_str}')
    except Exception as e:
        print(f'  ERROR with Hansen GFC: {e}')
        all_products['hansen_gfc'] = {pk: None for pk in PERIOD_KEYS}

    # ------------------------------------------------------------------
    # Directional agreement
    # ------------------------------------------------------------------
    print('\n' + '-' * 70)
    print('DIRECTIONAL AGREEMENT ANALYSIS')
    print('-' * 70)

    directional = compute_directional_agreement(all_products)

    for interval_label, info in directional.items():
        yr_from = info['from_year']
        yr_to = info['to_year']
        n_total = info['n_products_with_data']
        majority = info['majority_direction']
        print(f'\n  {interval_label} ({yr_from} -> {yr_to}):')
        print(f'    Decline:  {info["decline"]}/{n_total}')
        print(f'    Increase: {info["increase"]}/{n_total}')
        print(f'    Stable:   {info["stable"]}/{n_total}')
        print(f'    Majority: {majority}')
        for pname, pdir in info['product_directions'].items():
            print(f'      {pname:<18}: {pdir}')

    # ------------------------------------------------------------------
    # Concordance table
    # ------------------------------------------------------------------
    concordance_table = build_concordance_table(all_products)

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    print('\n' + '-' * 70)
    print('CONCORDANCE SUMMARY TABLE')
    print('-' * 70)

    # Header
    product_names = list(all_products.keys())
    header = f'{"Period":<20}'
    for pn in product_names:
        header += f'{pn:>18}'
    print(header)
    print('-' * (20 + 18 * len(product_names)))

    for pk in PERIOD_KEYS:
        year = TARGET_YEARS[pk]
        row = f'T={year} ({pk[:8]:<8})'
        for pn in product_names:
            val = all_products[pn].get(pk)
            if val is not None:
                row += f'{val:>18,.0f}'
            else:
                row += f'{"N/A":>18}'
        print(row)

    # ------------------------------------------------------------------
    # Compute overall directional agreement score
    # ------------------------------------------------------------------
    total_decline_agree = 0
    total_intervals_with_data = 0

    for interval_label, info in directional.items():
        if info['n_products_with_data'] >= 2:
            total_intervals_with_data += 1
            majority = info['majority_direction']
            # Count fraction of products agreeing with majority
            majority_count = max(info['decline'], info['increase'], info['stable'])
            if majority_count > info['n_products_with_data'] / 2:
                total_decline_agree += 1

    overall_agreement_score = (
        total_decline_agree / total_intervals_with_data * 100
        if total_intervals_with_data > 0 else None
    )

    # ------------------------------------------------------------------
    # Write output JSON
    # ------------------------------------------------------------------
    print('\n' + '-' * 70)
    print('WRITING OUTPUT')
    print('-' * 70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    output = {
        'metadata': {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'purpose': (
                'Multi-product cross-validation of study forest area estimates '
                'against 5 independent land-cover products'
            ),
            'study_bbox': BBOX,
            'target_years': {pk: TARGET_YEARS[pk] for pk in PERIOD_KEYS},
            'period_labels': {pk: PERIODS[pk]['label'] for pk in PERIOD_KEYS},
            'n_products': len(all_products),
            'product_list': list(all_products.keys()),
            'product_metadata': product_metadata,
            'olofsson_source': OLOFSSON_PATH,
            'mapbiomas_source': MAPBIOMAS_PATH,
        },
        'products': {
            pname: {
                pk: areas.get(pk) for pk in PERIOD_KEYS
            }
            for pname, areas in all_products.items()
        },
        'directional_agreement': directional,
        'concordance_table': concordance_table,
        'summary': {
            'overall_majority_agreement_pct': (
                round(overall_agreement_score, 1)
                if overall_agreement_score is not None else None
            ),
            'intervals_with_majority_consensus': total_decline_agree,
            'total_intervals_assessed': total_intervals_with_data,
            'note': (
                'Directional agreement counts how many of the 3 intervals '
                'have a clear majority (>50%) of products agreeing on the '
                'direction of change.'
            ),
        },
    }

    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f'  Written: {OUTPUT_PATH}')

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)

    n_products_ok = sum(
        1 for pname, areas in all_products.items()
        if any(v is not None for v in areas.values())
    )
    print(f'  Products with data        : {n_products_ok}/{len(all_products)}')

    for interval_label, info in directional.items():
        yr_from = info['from_year']
        yr_to = info['to_year']
        majority = info['majority_direction']
        n_agree = max(info['decline'], info['increase'], info['stable'])
        n_total = info['n_products_with_data']
        print(f'  {interval_label} ({yr_from}->{yr_to})  : '
              f'{majority} ({n_agree}/{n_total} products)')

    if overall_agreement_score is not None:
        print(f'  Majority consensus        : '
              f'{total_decline_agree}/{total_intervals_with_data} intervals '
              f'({overall_agreement_score:.1f}%)')

    print(f'  Output file               : {OUTPUT_PATH}')
    print(f'  Elapsed                   : {elapsed:.1f} s')
    print('=' * 70)


if __name__ == '__main__':
    main()
