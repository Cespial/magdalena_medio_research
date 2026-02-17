#!/usr/bin/env python3
"""
08c_water_yield_terraclimate.py
================================
Phase 3.2 enhancement: Water yield estimation using TerraClimate + ERA5-Land
cross-validation to replace non-functional MODIS MOD16A2.

Problem:
  The original 08_ecosystem_services.py used MODIS MOD16A2 for actual
  evapotranspiration (AET), but returned null/zero for T1-T3 (2013, 2016,
  2020) due to coverage gaps in the Magdalena Medio region.

Solution:
  Uses two independent, continuous global gridded datasets:
    1. TerraClimate (IDAHO_EPSCOR/TERRACLIMATE)
       - Monthly, ~4km, 1958-present
       - Bands: 'aet' (actual ET, scale 0.1 mm), 'pet' (potential ET, scale 0.1 mm)
    2. ERA5-Land (ECMWF/ERA5_LAND/MONTHLY_AGGR)
       - Monthly, ~11km
       - Band: 'total_evaporation_sum' (m, negative = actual ET)

  Precipitation from CHIRPS daily (UCSB-CHG/CHIRPS/DAILY).

  For each period:
    - Annual P from CHIRPS
    - AET from TerraClimate (sum monthly 'aet' * 0.1)
    - AET from ERA5-Land (sum monthly 'total_evaporation_sum' * -1000)
    - Water Yield = P - AET (for each source)
    - LULC-weighted water yield using Kc coefficients (FAO 56)
    - Baseflow recharge = P * infiltration_coeff per LULC class

  Cross-validation: comparison of TerraClimate vs ERA5-Land AET estimates
  provides a measure of inter-source agreement.

Outputs:
  - outputs/phase3_stats/water_yield_terraclimate.json
  - Updates outputs/phase3_stats/ecosystem_services_results.json
"""

import ee
import json
import math
import os
import sys
from datetime import datetime

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from gee_config import PERIODS, LULC_CLASSES, STUDY_AREA_BBOX
from scripts.utils import get_study_area

OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'outputs', 'phase3_stats')
OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'water_yield_terraclimate.json')
ES_RESULTS_PATH = os.path.join(OUTPUT_DIR, 'ecosystem_services_results.json')
OLOFSSON_PATH = os.path.join(OUTPUT_DIR, 'olofsson_area_estimates.json')

# ---------------------------------------------------------------------------
# LULC coefficients
# ---------------------------------------------------------------------------

# Kc crop coefficients (FAO 56) by LULC class
KC_VALUES = {
    1: 1.00,   # Bosque denso - alta ET
    2: 0.85,   # Bosque secundario
    3: 0.60,   # Pasturas
    4: 0.70,   # Cultivos
    5: 1.20,   # Agua - evaporacion abierta
    6: 0.30,   # Urbano - impermeabilizado
    7: 0.15,   # Suelo desnudo
}

# Infiltration / baseflow recharge coefficients
RECHARGE_VALUES = {
    1: 0.35,   # Bosque denso - alta infiltracion
    2: 0.30,   # Bosque secundario
    3: 0.15,   # Pasturas - escorrentia
    4: 0.18,   # Cultivos
    5: 0.00,   # Agua - no aplica
    6: 0.05,   # Urbano - casi impermeable
    7: 0.10,   # Suelo desnudo
}

CLASS_NAMES = {
    1: 'Bosque denso',
    2: 'Bosque secundario',
    3: 'Pasturas',
    4: 'Cultivos',
    5: 'Agua',
    6: 'Urbano',
    7: 'Suelo desnudo',
}


# ---------------------------------------------------------------------------
# GEE computation functions
# ---------------------------------------------------------------------------

def compute_annual_precipitation(year, region):
    """
    Compute annual precipitation (mm) from CHIRPS daily.

    Args:
        year: int, calendar year
        region: ee.Geometry

    Returns:
        ee.Image with annual precipitation in mm
    """
    chirps = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filterDate(f'{year}-01-01', f'{year}-12-31')
              .filterBounds(region)
              .sum()
              .clip(region)
              .rename('precipitation'))
    return chirps


def compute_aet_terraclimate(year, region):
    """
    Compute annual actual ET from TerraClimate.
    Band 'aet' is in 0.1 mm units; multiply by 0.1 to get mm.

    Args:
        year: int, calendar year
        region: ee.Geometry

    Returns:
        ee.Image with annual AET in mm
    """
    tc = (ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
          .filterDate(f'{year}-01-01', f'{year}-12-31')
          .filterBounds(region)
          .select('aet')
          .sum()
          .multiply(0.1)
          .clip(region)
          .rename('aet_terraclimate'))
    return tc


def compute_pet_terraclimate(year, region):
    """
    Compute annual potential ET from TerraClimate.
    Band 'pet' is in 0.1 mm units.

    Args:
        year: int, calendar year
        region: ee.Geometry

    Returns:
        ee.Image with annual PET in mm
    """
    tc = (ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
          .filterDate(f'{year}-01-01', f'{year}-12-31')
          .filterBounds(region)
          .select('pet')
          .sum()
          .multiply(0.1)
          .clip(region)
          .rename('pet_terraclimate'))
    return tc


def compute_aet_era5(year, region):
    """
    Compute annual actual ET from ERA5-Land monthly aggregates.
    Band 'total_evaporation_sum' is in meters and is negative
    (convention: negative = flux away from surface).
    Multiply by -1000 to convert to positive mm.

    Args:
        year: int, calendar year
        region: ee.Geometry

    Returns:
        ee.Image with annual AET in mm
    """
    era5 = (ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
            .filterDate(f'{year}-01-01', f'{year}-12-31')
            .filterBounds(region)
            .select('total_evaporation_sum')
            .sum()
            .multiply(-1000)
            .clip(region)
            .rename('aet_era5'))
    return era5


def compute_mean_over_region(image, region, scale=4000):
    """
    Compute spatial mean of an image over a region.

    Args:
        image: ee.Image (single band)
        region: ee.Geometry
        scale: int, resolution in meters (default 4km for TerraClimate)

    Returns:
        float or None
    """
    stats = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=scale,
        maxPixels=1e13
    )
    result = stats.getInfo()
    # Return the first (and only) band value
    if result:
        vals = list(result.values())
        if vals and vals[0] is not None:
            return round(vals[0], 2)
    return None


def compute_lulc_weighted_water_yield(precip_img, aet_img, region, olofsson_areas):
    """
    Compute LULC-weighted water yield using Kc coefficients and
    Olofsson-adjusted area proportions.

    Water yield per class = (P - AET * Kc_i)
    Weighted total = sum( WY_i * area_fraction_i )

    Args:
        precip_img: ee.Image of annual precipitation (mm)
        aet_img: ee.Image of annual AET (mm)
        region: ee.Geometry
        olofsson_areas: dict from olofsson_area_estimates.json per_class

    Returns:
        float, LULC-weighted water yield (mm)
    """
    # Get region-mean P and AET
    p_mean = compute_mean_over_region(precip_img, region, scale=5000)
    aet_mean = compute_mean_over_region(aet_img, region, scale=5000)

    if p_mean is None or aet_mean is None:
        return None

    # Compute total area for weighting
    total_area = 0.0
    class_areas = {}
    for cls_id in range(1, 8):
        cls_str = str(cls_id)
        if cls_str in olofsson_areas:
            area_ha = olofsson_areas[cls_str].get('adjusted_area_ha', 0.0)
        else:
            area_ha = 0.0
        class_areas[cls_id] = area_ha
        total_area += area_ha

    if total_area == 0:
        return None

    # Compute weighted water yield
    weighted_wy = 0.0
    for cls_id in range(1, 8):
        kc = KC_VALUES[cls_id]
        area_frac = class_areas[cls_id] / total_area
        wy_class = p_mean - (aet_mean * kc)
        weighted_wy += wy_class * area_frac

    return round(weighted_wy, 2)


def compute_baseflow_recharge(precip_img, region, olofsson_areas):
    """
    Compute LULC-weighted baseflow recharge.
    Baseflow per class = P * infiltration_coeff_i
    Weighted total = sum( BF_i * area_fraction_i )

    Args:
        precip_img: ee.Image of annual precipitation (mm)
        region: ee.Geometry
        olofsson_areas: dict from olofsson_area_estimates.json per_class

    Returns:
        float, LULC-weighted baseflow (mm)
    """
    p_mean = compute_mean_over_region(precip_img, region, scale=5000)
    if p_mean is None:
        return None

    total_area = 0.0
    class_areas = {}
    for cls_id in range(1, 8):
        cls_str = str(cls_id)
        if cls_str in olofsson_areas:
            area_ha = olofsson_areas[cls_str].get('adjusted_area_ha', 0.0)
        else:
            area_ha = 0.0
        class_areas[cls_id] = area_ha
        total_area += area_ha

    if total_area == 0:
        return None

    weighted_bf = 0.0
    for cls_id in range(1, 8):
        rc = RECHARGE_VALUES[cls_id]
        area_frac = class_areas[cls_id] / total_area
        bf_class = p_mean * rc
        weighted_bf += bf_class * area_frac

    return round(weighted_bf, 2)


def compute_aridity_index(pet_img, precip_img, region):
    """
    Compute aridity index = PET / P.
    Values > 1 indicate arid conditions; < 1 indicate humid conditions.

    Args:
        pet_img: ee.Image of annual PET (mm)
        precip_img: ee.Image of annual P (mm)
        region: ee.Geometry

    Returns:
        float or None
    """
    pet_mean = compute_mean_over_region(pet_img, region, scale=5000)
    p_mean = compute_mean_over_region(precip_img, region, scale=5000)
    if pet_mean is None or p_mean is None or p_mean == 0:
        return None
    return round(pet_mean / p_mean, 4)


# ---------------------------------------------------------------------------
# Period processing
# ---------------------------------------------------------------------------

def process_period(period_key, period_info, region, olofsson_per_class):
    """
    Process a single period: compute all water yield components.

    Args:
        period_key: str (e.g. 'pre_acuerdo')
        period_info: dict from PERIODS
        region: ee.Geometry
        olofsson_per_class: dict of per_class areas

    Returns:
        dict with all water yield metrics for this period
    """
    year = period_info['map_year']
    label = period_info['label']

    print(f"\n  Processing {label} (year={year})...")

    result = {
        'year': year,
        'label': label,
    }

    # --- 1. Annual precipitation (CHIRPS) ---
    print(f"    [1/6] CHIRPS precipitation {year}...")
    try:
        precip_img = compute_annual_precipitation(year, region)
        p_mm = compute_mean_over_region(precip_img, region, scale=5000)
        result['precipitation_mm'] = p_mm
        print(f"          P = {p_mm} mm/yr")
    except Exception as e:
        print(f"          ERROR computing precipitation: {e}")
        result['precipitation_mm'] = None
        precip_img = None

    # --- 2. AET from TerraClimate ---
    print(f"    [2/6] TerraClimate AET {year}...")
    try:
        aet_tc_img = compute_aet_terraclimate(year, region)
        aet_tc_mm = compute_mean_over_region(aet_tc_img, region, scale=4000)
        result['aet_terraclimate_mm'] = aet_tc_mm
        print(f"          AET_tc = {aet_tc_mm} mm/yr")
    except Exception as e:
        print(f"          ERROR computing TerraClimate AET: {e}")
        result['aet_terraclimate_mm'] = None
        aet_tc_img = None

    # --- 3. AET from ERA5-Land ---
    print(f"    [3/6] ERA5-Land AET {year}...")
    try:
        aet_era5_img = compute_aet_era5(year, region)
        aet_era5_mm = compute_mean_over_region(aet_era5_img, region, scale=11000)
        result['aet_era5_mm'] = aet_era5_mm
        print(f"          AET_era5 = {aet_era5_mm} mm/yr")
    except Exception as e:
        print(f"          ERROR computing ERA5-Land AET: {e}")
        result['aet_era5_mm'] = None
        aet_era5_img = None

    # --- 4. Water Yield = P - AET ---
    print(f"    [4/6] Water yield (P - AET)...")
    try:
        if result['precipitation_mm'] is not None and result['aet_terraclimate_mm'] is not None:
            wy_tc = round(result['precipitation_mm'] - result['aet_terraclimate_mm'], 2)
            result['water_yield_terraclimate_mm'] = wy_tc
            print(f"          WY_tc = {wy_tc} mm/yr")
        else:
            result['water_yield_terraclimate_mm'] = None
            print(f"          WY_tc = None (missing inputs)")

        if result['precipitation_mm'] is not None and result['aet_era5_mm'] is not None:
            wy_era5 = round(result['precipitation_mm'] - result['aet_era5_mm'], 2)
            result['water_yield_era5_mm'] = wy_era5
            print(f"          WY_era5 = {wy_era5} mm/yr")
        else:
            result['water_yield_era5_mm'] = None
            print(f"          WY_era5 = None (missing inputs)")
    except Exception as e:
        print(f"          ERROR computing water yield: {e}")
        result['water_yield_terraclimate_mm'] = None
        result['water_yield_era5_mm'] = None

    # --- 5. PET and aridity index ---
    print(f"    [5/6] PET and aridity index...")
    try:
        pet_img = compute_pet_terraclimate(year, region)
        pet_mm = compute_mean_over_region(pet_img, region, scale=4000)
        result['pet_terraclimate_mm'] = pet_mm
        print(f"          PET = {pet_mm} mm/yr")

        if precip_img is not None:
            ai = compute_aridity_index(pet_img, precip_img, region)
            result['aridity_index'] = ai
            print(f"          AI (PET/P) = {ai}")
        else:
            result['aridity_index'] = None
    except Exception as e:
        print(f"          ERROR computing PET/AI: {e}")
        result['pet_terraclimate_mm'] = None
        result['aridity_index'] = None

    # --- 6. LULC-weighted water yield and baseflow ---
    print(f"    [6/6] LULC-weighted WY and baseflow recharge...")
    try:
        if precip_img is not None and aet_tc_img is not None:
            lulc_wy = compute_lulc_weighted_water_yield(
                precip_img, aet_tc_img, region, olofsson_per_class
            )
            result['lulc_weighted_wy_mm'] = lulc_wy
            print(f"          LULC-weighted WY = {lulc_wy} mm/yr")
        else:
            result['lulc_weighted_wy_mm'] = None
            print(f"          LULC-weighted WY = None (missing inputs)")

        if precip_img is not None:
            bf = compute_baseflow_recharge(precip_img, region, olofsson_per_class)
            result['baseflow_mm'] = bf
            print(f"          Baseflow recharge = {bf} mm/yr")
        else:
            result['baseflow_mm'] = None
            print(f"          Baseflow recharge = None (missing inputs)")
    except Exception as e:
        print(f"          ERROR computing LULC-weighted metrics: {e}")
        result['lulc_weighted_wy_mm'] = None
        result['baseflow_mm'] = None

    # --- Cross-validation: AET agreement ---
    try:
        if (result['aet_terraclimate_mm'] is not None
                and result['aet_era5_mm'] is not None):
            aet_diff = round(
                result['aet_terraclimate_mm'] - result['aet_era5_mm'], 2
            )
            aet_mean_both = round(
                (result['aet_terraclimate_mm'] + result['aet_era5_mm']) / 2, 2
            )
            if aet_mean_both != 0:
                aet_rel_diff_pct = round(
                    abs(aet_diff) / aet_mean_both * 100, 2
                )
            else:
                aet_rel_diff_pct = None
            result['cross_validation'] = {
                'aet_difference_mm': aet_diff,
                'aet_mean_mm': aet_mean_both,
                'aet_relative_diff_pct': aet_rel_diff_pct,
            }
            print(f"          Cross-val: AET diff = {aet_diff} mm "
                  f"({aet_rel_diff_pct}% relative)")
        else:
            result['cross_validation'] = None
    except Exception as e:
        print(f"          ERROR in cross-validation: {e}")
        result['cross_validation'] = None

    return result


# ---------------------------------------------------------------------------
# Temporal trend computation
# ---------------------------------------------------------------------------

def compute_temporal_trend(period_results):
    """
    Compute temporal trends across the 4 periods.

    Args:
        period_results: dict of period_key -> result dict

    Returns:
        dict with trend metrics
    """
    ordered_keys = ['pre_acuerdo', 'transicion', 'post_acuerdo_1', 'post_acuerdo_2']
    years = []
    wy_tc_values = []
    wy_era5_values = []
    precip_values = []
    baseflow_values = []

    for pk in ordered_keys:
        if pk not in period_results:
            continue
        r = period_results[pk]
        years.append(r['year'])
        precip_values.append(r.get('precipitation_mm'))
        wy_tc_values.append(r.get('water_yield_terraclimate_mm'))
        wy_era5_values.append(r.get('water_yield_era5_mm'))
        baseflow_values.append(r.get('baseflow_mm'))

    trend = {}

    # Precipitation trend
    valid_p = [(y, v) for y, v in zip(years, precip_values) if v is not None]
    if len(valid_p) >= 2:
        first_p = valid_p[0][1]
        last_p = valid_p[-1][1]
        trend['precipitation_change_mm'] = round(last_p - first_p, 2)
        if first_p != 0:
            trend['precipitation_change_pct'] = round(
                (last_p - first_p) / first_p * 100, 2
            )

    # Water yield TerraClimate trend
    valid_tc = [(y, v) for y, v in zip(years, wy_tc_values) if v is not None]
    if len(valid_tc) >= 2:
        first_tc = valid_tc[0][1]
        last_tc = valid_tc[-1][1]
        trend['wy_terraclimate_change_mm'] = round(last_tc - first_tc, 2)
        if first_tc != 0:
            trend['wy_terraclimate_change_pct'] = round(
                (last_tc - first_tc) / first_tc * 100, 2
            )

    # Water yield ERA5 trend
    valid_era5 = [(y, v) for y, v in zip(years, wy_era5_values) if v is not None]
    if len(valid_era5) >= 2:
        first_era5 = valid_era5[0][1]
        last_era5 = valid_era5[-1][1]
        trend['wy_era5_change_mm'] = round(last_era5 - first_era5, 2)
        if first_era5 != 0:
            trend['wy_era5_change_pct'] = round(
                (last_era5 - first_era5) / first_era5 * 100, 2
            )

    # Baseflow trend
    valid_bf = [(y, v) for y, v in zip(years, baseflow_values) if v is not None]
    if len(valid_bf) >= 2:
        first_bf = valid_bf[0][1]
        last_bf = valid_bf[-1][1]
        trend['baseflow_change_mm'] = round(last_bf - first_bf, 2)
        if first_bf != 0:
            trend['baseflow_change_pct'] = round(
                (last_bf - first_bf) / first_bf * 100, 2
            )

    # Period-to-period changes (TerraClimate WY)
    changes = []
    for i in range(len(valid_tc) - 1):
        y1, v1 = valid_tc[i]
        y2, v2 = valid_tc[i + 1]
        changes.append({
            'from_year': y1,
            'to_year': y2,
            'change_mm': round(v2 - v1, 2),
        })
    if changes:
        trend['period_changes_wy_tc'] = changes

    return trend


# ---------------------------------------------------------------------------
# Merge into ecosystem_services_results.json
# ---------------------------------------------------------------------------

def merge_into_ecosystem_services(period_results):
    """
    Read the existing ecosystem_services_results.json and merge
    updated water yield data into it.

    Args:
        period_results: dict of period_key -> result dict
    """
    if not os.path.exists(ES_RESULTS_PATH):
        print(f"\n  WARNING: {ES_RESULTS_PATH} not found; skipping merge.")
        return

    try:
        with open(ES_RESULTS_PATH, 'r') as f:
            existing = json.load(f)
    except Exception as e:
        print(f"\n  ERROR reading existing ES results: {e}")
        return

    updated = False
    for period_key, wy_result in period_results.items():
        if period_key in existing:
            # Update water yield fields
            wy_tc = wy_result.get('water_yield_terraclimate_mm')
            if wy_tc is not None:
                existing[period_key]['water_yield_mm'] = wy_tc
                updated = True
            bf = wy_result.get('baseflow_mm')
            if bf is not None:
                existing[period_key]['baseflow_mm'] = bf
                updated = True

    if updated:
        # Update metadata
        if '_metadata' in existing:
            existing['_metadata']['water_yield_source'] = (
                'TerraClimate AET (IDAHO_EPSCOR/TERRACLIMATE), '
                'cross-validated with ERA5-Land'
            )
            existing['_metadata']['water_yield_updated'] = (
                datetime.now().strftime('%Y-%m-%d %H:%M')
            )

        with open(ES_RESULTS_PATH, 'w') as f:
            json.dump(existing, f, indent=2, ensure_ascii=False)
        print(f"\n  Updated: {ES_RESULTS_PATH}")
    else:
        print(f"\n  No water yield data to merge into ES results.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("08c  WATER YIELD: TerraClimate + ERA5-Land Cross-Validation")
    print(f"     {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 65)

    print("\nReplaces non-functional MODIS MOD16A2 water yield (null for T1-T3)")
    print("Data sources:")
    print("  - Precipitation: CHIRPS Daily (UCSB-CHG/CHIRPS/DAILY)")
    print("  - AET source 1:  TerraClimate (IDAHO_EPSCOR/TERRACLIMATE)")
    print("  - AET source 2:  ERA5-Land (ECMWF/ERA5_LAND/MONTHLY_AGGR)")
    print("  - LULC areas:    Olofsson et al. (2014) stratified estimates")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    region = get_study_area()

    # ------------------------------------------------------------------
    # Load Olofsson area estimates
    # ------------------------------------------------------------------
    print(f"\nLoading Olofsson area estimates...")
    try:
        with open(OLOFSSON_PATH, 'r') as f:
            olofsson = json.load(f)
        olofsson_periods = olofsson['periods']
        print(f"  Loaded {len(olofsson_periods)} periods from Olofsson estimates")
    except Exception as e:
        print(f"  ERROR loading Olofsson estimates: {e}")
        print("  Will proceed without LULC weighting.")
        olofsson_periods = {}

    # ------------------------------------------------------------------
    # Print Kc and recharge coefficient tables
    # ------------------------------------------------------------------
    print(f"\nKc coefficients (FAO 56):")
    print(f"  {'Class':<22s} {'Kc':>6s} {'Recharge':>10s}")
    for cls_id in range(1, 8):
        print(f"  {CLASS_NAMES[cls_id]:<22s} {KC_VALUES[cls_id]:6.2f} "
              f"{RECHARGE_VALUES[cls_id]:10.2f}")

    # ------------------------------------------------------------------
    # Process each period
    # ------------------------------------------------------------------
    period_results = {}
    period_keys = ['pre_acuerdo', 'transicion', 'post_acuerdo_1', 'post_acuerdo_2']

    for pk in period_keys:
        if pk not in PERIODS:
            print(f"\n  WARNING: Period '{pk}' not found in PERIODS config; skipping.")
            continue

        period_info = PERIODS[pk]

        # Get Olofsson per_class for this period
        if pk in olofsson_periods:
            olofsson_per_class = olofsson_periods[pk].get('per_class', {})
        else:
            olofsson_per_class = {}

        try:
            result = process_period(pk, period_info, region, olofsson_per_class)
            period_results[pk] = result
        except Exception as e:
            print(f"\n  ERROR processing period {pk}: {e}")
            period_results[pk] = {
                'year': period_info['map_year'],
                'label': period_info['label'],
                'error': str(e),
            }

    # ------------------------------------------------------------------
    # Compute temporal trends
    # ------------------------------------------------------------------
    print("\n" + "-" * 65)
    print("  TEMPORAL TRENDS")
    print("-" * 65)

    try:
        temporal_trend = compute_temporal_trend(period_results)
        for key, val in temporal_trend.items():
            if not isinstance(val, list):
                print(f"  {key}: {val}")
        if 'period_changes_wy_tc' in temporal_trend:
            for ch in temporal_trend['period_changes_wy_tc']:
                print(f"  WY change {ch['from_year']}-{ch['to_year']}: "
                      f"{ch['change_mm']:+.2f} mm")
    except Exception as e:
        print(f"  ERROR computing temporal trends: {e}")
        temporal_trend = {'error': str(e)}

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    print("\n" + "-" * 65)
    print("  SUMMARY TABLE")
    print("-" * 65)
    print(f"  {'Period':<16s} {'Year':>5s} {'P(mm)':>8s} {'AET_tc':>8s} "
          f"{'AET_e5':>8s} {'WY_tc':>8s} {'WY_e5':>8s} {'BF':>8s}")
    print(f"  {'-'*16:<16s} {'-----':>5s} {'--------':>8s} {'--------':>8s} "
          f"{'--------':>8s} {'--------':>8s} {'--------':>8s} {'--------':>8s}")

    for pk in period_keys:
        if pk not in period_results:
            continue
        r = period_results[pk]

        def fmt(v):
            return f"{v:8.1f}" if v is not None else "    null"

        print(f"  {pk:<16s} {r['year']:5d} {fmt(r.get('precipitation_mm'))} "
              f"{fmt(r.get('aet_terraclimate_mm'))} {fmt(r.get('aet_era5_mm'))} "
              f"{fmt(r.get('water_yield_terraclimate_mm'))} "
              f"{fmt(r.get('water_yield_era5_mm'))} {fmt(r.get('baseflow_mm'))}")

    # ------------------------------------------------------------------
    # Build output JSON
    # ------------------------------------------------------------------
    output = {
        'periods': period_results,
        'temporal_trend': temporal_trend,
        '_metadata': {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'script': '08c_water_yield_terraclimate.py',
            'description': (
                'Water yield estimation using TerraClimate + ERA5-Land '
                'cross-validation. Replaces non-functional MODIS MOD16A2 '
                'approach (null for T1-T3).'
            ),
            'precipitation_source': 'CHIRPS Daily (UCSB-CHG/CHIRPS/DAILY)',
            'aet_source_primary': 'TerraClimate (IDAHO_EPSCOR/TERRACLIMATE), ~4km monthly',
            'aet_source_secondary': 'ERA5-Land (ECMWF/ERA5_LAND/MONTHLY_AGGR), ~11km monthly',
            'lulc_areas': 'Olofsson et al. (2014) stratified area estimates',
            'kc_reference': 'FAO Irrigation and Drainage Paper 56',
            'kc_coefficients': KC_VALUES,
            'recharge_coefficients': RECHARGE_VALUES,
            'water_yield_formula': 'WY = P - AET (simple water balance)',
            'lulc_weighted_formula': 'WY_weighted = sum(area_frac_i * (P - AET * Kc_i))',
            'baseflow_formula': 'BF = sum(area_frac_i * P * recharge_coeff_i)',
            'cross_validation_note': (
                'AET from TerraClimate and ERA5-Land are compared; '
                'relative difference provides inter-source agreement metric.'
            ),
            'units': {
                'precipitation_mm': 'mm/year',
                'aet_mm': 'mm/year',
                'water_yield_mm': 'mm/year',
                'baseflow_mm': 'mm/year',
                'aridity_index': 'dimensionless (PET/P)',
            },
        },
    }

    # Write water yield results
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    print(f"\nResults written to: {OUTPUT_PATH}")

    # ------------------------------------------------------------------
    # Merge into ecosystem_services_results.json
    # ------------------------------------------------------------------
    print("\nMerging water yield into ecosystem_services_results.json...")
    merge_into_ecosystem_services(period_results)

    print("\nDone.")
    print("Next step: 09_climate_analysis.py")

    return output


if __name__ == '__main__':
    results = main()
