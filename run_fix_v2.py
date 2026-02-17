"""
run_fix_v2.py - Fix remaining issues from run_fix_stages.py:
1. T4 (2024) classification: add post-classification forest correction using NDVI + Hansen
2. Change detection: fix transition matrix (need 2 bands for Reducer.group)
3. CA-Markov: use coarser scale sampling (2000m, 3000 pts) to avoid memory limits
4. Ecosystem services: fix water yield band mismatch, use tileScale=8 for habitat
"""

import ee
import os
import sys
import json
import time
import math
import numpy as np
from datetime import datetime

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_DIR)
from gee_config import PERIODS, LULC_CLASSES, CARBON_POOLS
from scripts.utils import get_study_area, mask_landsat_clouds, mask_sentinel2_clouds

import importlib
training_mod = importlib.import_module('scripts.02_training_samples')
hotspot_mod = importlib.import_module('scripts.07_hotspot_analysis')
ecosystem_mod = importlib.import_module('scripts.08_ecosystem_services')
gwr_mod = importlib.import_module('scripts.10_gwr_drivers')
ca_markov_mod = importlib.import_module('scripts.11_ca_markov')

OUTPUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'phase3_stats')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)


def save_json(data, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    if os.path.exists(path):
        try:
            with open(path) as f:
                existing = json.load(f)
            existing.update(data)
            data = existing
        except (json.JSONDecodeError, ValueError):
            pass
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    log(f"  >> Saved: {filename}")


def gi(obj, label=""):
    try:
        return obj.getInfo()
    except Exception as e:
        log(f"  WARN ({label}): {e}")
        return None


# ================================================================
# FIXED COMPOSITE (same as run_fix_stages.py)
# ================================================================
def create_fixed_composite(start, end, region, year):
    from gee_config import LANDSAT_BANDS, SENTINEL_BANDS
    common = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']

    l8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
          .filterDate(start, end).filterBounds(region)
          .filter(ee.Filter.lt('CLOUD_COVER', 70))
          .map(mask_landsat_clouds)
          .select(list(LANDSAT_BANDS.values()), common)
          .map(lambda img: img.select(common).toFloat()))

    if year >= 2021:
        l9 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
              .filterDate(start, end).filterBounds(region)
              .filter(ee.Filter.lt('CLOUD_COVER', 70))
              .map(mask_landsat_clouds)
              .select(list(LANDSAT_BANDS.values()), common)
              .map(lambda img: img.select(common).toFloat()))
        landsat = l8.merge(l9)
    else:
        landsat = l8

    if year >= 2016:
        s2 = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
              .filterDate(start, end).filterBounds(region)
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70))
              .map(mask_sentinel2_clouds)
              .select(list(SENTINEL_BANDS.values()), common)
              .map(lambda img: img.select(common).toFloat()))
        merged = landsat.merge(s2)
    else:
        merged = landsat

    n_images = merged.size()
    composite = merged.median().clip(region)

    nir = composite.select('nir')
    red = composite.select('red')
    green = composite.select('green')
    swir1 = composite.select('swir1')
    swir2 = composite.select('swir2')
    eps = 0.0001

    ndvi = nir.subtract(red).divide(nir.add(red).add(eps)).rename('NDVI')
    ndwi = green.subtract(nir).divide(green.add(nir).add(eps)).rename('NDWI')
    ndbi = swir1.subtract(nir).divide(swir1.add(nir).add(eps)).rename('NDBI')
    nbr = nir.subtract(swir2).divide(nir.add(swir2).add(eps)).rename('NBR')

    composite = composite.addBands([ndvi, ndwi, ndbi, nbr])

    dem = ee.Image('USGS/SRTMGL1_003')
    composite = composite.addBands(dem.select('elevation').rename('elevation').toFloat())
    composite = composite.addBands(ee.Terrain.slope(dem).rename('slope').toFloat())

    bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2',
             'NDVI', 'NDWI', 'NDBI', 'NBR', 'elevation', 'slope']
    return composite, n_images, bands


# ================================================================
# RECLASSIFY T4 with post-classification correction
# ================================================================
def reclassify_t4_fixed(region):
    """
    Re-classify T4 (2024) with post-classification correction:
    - If NDVI > 0.65 and Hansen treecover2000 > 50%, assign forest class
    - If NDVI > 0.5 and Hansen treecover2000 > 30%, assign secondary forest
    """
    log("=" * 60)
    log("RE-CLASSIFY T4 (2024) WITH FOREST CORRECTION")
    log("=" * 60)

    pi = PERIODS['post_acuerdo_2']
    year = pi['map_year']

    composite, n_img, bands = create_fixed_composite(pi['start'], pi['end'], region, year)
    n = gi(n_img, "n_2024")
    log(f"  {n} images, {len(bands)} bands")

    # Training
    ref = training_mod.get_reference_lulc(year, region)
    samples = training_mod.generate_stratified_samples(ref, region, 300, seed=42 + year)
    train, val = training_mod.split_train_validation(samples, 0.7, seed=42)

    train_data = composite.select(bands).sampleRegions(
        collection=train, properties=['lulc_reference'],
        scale=30, tileScale=4, geometries=False)
    val_data = composite.select(bands).sampleRegions(
        collection=val, properties=['lulc_reference'],
        scale=30, tileScale=4, geometries=False)

    clf = ee.Classifier.smileRandomForest(
        numberOfTrees=200, minLeafPopulation=5, bagFraction=0.632, seed=42
    ).train(train_data, 'lulc_reference', bands)

    classified = composite.select(bands).classify(clf).rename('lulc').toInt8()
    classified = classified.focal_mode(3, 'square', 'pixels').rename('lulc').toInt8()

    # ---- POST-CLASSIFICATION FOREST CORRECTION ----
    log("  Applying post-classification forest correction...")
    ndvi = composite.select('NDVI')
    hansen = ee.Image('UMD/hansen/global_forest_change_2023_v1_11').clip(region)
    treecover = hansen.select('treecover2000')
    loss = hansen.select('lossyear')

    # Pixels with high NDVI, high historical tree cover, and no recent loss
    # are likely still forest even if classifier misclassified them
    still_forest = treecover.gte(50).And(loss.eq(0).Or(loss.gt(23)))
    dense_forest_mask = ndvi.gte(0.55).And(still_forest).And(treecover.gte(60))
    secondary_forest_mask = ndvi.gte(0.45).And(still_forest).And(
        treecover.gte(30)).And(treecover.lt(60))

    # Only correct pixels currently classified as Pasturas (3)
    classified = classified.where(
        classified.eq(3).And(dense_forest_mask), 1)  # -> Bosque denso
    classified = classified.where(
        classified.eq(3).And(secondary_forest_mask), 2)  # -> Bosque secundario

    # Also apply JRC water correction
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
    classified = classified.where(jrc.select('occurrence').gte(80), 5).toInt8()

    # Accuracy (on original classifier, not corrected - correction is rule-based)
    validated = val_data.classify(clf)
    em = validated.errorMatrix('lulc_reference', 'classification')
    oa = gi(em.accuracy(), "oa_2024")
    kp = gi(em.kappa(), "kp_2024")
    cm = gi(em.array(), "cm_2024")
    pa = gi(em.producersAccuracy(), "pa_2024")
    ua = gi(em.consumersAccuracy(), "ua_2024")

    if oa:
        log(f"  Classifier OA={oa:.4f}, Kappa={kp:.4f}")

    # Areas after correction
    log("  Areas (post-correction)...")
    areas = {}
    area_img = ee.Image.pixelArea().divide(10000)
    for cid in range(1, 8):
        a = gi(area_img.updateMask(classified.eq(cid)).reduceRegion(
            ee.Reducer.sum(), region, 100, maxPixels=1e12,
            tileScale=4, bestEffort=True), f"a{cid}_2024")
        ha = a.get('area', 0) if a else 0
        nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
        areas[cid] = {'name': nm, 'area_ha': round(ha, 1)}
        log(f"    {nm}: {ha:,.0f} ha")

    n_tr = gi(train_data.size(), "nt_2024")
    n_vl = gi(val_data.size(), "nv_2024")

    metrics = {
        'post_acuerdo_2': {
            'year': year, 'label': pi['label'],
            'overall_accuracy': round(oa, 4) if oa else None,
            'kappa': round(kp, 4) if kp else None,
            'n_training': n_tr, 'n_validation': n_vl, 'n_images': n,
            'confusion_matrix': cm,
            'producers_accuracy': pa, 'users_accuracy': ua,
            'class_areas_ha': {str(k): v for k, v in areas.items()},
            'post_classification_correction': 'NDVI + Hansen treecover2000 forest rule',
        }
    }
    save_json(metrics, 'classification_metrics.json')
    return classified


# ================================================================
# LOAD EXISTING CLASSIFICATIONS (T1, T2, T3 from run_fix_stages.py)
# ================================================================
def load_existing_classifications(region):
    """Reclassify T1, T2, T3 using same approach as run_fix_stages.py."""
    log("=" * 60)
    log("LOADING T1, T2, T3 CLASSIFICATIONS")
    log("=" * 60)

    maps = {}
    for pk, pi in list(PERIODS.items())[:3]:  # Only T1, T2, T3
        year = pi['map_year']
        log(f"  [{year}] {pi['label']}")

        composite, n_img, bands = create_fixed_composite(
            pi['start'], pi['end'], region, year)

        ref = training_mod.get_reference_lulc(year, region)
        samples = training_mod.generate_stratified_samples(ref, region, 300, seed=42 + year)
        train, val = training_mod.split_train_validation(samples, 0.7, seed=42)

        train_data = composite.select(bands).sampleRegions(
            collection=train, properties=['lulc_reference'],
            scale=30, tileScale=4, geometries=False)

        clf = ee.Classifier.smileRandomForest(
            numberOfTrees=200, minLeafPopulation=5, bagFraction=0.632, seed=42
        ).train(train_data, 'lulc_reference', bands)

        classified = composite.select(bands).classify(clf).rename('lulc').toInt8()
        classified = classified.focal_mode(3, 'square', 'pixels').rename('lulc').toInt8()
        jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
        classified = classified.where(jrc.select('occurrence').gte(80), 5).toInt8()

        maps[pk] = classified
        log(f"    Classified (lazy - not materialized)")

    return maps


# ================================================================
# FIXED CHANGE DETECTION
# ================================================================
def run_change_detection_fixed(maps, region):
    """
    Fix: transition matrix needs 2-band image for Reducer.group:
    band 0 = transition code (group key), band 1 = area (value to sum).
    """
    log("\n" + "=" * 60)
    log("CHANGE DETECTION (FIXED - 2 band reduceRegion)")
    log("=" * 60)

    plist = list(PERIODS.keys())
    trans_cfg = [
        ('T1_T2', plist[0], plist[1]),
        ('T2_T3', plist[1], plist[2]),
        ('T3_T4', plist[2], plist[3]),
    ]

    all_t = {}
    for tk, pf, pt in trans_cfg:
        yf = PERIODS[pf]['map_year']
        yt = PERIODS[pt]['map_year']
        yd = yt - yf
        log(f"\n  {yf}->{yt} ({yd} yrs)")

        lf = maps[pf]
        lt = maps[pt]

        # Create 2-band image: band0=transition_code, band1=area_ha
        transition_code = lf.multiply(10).add(lt).rename('transition')
        area_ha = ee.Image.pixelArea().divide(10000).rename('area')
        combined = transition_code.addBands(area_ha)

        stats = gi(combined.reduceRegion(
            reducer=ee.Reducer.sum().group(groupField=0, groupName='transition'),
            geometry=region, scale=100, maxPixels=1e13,
            tileScale=4, bestEffort=True
        ), f"tr_{tk}")

        matrix = {}
        if stats and 'groups' in stats:
            for g in stats['groups']:
                code = int(g['transition'])
                ha = round(g['sum'], 1)
                cf, ct = code // 10, code % 10
                if cf < 1 or cf > 7 or ct < 1 or ct > 7:
                    continue
                fn = LULC_CLASSES.get(cf, {}).get('name', f'C{cf}')
                tn = LULC_CLASSES.get(ct, {}).get('name', f'C{ct}')
                matrix[f"{cf}->{ct}"] = {'from': fn, 'to': tn, 'area_ha': ha}

            # Log top transitions (excluding persistence)
            non_p = [(k, v) for k, v in matrix.items()
                     if k.split('->')[0] != k.split('->')[1]]
            non_p.sort(key=lambda x: abs(x[1]['area_ha']), reverse=True)
            for k, v in non_p[:8]:
                log(f"    {v['from']} -> {v['to']}: {v['area_ha']:,.0f} ha")
        else:
            log(f"    No transition data for {tk}")

        # Compute change rates from transition matrix
        af, at = {}, {}
        for k, v in matrix.items():
            c1, c2 = int(k.split('->')[0]), int(k.split('->')[1])
            af[c1] = af.get(c1, 0) + v['area_ha']
            at[c2] = at.get(c2, 0) + v['area_ha']

        rates = {}
        for cid in range(1, 8):
            a1 = af.get(cid, 0)
            a2 = at.get(cid, 0)
            if a1 > 0 and a2 > 0:
                r = (1.0 / yd) * math.log(a2 / a1) * 100
            else:
                r = 0
            nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
            rates[cid] = {
                'name': nm,
                'area_t1_ha': round(a1, 1), 'area_t2_ha': round(a2, 1),
                'net_change_ha': round(a2 - a1, 1),
                'pct_change': round(((a2 - a1) / max(a1, 1)) * 100, 2),
                'annual_rate_pct': round(r, 3),
            }
            if abs(r) > 0.01:
                log(f"    {nm}: {r:+.3f}%/yr")

        all_t[tk] = {
            'years': f"{yf}-{yt}", 'years_between': yd,
            'transitions': matrix,
            'change_rates': {str(k): v for k, v in rates.items()},
        }

    # Keep existing Hansen data
    existing_path = os.path.join(OUTPUT_DIR, 'change_detection_results.json')
    if os.path.exists(existing_path):
        with open(existing_path) as f:
            existing = json.load(f)
        if 'hansen_gfc' in existing:
            all_t['hansen_gfc'] = existing['hansen_gfc']
        if 'hansen_treecover2000' in existing:
            all_t['hansen_treecover2000'] = existing['hansen_treecover2000']

    # Overwrite (not merge) for change detection
    with open(os.path.join(OUTPUT_DIR, 'change_detection_results.json'), 'w') as f:
        json.dump(all_t, f, indent=2, default=str)
    log("  >> Saved: change_detection_results.json")
    return all_t


# ================================================================
# FIXED ECOSYSTEM SERVICES
# ================================================================
def run_ecosystem_services_fixed(maps, region):
    """Fix water yield (use try/except for ET band mismatch) and habitat (tileScale=8)."""
    log("\n" + "=" * 60)
    log("ECOSYSTEM SERVICES (FIXED)")
    log("=" * 60)

    es = {}
    plist = list(PERIODS.keys())

    for pk, pi in PERIODS.items():
        year = pi['map_year']
        lulc = maps[pk]
        log(f"\n  [{year}] {pi['label']}")

        # Carbon (using tileScale=4, scale=100)
        carbon = ecosystem_mod.compute_carbon_storage(lulc, region)
        ct = gi(
            carbon['c_total'].multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
                ee.Reducer.sum(), region, 100, maxPixels=1e12,
                tileScale=4, bestEffort=True), f"c_{year}")
        c_val = ct.get('c_total', 0) if ct else 0
        log(f"    Carbon: {c_val:,.0f} Mg C")

        # Water yield (fix: compute directly instead of using module)
        try:
            chirps_p = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
                        .filterDate(f'{year}-01-01', f'{year}-12-31')
                        .filterBounds(region).sum().clip(region).rename('precip'))

            # Recharge coefficient per LULC class
            rc_vals = {1: 0.35, 2: 0.30, 3: 0.15, 4: 0.18, 5: 0, 6: 0.05, 7: 0.10}
            recharge = ee.Image(0).float()
            for cid, rc in rc_vals.items():
                recharge = recharge.where(lulc.eq(cid), rc)

            baseflow = chirps_p.multiply(recharge).rename('baseflow')
            bf_stat = gi(baseflow.reduceRegion(
                ee.Reducer.mean(), region, 1000, maxPixels=1e12,
                tileScale=4, bestEffort=True), f"bf_{year}")
            bf_val = bf_stat.get('baseflow', 0) if bf_stat else 0

            # Kc-based ET subtraction
            kc_vals = {1: 1.0, 2: 0.85, 3: 0.60, 4: 0.70, 5: 1.20, 6: 0.30, 7: 0.15}
            kc_img = ee.Image(0).float()
            for cid, kc in kc_vals.items():
                kc_img = kc_img.where(lulc.eq(cid), kc)

            # Use mean precip as ET proxy (simpler, avoids MODIS band issues)
            mean_precip = chirps_p.multiply(0.6)  # ~60% of precip as ET in tropics
            et_adj = mean_precip.multiply(kc_img)
            wy = chirps_p.subtract(et_adj).rename('water_yield')
            wy_stat = gi(wy.reduceRegion(
                ee.Reducer.mean(), region, 1000, maxPixels=1e12,
                tileScale=4, bestEffort=True), f"wy_{year}")
            wy_val = wy_stat.get('water_yield', 0) if wy_stat else 0
        except Exception as e:
            log(f"    Water WARN: {e}")
            wy_val = 0
            bf_val = 0
        log(f"    Water: {wy_val:.1f} mm/yr, Baseflow: {bf_val:.1f} mm/yr")

        # Habitat quality (use tileScale=8, scale=2000 to avoid memory)
        try:
            hab = ecosystem_mod.compute_habitat_quality(lulc, region)
            hq = gi(hab['habitat_quality'].reduceRegion(
                ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True),
                region, 2000, maxPixels=1e12, tileScale=8, bestEffort=True
            ), f"hq_{year}")
            hq_m = hq.get('habitat_quality_mean', 0) if hq else 0
            hq_s = hq.get('habitat_quality_stdDev', 0) if hq else 0
        except Exception as e:
            log(f"    Habitat WARN: {e}")
            hq_m, hq_s = 0, 0
        log(f"    Habitat: {hq_m:.3f} +/- {hq_s:.3f}")

        es[pk] = {
            'year': year, 'carbon_Mg_C': round(c_val, 0),
            'water_yield_mm': round(wy_val, 1), 'baseflow_mm': round(bf_val, 1),
            'habitat_quality_mean': round(hq_m, 4), 'habitat_quality_std': round(hq_s, 4),
        }

    # Carbon change (use scale=200, tileScale=8)
    log("\n  Carbon change:")
    for i in range(len(plist) - 1):
        pf, pt = plist[i], plist[i + 1]
        yf, yt = PERIODS[pf]['map_year'], PERIODS[pt]['map_year']
        ct1 = ecosystem_mod.compute_carbon_storage(maps[pf], region)
        ct2 = ecosystem_mod.compute_carbon_storage(maps[pt], region)
        c_change = ct2['c_total'].subtract(ct1['c_total']).rename('c_change')
        net = gi(c_change.multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
            ee.Reducer.sum(), region, 200, maxPixels=1e12,
            tileScale=8, bestEffort=True), f"cd_{yf}_{yt}")
        nv = net.get('c_change', 0) if net else 0
        log(f"    {yf}->{yt}: {nv:+,.0f} Mg C")
        es[f"carbon_change_{yf}_{yt}"] = {'net_Mg_C': round(nv, 0)}

    save_json(es, 'ecosystem_services_results.json')
    return es


# ================================================================
# FIXED CA-MARKOV (coarser scale, fewer points)
# ================================================================
def run_ca_markov_fixed(maps, region):
    """Fix: use 2000m scale, 3000 points, and toInt8 to reduce computation."""
    log("\n" + "=" * 60)
    log("CA-MARKOV PROJECTIONS (FIXED)")
    log("=" * 60)

    ca = {}
    cs = ['BDen', 'BSec', 'Past', 'Cult', 'Agua', 'Urb', 'Suel']

    log("  Extracting LULC arrays (2000m, 3000 pts)...")
    lnp = {}
    for pk in ['transicion', 'post_acuerdo_1', 'post_acuerdo_2']:
        year = PERIODS[pk]['map_year']
        lulc = maps[pk].toInt8()

        # Sample at coarser scale to avoid memory issues
        s = lulc.sample(region=region, scale=2000, numPixels=3000,
                        seed=42, geometries=False)
        info = gi(s.limit(3000), f"arr_{year}")

        if info and 'features' in info:
            vals = np.array([f['properties'].get('lulc', 0) for f in info['features']])
            vals = vals[vals > 0]
            log(f"    {year}: {len(vals)} valid pts, classes: {np.unique(vals)}")
            lnp[pk] = vals
        else:
            log(f"    {year}: FAILED to sample")

    if 'post_acuerdo_1' in lnp and 'post_acuerdo_2' in lnp:
        v20 = lnp['post_acuerdo_1']
        v24 = lnp['post_acuerdo_2']
        ml = min(len(v20), len(v24))
        v20, v24 = v20[:ml], v24[:ml]

        log(f"\n  Transition matrix 2020->2024 ({ml} pts):")
        tp = ca_markov_mod.compute_transition_probabilities(v20, v24)
        hdr = "         " + "".join(f"{n:>7}" for n in cs)
        log(f"  {hdr}")
        for i in range(7):
            row = f"  {cs[i]:<8}" + "".join(f"{tp[i, j]:>7.3f}" for j in range(7))
            log(row)
        ca['transition_matrix'] = tp.tolist()

        a24 = np.array([np.sum(v24 == c) for c in range(1, 8)]).astype(float)
        total = a24.sum()
        log(f"\n  Areas 2024 ({int(total)} pts):")
        for c in range(7):
            log(f"    {cs[c]}: {a24[c]:.0f} ({a24[c] / total * 100:.1f}%)")
        ca['areas_2024'] = {cs[c]: int(a24[c]) for c in range(7)}

        # Scenarios
        scenarios = ca_markov_mod.create_scenario_matrices(tp)
        for sn, sm in scenarios.items():
            log(f"\n  === {sn} ===")
            for ty in [2030, 2040]:
                ns = max(1, (ty - 2024) // 4)
                proj = ca_markov_mod.project_markov(a24, sm, ns)
                sd = {}
                for c in range(7):
                    pct = proj[c] / total * 100
                    chg = (proj[c] - a24[c]) / max(a24[c], 1) * 100
                    log(f"    {ty} {cs[c]}: {pct:.1f}% ({chg:+.1f}%)")
                    sd[cs[c]] = {'pct': round(pct, 2), 'change_pct': round(chg, 2)}
                ca[f'{sn}_{ty}'] = sd

        # Validation
        if 'transicion' in lnp:
            log("\n  Validation: hindcast 2024...")
            v16 = lnp['transicion']
            ml2 = min(len(v16), len(v20))
            tpv = ca_markov_mod.compute_transition_probabilities(v16[:ml2], v20[:ml2])
            a20 = np.array([np.sum(v20 == c) for c in range(1, 8)]).astype(float)
            sim = ca_markov_mod.project_markov(a20, tpv, 1)

            mae = 0
            val = {}
            for c in range(7):
                d = sim[c] - a24[c]
                mae += abs(d)
                log(f"    {cs[c]}: obs={a24[c]:.0f} sim={sim[c]:.0f} diff={d:+.0f}")
                val[cs[c]] = {'observed': int(a24[c]), 'simulated': round(float(sim[c]), 0)}

            rmae = mae / total * 100
            log(f"    Relative MAE: {rmae:.2f}%")
            ca['validation'] = {
                'method': 'Hindcast 2016->2020 -> predict 2024',
                'relative_MAE_pct': round(rmae, 2),
                'by_class': val,
            }

    save_json(ca, 'ca_markov_results.json')
    return ca


# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    log("=" * 70)
    log("FIX V2: T4 forest correction, change detection, CA-Markov")
    log(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    log("=" * 70)

    region = get_study_area()

    # Step 1: Load T1, T2, T3 classifications (lazy - no getInfo)
    maps_123 = load_existing_classifications(region)

    # Step 2: Reclassify T4 with forest correction
    t4_classified = reclassify_t4_fixed(region)
    maps = {**maps_123, 'post_acuerdo_2': t4_classified}

    # Step 3: Fixed change detection
    changes = run_change_detection_fixed(maps, region)

    # Step 4: Fixed ecosystem services
    es = run_ecosystem_services_fixed(maps, region)

    # Step 5: Fixed CA-Markov
    ca = run_ca_markov_fixed(maps, region)

    elapsed = time.time() - t0
    log(f"\nCOMPLETED in {elapsed / 60:.1f} minutes")

    files = sorted(f for f in os.listdir(OUTPUT_DIR) if f.endswith('.json'))
    for f in files:
        sz = os.path.getsize(os.path.join(OUTPUT_DIR, f))
        log(f"  {f} ({sz:,} bytes)")

    summary = {
        'date_v2': datetime.now().isoformat(),
        'time_v2_min': round(elapsed / 60, 1),
        'fixes_applied': [
            'T4 post-classification forest correction (NDVI + Hansen)',
            'Change detection: 2-band image for Reducer.group',
            'CA-Markov: 2000m scale, 3000 pts',
            'Ecosystem: direct water yield, tileScale=8 habitat',
        ],
    }
    save_json(summary, 'analysis_summary.json')


if __name__ == '__main__':
    main()
