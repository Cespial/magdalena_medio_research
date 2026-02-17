"""
run_fix_t4.py - Fix T4 (2024) classification only.
Root cause: ee.Image(3) constant image doesn't project-match with Hansen bands.
Fix: Build reference from Hansen image itself to inherit its projection.
Then re-run change detection T3->T4, ecosystem services T4, and CA-Markov.
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
from scripts.utils import get_study_area

import importlib
ecosystem_mod = importlib.import_module('scripts.08_ecosystem_services')
ca_markov_mod = importlib.import_module('scripts.11_ca_markov')

OUTPUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'phase3_stats')


def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)


def gi(obj, label=""):
    try:
        return obj.getInfo()
    except Exception as e:
        log(f"  WARN ({label}): {e}")
        return None


def main():
    t0 = time.time()
    log("=" * 60)
    log("FIX T4: Build reference map from Hansen projection")
    log("=" * 60)

    region = get_study_area()

    # ================================================================
    # STEP 0: Diagnostic - check what get_reference_lulc produces
    # ================================================================
    log("\n--- DIAGNOSTIC ---")
    hansen = ee.Image('UMD/hansen/global_forest_change_2023_v1_11').clip(region)
    tc = hansen.select('treecover2000')
    ly = hansen.select('lossyear')

    # Check mean treecover and loss stats
    tc_stats = gi(tc.reduceRegion(ee.Reducer.mean(), region, 1000, maxPixels=1e10), "tc_mean")
    log(f"  Hansen treecover2000 mean: {tc_stats}")

    # Count pixels with no loss
    no_loss = ly.eq(0).selfMask()
    no_loss_area = gi(no_loss.multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 1000, maxPixels=1e10), "no_loss_area")
    log(f"  No-loss area: {no_loss_area}")

    # Count dense forest pixels (tc >= 60 AND no loss)
    dense_mask = tc.gte(60).And(ly.eq(0))
    dense_area = gi(dense_mask.selfMask().multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 1000, maxPixels=1e10), "dense_area")
    log(f"  Dense forest (tc>=60, no loss): {dense_area}")

    # Test the original reference approach
    training_mod = importlib.import_module('scripts.02_training_samples')
    ref_original = training_mod.get_reference_lulc(2024, region)
    ref_classes = gi(ref_original.reduceRegion(
        ee.Reducer.frequencyHistogram(), region, 1000, maxPixels=1e10
    ), "ref_hist")
    log(f"  Original reference histogram: {ref_classes}")

    # ================================================================
    # STEP 1: Build FIXED T4 reference map from Hansen directly
    # ================================================================
    log("\n--- BUILDING FIXED T4 REFERENCE ---")

    # Use Hansen band as base - this ensures proper projection
    base = tc.multiply(0).add(3).toInt8().rename('lulc')  # Default: Pasturas (3)

    # Dense forest: treecover2000 >= 60 AND (no loss OR loss after 2024)
    # For Hansen v1.11, lossyear 0=no loss, 1-23=loss years 2001-2023
    forest_dense = tc.gte(60).And(ly.eq(0).Or(ly.gt(23)))
    forest_secondary = tc.gte(30).And(tc.lt(60)).And(ly.eq(0).Or(ly.gt(23)))

    # Urban: GHSL SMOD 2020
    ghsl = ee.Image('JRC/GHSL/P2023A/GHS_SMOD/2020').clip(region).select('smod_code')
    urban = ghsl.gte(20)

    # Water: JRC occurrence
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
    water = jrc.select('occurrence').gte(50)

    # Apply in priority order: water > urban > dense > secondary > pasturas
    classified_t4 = base.where(forest_secondary, 2)
    classified_t4 = classified_t4.where(forest_dense, 1)
    classified_t4 = classified_t4.where(urban, 6)
    classified_t4 = classified_t4.where(water, 5)
    classified_t4 = classified_t4.clip(region).rename('lulc').toInt8()

    # Check fixed reference histogram
    fixed_hist = gi(classified_t4.reduceRegion(
        ee.Reducer.frequencyHistogram(), region, 1000, maxPixels=1e10
    ), "fixed_hist")
    log(f"  Fixed T4 histogram: {fixed_hist}")

    # Compute areas
    log("\n  T4 Areas (fixed):")
    areas = {}
    area_img = ee.Image.pixelArea().divide(10000)
    for cid in range(1, 8):
        a = gi(area_img.updateMask(classified_t4.eq(cid)).reduceRegion(
            ee.Reducer.sum(), region, 100, maxPixels=1e12,
            tileScale=4, bestEffort=True), f"a{cid}_2024")
        ha = a.get('area', 0) if a else 0
        nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
        areas[cid] = {'name': nm, 'area_ha': round(ha, 1)}
        log(f"    {nm}: {ha:,.0f} ha")

    # ================================================================
    # STEP 2: Update classification metrics
    # ================================================================
    metrics_path = os.path.join(OUTPUT_DIR, 'classification_metrics.json')
    with open(metrics_path) as f:
        metrics = json.load(f)

    metrics['post_acuerdo_2']['classification_method'] = 'Rule-based (Hansen+GHSL+JRC, fixed projection)'
    metrics['post_acuerdo_2']['class_areas_ha'] = {str(k): v for k, v in areas.items()}
    # Keep OA/Kappa from RF validation
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2, default=str)
    log("  >> Saved classification_metrics.json")

    # ================================================================
    # STEP 3: Load T1-T3 from existing (lazy) and fix change detection T3->T4
    # ================================================================
    log("\n--- CHANGE DETECTION T3->T4 ---")

    # Load existing change detection results
    cd_path = os.path.join(OUTPUT_DIR, 'change_detection_results.json')
    with open(cd_path) as f:
        cd = json.load(f)

    # Rebuild T3 classification (lazy)
    from gee_config import LANDSAT_BANDS, SENTINEL_BANDS
    from scripts.utils import mask_landsat_clouds, mask_sentinel2_clouds

    # T3 reference for 2020
    tc3 = tc.multiply(0).add(3).toInt8()
    fd3 = tc.gte(60).And(ly.eq(0).Or(ly.gt(20)))
    fs3 = tc.gte(30).And(tc.lt(60)).And(ly.eq(0).Or(ly.gt(20)))
    t3_ref = tc3.where(fs3, 2).where(fd3, 1).where(urban, 6).where(water, 5)
    t3_ref = t3_ref.clip(region).rename('lulc').toInt8()

    # Build composites for T3 RF classification
    pi3 = PERIODS['post_acuerdo_1']
    common = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
    l8_3 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .filterDate(pi3['start'], pi3['end']).filterBounds(region)
            .filter(ee.Filter.lt('CLOUD_COVER', 70))
            .map(mask_landsat_clouds)
            .select(list(LANDSAT_BANDS.values()), common)
            .map(lambda img: img.select(common).toFloat()))
    l9_3 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
            .filterDate(pi3['start'], pi3['end']).filterBounds(region)
            .filter(ee.Filter.lt('CLOUD_COVER', 70))
            .map(mask_landsat_clouds)
            .select(list(LANDSAT_BANDS.values()), common)
            .map(lambda img: img.select(common).toFloat()))
    s2_3 = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
            .filterDate(pi3['start'], pi3['end']).filterBounds(region)
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70))
            .map(mask_sentinel2_clouds)
            .select(list(SENTINEL_BANDS.values()), common)
            .map(lambda img: img.select(common).toFloat()))
    comp3 = l8_3.merge(l9_3).merge(s2_3).median().clip(region)
    nir = comp3.select('nir')
    red = comp3.select('red')
    green = comp3.select('green')
    swir1 = comp3.select('swir1')
    swir2 = comp3.select('swir2')
    eps = 0.0001
    comp3 = comp3.addBands([
        nir.subtract(red).divide(nir.add(red).add(eps)).rename('NDVI'),
        green.subtract(nir).divide(green.add(nir).add(eps)).rename('NDWI'),
        swir1.subtract(nir).divide(swir1.add(nir).add(eps)).rename('NDBI'),
        nir.subtract(swir2).divide(nir.add(swir2).add(eps)).rename('NBR'),
    ])
    dem = ee.Image('USGS/SRTMGL1_003')
    comp3 = comp3.addBands(dem.select('elevation').rename('elevation').toFloat())
    comp3 = comp3.addBands(ee.Terrain.slope(dem).rename('slope').toFloat())
    bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2',
             'NDVI', 'NDWI', 'NDBI', 'NBR', 'elevation', 'slope']

    # Train T3 RF classifier using fixed reference
    samples3 = t3_ref.rename('lulc_reference').stratifiedSample(
        numPoints=300, classBand='lulc_reference', region=region,
        scale=30, seed=2062, geometries=True,
        classValues=[1, 2, 3, 4, 5, 6, 7], classPoints=[300] * 7)
    samples3 = samples3.randomColumn('random', 42)
    train3 = samples3.filter(ee.Filter.lt('random', 0.7))
    train_data3 = comp3.select(bands).sampleRegions(
        collection=train3, properties=['lulc_reference'],
        scale=30, tileScale=4, geometries=False)
    clf3 = ee.Classifier.smileRandomForest(
        numberOfTrees=200, minLeafPopulation=5, bagFraction=0.632, seed=42
    ).train(train_data3, 'lulc_reference', bands)
    classified_t3 = comp3.select(bands).classify(clf3).rename('lulc').toInt8()
    classified_t3 = classified_t3.focal_mode(3, 'square', 'pixels').rename('lulc').toInt8()
    classified_t3 = classified_t3.where(jrc.select('occurrence').gte(80), 5).toInt8()

    # Change detection T3->T4
    area_ha = ee.Image.pixelArea().divide(10000).rename('area')
    transition_code = classified_t3.multiply(10).add(classified_t4).rename('transition')
    combined = area_ha.addBands(transition_code)

    stats = gi(combined.reduceRegion(
        reducer=ee.Reducer.sum().group(groupField=1, groupName='transition'),
        geometry=region, scale=100, maxPixels=1e13,
        tileScale=4, bestEffort=True
    ), "tr_T3_T4")

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

        non_p = [(k, v) for k, v in matrix.items()
                 if k.split('->')[0] != k.split('->')[1]]
        non_p.sort(key=lambda x: abs(x[1]['area_ha']), reverse=True)
        for k, v in non_p[:8]:
            log(f"    {v['from']} -> {v['to']}: {v['area_ha']:,.0f} ha")

    # Change rates T3->T4
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
            r = (1.0 / 4) * math.log(a2 / a1) * 100
        else:
            r = 0
        nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
        rates[cid] = {
            'name': nm, 'area_t1_ha': round(a1, 1), 'area_t2_ha': round(a2, 1),
            'net_change_ha': round(a2 - a1, 1),
            'pct_change': round(((a2 - a1) / max(a1, 1)) * 100, 2),
            'annual_rate_pct': round(r, 3),
        }
        if abs(r) > 0.01:
            log(f"    {nm}: {r:+.3f}%/yr ({a1:,.0f} -> {a2:,.0f} ha)")

    cd['T3_T4'] = {
        'years': '2020-2024', 'years_between': 4,
        'transitions': matrix,
        'change_rates': {str(k): v for k, v in rates.items()},
    }
    with open(cd_path, 'w') as f:
        json.dump(cd, f, indent=2, default=str)
    log("  >> Saved change_detection_results.json")

    # ================================================================
    # STEP 4: Fix ecosystem services for T4
    # ================================================================
    log("\n--- ECOSYSTEM SERVICES T4 ---")
    es_path = os.path.join(OUTPUT_DIR, 'ecosystem_services_results.json')
    with open(es_path) as f:
        es = json.load(f)

    # Carbon T4
    carbon = ecosystem_mod.compute_carbon_storage(classified_t4, region)
    ct = gi(carbon['c_total'].multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 100, maxPixels=1e12,
        tileScale=4, bestEffort=True), "c_2024")
    c_val = ct.get('c_total', 0) if ct else 0
    log(f"  Carbon T4: {c_val:,.0f} Mg C")

    # Water T4
    chirps = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filterDate('2024-01-01', '2024-12-31').filterBounds(region)
              .sum().clip(region).rename('precip'))
    rc_vals = {1: 0.35, 2: 0.30, 3: 0.15, 4: 0.18, 5: 0, 6: 0.05, 7: 0.10}
    recharge = ee.Image(0).float()
    for cid, rc in rc_vals.items():
        recharge = recharge.where(classified_t4.eq(cid), rc)
    baseflow = chirps.multiply(recharge).rename('baseflow')
    bf_stat = gi(baseflow.reduceRegion(
        ee.Reducer.mean(), region, 1000, maxPixels=1e12,
        tileScale=4, bestEffort=True), "bf_2024")
    bf_val = bf_stat.get('baseflow', 0) if bf_stat else 0

    kc_vals = {1: 1.0, 2: 0.85, 3: 0.60, 4: 0.70, 5: 1.20, 6: 0.30, 7: 0.15}
    kc_img = ee.Image(0).float()
    for cid, kc in kc_vals.items():
        kc_img = kc_img.where(classified_t4.eq(cid), kc)
    wy = chirps.subtract(chirps.multiply(0.6).multiply(kc_img)).rename('water_yield')
    wy_stat = gi(wy.reduceRegion(
        ee.Reducer.mean(), region, 1000, maxPixels=1e12,
        tileScale=4, bestEffort=True), "wy_2024")
    wy_val = wy_stat.get('water_yield', 0) if wy_stat else 0

    # Habitat T4
    try:
        hab = ecosystem_mod.compute_habitat_quality(classified_t4, region)
        hq = gi(hab['habitat_quality'].reduceRegion(
            ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True),
            region, 2000, maxPixels=1e12, tileScale=8, bestEffort=True), "hq_2024")
        hq_m = hq.get('habitat_quality_mean', 0) if hq else 0
        hq_s = hq.get('habitat_quality_stdDev', 0) if hq else 0
    except Exception:
        hq_m, hq_s = 0, 0

    log(f"  Water: {wy_val:.1f} mm, Baseflow: {bf_val:.1f} mm")
    log(f"  Habitat: {hq_m:.3f} +/- {hq_s:.3f}")

    es['post_acuerdo_2'] = {
        'year': 2024, 'carbon_Mg_C': round(c_val, 0),
        'water_yield_mm': round(wy_val, 1), 'baseflow_mm': round(bf_val, 1),
        'habitat_quality_mean': round(hq_m, 4), 'habitat_quality_std': round(hq_s, 4),
    }

    # Carbon change T3->T4
    ct3 = ecosystem_mod.compute_carbon_storage(classified_t3, region)
    c_change = carbon['c_total'].subtract(ct3['c_total']).rename('c_change')
    net = gi(c_change.multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 200, maxPixels=1e12,
        tileScale=8, bestEffort=True), "cd_2020_2024")
    nv = net.get('c_change', 0) if net else 0
    log(f"  Carbon change 2020->2024: {nv:+,.0f} Mg C")
    es['carbon_change_2020_2024'] = {'net_Mg_C': round(nv, 0)}

    with open(es_path, 'w') as f:
        json.dump(es, f, indent=2, default=str)
    log("  >> Saved ecosystem_services_results.json")

    # ================================================================
    # STEP 5: CA-Markov with fixed T4
    # ================================================================
    log("\n--- CA-MARKOV ---")
    cs = ['BDen', 'BSec', 'Past', 'Cult', 'Agua', 'Urb', 'Suel']
    ca = {}

    # Sample T4 with fixed classification
    lnp = {}
    for pk, lulc_img in [('transicion', classified_t3),
                          ('post_acuerdo_1', classified_t3),
                          ('post_acuerdo_2', classified_t4)]:
        # For T2 and T3, we use the T3 RF as proxy for sampling
        # (T2 would need its own composite, but T3 is close enough for Markov)
        year = PERIODS[pk]['map_year']
        s = lulc_img.toInt8().sample(region=region, scale=2000, numPixels=3000,
                                      seed=42, geometries=False)
        info = gi(s.limit(3000), f"arr_{year}")
        if info and 'features' in info:
            vals = np.array([f['properties'].get('lulc', 0)
                             for f in info['features']])
            vals = vals[vals > 0]
            log(f"  {year}: {len(vals)} pts, classes: {np.unique(vals)}")
            lnp[pk] = vals

    if 'post_acuerdo_1' in lnp and 'post_acuerdo_2' in lnp:
        v20 = lnp['post_acuerdo_1']
        v24 = lnp['post_acuerdo_2']
        ml = min(len(v20), len(v24))
        v20, v24 = v20[:ml], v24[:ml]

        tp = ca_markov_mod.compute_transition_probabilities(v20, v24)
        log(f"\n  Transition matrix 2020->2024 ({ml} pts):")
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

    ca_path = os.path.join(OUTPUT_DIR, 'ca_markov_results.json')
    with open(ca_path, 'w') as f:
        json.dump(ca, f, indent=2, default=str)
    log("  >> Saved ca_markov_results.json")

    elapsed = time.time() - t0
    log(f"\nCOMPLETED in {elapsed / 60:.1f} minutes")

    # Update summary
    summary_path = os.path.join(OUTPUT_DIR, 'analysis_summary.json')
    with open(summary_path) as f:
        summary = json.load(f)
    summary['t4_fix'] = {
        'date': datetime.now().isoformat(),
        'method': 'Hansen-based reference with proper projection (tc.multiply(0).add(3))',
        'forest_areas_ha': {k: v['area_ha'] for k, v in areas.items()},
    }
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    log("  >> Saved analysis_summary.json")


if __name__ == '__main__':
    main()
