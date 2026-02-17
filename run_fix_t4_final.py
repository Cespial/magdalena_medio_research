"""
run_fix_t4_final.py - Fix T4 classification using Hansen GFC v1.12.

Root cause found: In Magdalena Medio, lossyear.eq(0) returns 0 for ALL pixels -
every pixel has experienced some tree cover loss between 2001-2023/2024.
For T1-T3, forest is captured by loss.gt(year_offset) (pixels lost AFTER that year).
For T4 (2024), loss.gt(24) returns 0 because v1.11 goes through 2023 only.

Fix: Use Hansen v1.12 (through 2024) with lossyear.gte(year_offset) instead of gt.
Pixels that lost forest IN 2024 still had forest entering 2024.
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
ecosystem_mod = importlib.import_module('scripts.08_ecosystem_services')
ca_markov_mod = importlib.import_module('scripts.11_ca_markov')
training_mod = importlib.import_module('scripts.02_training_samples')

OUTPUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'phase3_stats')

# Use v1.12 instead of v1.11
HANSEN_ASSET = 'UMD/hansen/global_forest_change_2024_v1_12'


def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)


def gi(obj, label=""):
    try:
        return obj.getInfo()
    except Exception as e:
        log(f"  WARN ({label}): {e}")
        return None


def get_reference_lulc_fixed(year, region):
    """
    Fixed reference map using Hansen v1.12 with gte instead of gt.
    Key change: loss.gte(year_offset) captures pixels losing forest IN the target year,
    meaning they still had forest entering that year.
    """
    hansen = ee.Image(HANSEN_ASSET).clip(region)
    treecover2000 = hansen.select('treecover2000')
    loss = hansen.select('lossyear')

    year_offset = year - 2000

    # gte instead of gt: include pixels that lost forest in the target year
    # (they still had forest entering that year)
    forest_dense = treecover2000.gte(60).And(
        loss.eq(0).Or(loss.gte(year_offset))
    )
    forest_secondary = treecover2000.gte(30).And(treecover2000.lt(60)).And(
        loss.eq(0).Or(loss.gte(year_offset))
    )

    # Use Hansen as base (for proper projection)
    reference = treecover2000.multiply(0).add(3).toInt8()  # Default: Pasturas
    reference = reference.where(forest_secondary, 2)
    reference = reference.where(forest_dense, 1)

    # Urban
    ghsl = ee.Image('JRC/GHSL/P2023A/GHS_SMOD/2020').clip(region).select('smod_code')
    reference = reference.where(ghsl.gte(20), 6)

    # Water
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
    reference = reference.where(jrc.select('occurrence').gte(50), 5)

    return reference.rename('lulc_reference')


def create_composite(start, end, region, year):
    from gee_config import LANDSAT_BANDS, SENTINEL_BANDS
    common = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']

    l8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
          .filterDate(start, end).filterBounds(region)
          .filter(ee.Filter.lt('CLOUD_COVER', 70))
          .map(mask_landsat_clouds)
          .select(list(LANDSAT_BANDS.values()), common)
          .map(lambda img: img.select(common).toFloat()))

    landsat = l8
    if year >= 2021:
        l9 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
              .filterDate(start, end).filterBounds(region)
              .filter(ee.Filter.lt('CLOUD_COVER', 70))
              .map(mask_landsat_clouds)
              .select(list(LANDSAT_BANDS.values()), common)
              .map(lambda img: img.select(common).toFloat()))
        landsat = l8.merge(l9)

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
    composite = composite.addBands([
        nir.subtract(red).divide(nir.add(red).add(eps)).rename('NDVI'),
        green.subtract(nir).divide(green.add(nir).add(eps)).rename('NDWI'),
        swir1.subtract(nir).divide(swir1.add(nir).add(eps)).rename('NDBI'),
        nir.subtract(swir2).divide(nir.add(swir2).add(eps)).rename('NBR'),
    ])
    dem = ee.Image('USGS/SRTMGL1_003')
    composite = composite.addBands(dem.select('elevation').rename('elevation').toFloat())
    composite = composite.addBands(ee.Terrain.slope(dem).rename('slope').toFloat())

    bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2',
             'NDVI', 'NDWI', 'NDBI', 'NBR', 'elevation', 'slope']
    return composite, n_images, bands


def classify_period(pk, region):
    """Classify a single period using fixed reference."""
    pi = PERIODS[pk]
    year = pi['map_year']
    log(f"\n  [{year}] {pi['label']}")

    composite, n_img, bands = create_composite(pi['start'], pi['end'], region, year)
    n = gi(n_img, f"n_{year}")
    log(f"    {n} images")

    ref = get_reference_lulc_fixed(year, region)

    # Check reference histogram
    hist = gi(ref.reduceRegion(
        ee.Reducer.frequencyHistogram(), region, 2000, maxPixels=1e10
    ), f"hist_{year}")
    log(f"    Ref histogram: {hist}")

    samples = ref.stratifiedSample(
        numPoints=300, classBand='lulc_reference', region=region,
        scale=30, seed=42 + year, geometries=True,
        classValues=[1, 2, 3, 4, 5, 6, 7], classPoints=[300] * 7)
    samples = samples.randomColumn('random', 42)
    train = samples.filter(ee.Filter.lt('random', 0.7))
    val = samples.filter(ee.Filter.gte('random', 0.7))

    train_data = composite.select(bands).sampleRegions(
        collection=train, properties=['lulc_reference'],
        scale=30, tileScale=4, geometries=False)
    val_data = composite.select(bands).sampleRegions(
        collection=val, properties=['lulc_reference'],
        scale=30, tileScale=4, geometries=False)

    clf = ee.Classifier.smileRandomForest(
        numberOfTrees=200, minLeafPopulation=5, bagFraction=0.632, seed=42
    ).train(train_data, 'lulc_reference', bands)

    imp = gi(ee.Dictionary(clf.explain().get('importance')), f"imp_{year}")
    if imp:
        top = sorted(imp.items(), key=lambda x: x[1], reverse=True)[:5]
        log(f"    Top: {', '.join(f'{k}:{v:.0f}' for k, v in top)}")

    classified = composite.select(bands).classify(clf).rename('lulc').toInt8()
    classified = classified.focal_mode(3, 'square', 'pixels').rename('lulc').toInt8()
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
    classified = classified.where(jrc.select('occurrence').gte(80), 5).toInt8()

    # Accuracy
    validated = val_data.classify(clf)
    em = validated.errorMatrix('lulc_reference', 'classification')
    oa = gi(em.accuracy(), f"oa_{year}")
    kp = gi(em.kappa(), f"kp_{year}")
    cm = gi(em.array(), f"cm_{year}")
    pa = gi(em.producersAccuracy(), f"pa_{year}")
    ua = gi(em.consumersAccuracy(), f"ua_{year}")
    n_tr = gi(train_data.size(), f"nt_{year}")
    n_vl = gi(val_data.size(), f"nv_{year}")

    if oa:
        log(f"    OA={oa:.4f} ({oa * 100:.1f}%), Kappa={kp:.4f}")

    # Areas
    log(f"    Areas:")
    areas = {}
    area_img = ee.Image.pixelArea().divide(10000)
    for cid in range(1, 8):
        a = gi(area_img.updateMask(classified.eq(cid)).reduceRegion(
            ee.Reducer.sum(), region, 100, maxPixels=1e12,
            tileScale=4, bestEffort=True), f"a{cid}_{year}")
        ha = a.get('area', 0) if a else 0
        nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
        areas[cid] = {'name': nm, 'area_ha': round(ha, 1)}
        log(f"      {nm}: {ha:,.0f} ha")

    metrics = {
        'year': year, 'label': pi['label'],
        'classification_method': 'Random Forest (Hansen v1.12 reference, gte fix)',
        'overall_accuracy': round(oa, 4) if oa else None,
        'kappa': round(kp, 4) if kp else None,
        'n_training': n_tr, 'n_validation': n_vl, 'n_images': n,
        'confusion_matrix': cm,
        'producers_accuracy': pa, 'users_accuracy': ua,
        'class_areas_ha': {str(k): v for k, v in areas.items()},
    }

    return classified, metrics, imp


def main():
    t0 = time.time()
    log("=" * 60)
    log("FIX T4 FINAL: Hansen v1.12 + gte year_offset")
    log("=" * 60)

    region = get_study_area()

    # ================================================================
    # STEP 0: Diagnostic with v1.12
    # ================================================================
    log("\n--- DIAGNOSTIC v1.12 ---")
    hansen = ee.Image(HANSEN_ASSET).clip(region)
    tc = hansen.select('treecover2000')
    ly = hansen.select('lossyear')

    tc_mean = gi(tc.reduceRegion(ee.Reducer.mean(), region, 1000, maxPixels=1e10), "tc")
    log(f"  treecover2000 mean: {tc_mean}")

    # Check lossyear histogram
    ly_hist = gi(ly.reduceRegion(
        ee.Reducer.frequencyHistogram(), region, 2000, maxPixels=1e10
    ), "ly_hist")
    log(f"  lossyear histogram: {ly_hist}")

    # Check no-loss area
    no_loss = ly.eq(0).selfMask().multiply(ee.Image.pixelArea().divide(10000))
    no_loss_ha = gi(no_loss.reduceRegion(ee.Reducer.sum(), region, 1000, maxPixels=1e10), "noloss")
    log(f"  No-loss area: {no_loss_ha}")

    # Check gte(24) area
    gte24 = ly.gte(24).selfMask().multiply(ee.Image.pixelArea().divide(10000))
    gte24_ha = gi(gte24.reduceRegion(ee.Reducer.sum(), region, 1000, maxPixels=1e10), "gte24")
    log(f"  lossyear >= 24 area: {gte24_ha}")

    # Check forest with fixed reference for 2024
    ref_2024 = get_reference_lulc_fixed(2024, region)
    ref_hist = gi(ref_2024.reduceRegion(
        ee.Reducer.frequencyHistogram(), region, 2000, maxPixels=1e10
    ), "ref_2024")
    log(f"  T4 reference histogram: {ref_hist}")

    # ================================================================
    # STEP 1: Re-classify T4 with v1.12
    # ================================================================
    log("\n--- CLASSIFYING T4 (2024) ---")
    t4_classified, t4_metrics, t4_imp = classify_period('post_acuerdo_2', region)

    # Also reclassify T3 for change detection consistency
    log("\n--- CLASSIFYING T3 (2020) ---")
    t3_classified, t3_metrics, t3_imp = classify_period('post_acuerdo_1', region)

    # Save metrics
    metrics_path = os.path.join(OUTPUT_DIR, 'classification_metrics.json')
    with open(metrics_path) as f:
        all_metrics = json.load(f)
    all_metrics['post_acuerdo_2'] = t4_metrics
    all_metrics['post_acuerdo_1'] = t3_metrics
    with open(metrics_path, 'w') as f:
        json.dump(all_metrics, f, indent=2, default=str)
    log("  >> Saved classification_metrics.json")

    # Save importance
    imp_path = os.path.join(OUTPUT_DIR, 'feature_importance.json')
    with open(imp_path) as f:
        all_imp = json.load(f)
    if t4_imp:
        all_imp['post_acuerdo_2'] = t4_imp
    if t3_imp:
        all_imp['post_acuerdo_1'] = t3_imp
    with open(imp_path, 'w') as f:
        json.dump(all_imp, f, indent=2, default=str)
    log("  >> Saved feature_importance.json")

    # ================================================================
    # STEP 2: Change detection T3->T4
    # ================================================================
    log("\n--- CHANGE DETECTION T3->T4 ---")
    area_ha = ee.Image.pixelArea().divide(10000).rename('area')
    transition_code = t3_classified.multiply(10).add(t4_classified).rename('transition')
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

        total_area = sum(v['area_ha'] for v in matrix.values())
        persist = sum(v['area_ha'] for k, v in matrix.items()
                      if k.split('->')[0] == k.split('->')[1])
        log(f"    Total: {total_area:,.0f} ha, Persistence: {persist / total_area * 100:.1f}%")

    # Change rates
    af, at = {}, {}
    for k, v in matrix.items():
        c1, c2 = int(k.split('->')[0]), int(k.split('->')[1])
        af[c1] = af.get(c1, 0) + v['area_ha']
        at[c2] = at.get(c2, 0) + v['area_ha']
    rates = {}
    for cid in range(1, 8):
        a1 = af.get(cid, 0)
        a2 = at.get(cid, 0)
        r = (1.0 / 4) * math.log(a2 / a1) * 100 if a1 > 0 and a2 > 0 else 0
        nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
        rates[cid] = {
            'name': nm, 'area_t1_ha': round(a1, 1), 'area_t2_ha': round(a2, 1),
            'net_change_ha': round(a2 - a1, 1),
            'pct_change': round(((a2 - a1) / max(a1, 1)) * 100, 2),
            'annual_rate_pct': round(r, 3),
        }
        if abs(r) > 0.01:
            log(f"    {nm}: {r:+.3f}%/yr ({a1:,.0f} -> {a2:,.0f} ha)")

    cd_path = os.path.join(OUTPUT_DIR, 'change_detection_results.json')
    with open(cd_path) as f:
        cd = json.load(f)
    cd['T3_T4'] = {
        'years': '2020-2024', 'years_between': 4,
        'transitions': matrix,
        'change_rates': {str(k): v for k, v in rates.items()},
    }
    with open(cd_path, 'w') as f:
        json.dump(cd, f, indent=2, default=str)
    log("  >> Saved change_detection_results.json")

    # ================================================================
    # STEP 3: Ecosystem services T4
    # ================================================================
    log("\n--- ECOSYSTEM SERVICES T4 ---")
    carbon = ecosystem_mod.compute_carbon_storage(t4_classified, region)
    ct = gi(carbon['c_total'].multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 100, maxPixels=1e12,
        tileScale=4, bestEffort=True), "c_2024")
    c_val = ct.get('c_total', 0) if ct else 0
    log(f"  Carbon: {c_val:,.0f} Mg C")

    chirps = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filterDate('2024-01-01', '2024-12-31').filterBounds(region)
              .sum().clip(region).rename('precip'))
    rc = {1: 0.35, 2: 0.30, 3: 0.15, 4: 0.18, 5: 0, 6: 0.05, 7: 0.10}
    recharge = ee.Image(0).float()
    for cid, rv in rc.items():
        recharge = recharge.where(t4_classified.eq(cid), rv)
    bf = gi(chirps.multiply(recharge).rename('bf').reduceRegion(
        ee.Reducer.mean(), region, 1000, maxPixels=1e12,
        tileScale=4, bestEffort=True), "bf")
    bf_val = bf.get('bf', 0) if bf else 0

    kc = {1: 1.0, 2: 0.85, 3: 0.60, 4: 0.70, 5: 1.20, 6: 0.30, 7: 0.15}
    kc_img = ee.Image(0).float()
    for cid, kv in kc.items():
        kc_img = kc_img.where(t4_classified.eq(cid), kv)
    wy = chirps.subtract(chirps.multiply(0.6).multiply(kc_img)).rename('wy')
    wy_stat = gi(wy.reduceRegion(
        ee.Reducer.mean(), region, 1000, maxPixels=1e12,
        tileScale=4, bestEffort=True), "wy")
    wy_val = wy_stat.get('wy', 0) if wy_stat else 0

    try:
        hab = ecosystem_mod.compute_habitat_quality(t4_classified, region)
        hq = gi(hab['habitat_quality'].reduceRegion(
            ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True),
            region, 2000, maxPixels=1e12, tileScale=8, bestEffort=True), "hq")
        hq_m = hq.get('habitat_quality_mean', 0) if hq else 0
        hq_s = hq.get('habitat_quality_stdDev', 0) if hq else 0
    except Exception:
        hq_m, hq_s = 0, 0

    log(f"  Water: {wy_val:.1f} mm, Baseflow: {bf_val:.1f} mm, Habitat: {hq_m:.3f}")

    # Carbon change T3->T4
    ct3 = ecosystem_mod.compute_carbon_storage(t3_classified, region)
    c_change = carbon['c_total'].subtract(ct3['c_total']).rename('c_change')
    net = gi(c_change.multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
        ee.Reducer.sum(), region, 200, maxPixels=1e12,
        tileScale=8, bestEffort=True), "cc")
    nv = net.get('c_change', 0) if net else 0
    log(f"  Carbon change T3->T4: {nv:+,.0f} Mg C")

    es_path = os.path.join(OUTPUT_DIR, 'ecosystem_services_results.json')
    with open(es_path) as f:
        es = json.load(f)
    es['post_acuerdo_2'] = {
        'year': 2024, 'carbon_Mg_C': round(c_val, 0),
        'water_yield_mm': round(wy_val, 1), 'baseflow_mm': round(bf_val, 1),
        'habitat_quality_mean': round(hq_m, 4), 'habitat_quality_std': round(hq_s, 4),
    }
    es['carbon_change_2020_2024'] = {'net_Mg_C': round(nv, 0)}
    with open(es_path, 'w') as f:
        json.dump(es, f, indent=2, default=str)
    log("  >> Saved ecosystem_services_results.json")

    # ================================================================
    # STEP 4: CA-Markov
    # ================================================================
    log("\n--- CA-MARKOV ---")
    cs = ['BDen', 'BSec', 'Past', 'Cult', 'Agua', 'Urb', 'Suel']
    ca = {}

    lnp = {}
    for pk, lulc_img, yr in [
        ('transicion', t3_classified, 2016),
        ('post_acuerdo_1', t3_classified, 2020),
        ('post_acuerdo_2', t4_classified, 2024)
    ]:
        s = lulc_img.toInt8().sample(region=region, scale=2000, numPixels=3000,
                                      seed=42, geometries=False)
        info = gi(s.limit(3000), f"arr_{yr}")
        if info and 'features' in info:
            vals = np.array([f['properties'].get('lulc', 0) for f in info['features']])
            vals = vals[vals > 0]
            log(f"  {yr}: {len(vals)} pts, classes: {np.unique(vals)}")
            lnp[pk] = vals

    if 'post_acuerdo_1' in lnp and 'post_acuerdo_2' in lnp:
        v20 = lnp['post_acuerdo_1']
        v24 = lnp['post_acuerdo_2']
        ml = min(len(v20), len(v24))
        v20, v24 = v20[:ml], v24[:ml]

        tp = ca_markov_mod.compute_transition_probabilities(v20, v24)
        log(f"\n  Transition matrix ({ml} pts):")
        for i in range(7):
            row = f"  {cs[i]:<8}" + "".join(f"{tp[i, j]:>7.3f}" for j in range(7))
            log(row)
        ca['transition_matrix'] = tp.tolist()

        a24 = np.array([np.sum(v24 == c) for c in range(1, 8)]).astype(float)
        total = a24.sum()
        for c in range(7):
            log(f"  {cs[c]}: {a24[c]:.0f} ({a24[c] / total * 100:.1f}%)")
        ca['areas_2024'] = {cs[c]: int(a24[c]) for c in range(7)}

        scenarios = ca_markov_mod.create_scenario_matrices(tp)
        for sn, sm in scenarios.items():
            for ty in [2030, 2040]:
                ns = max(1, (ty - 2024) // 4)
                proj = ca_markov_mod.project_markov(a24, sm, ns)
                sd = {}
                for c in range(7):
                    pct = proj[c] / total * 100
                    chg = (proj[c] - a24[c]) / max(a24[c], 1) * 100
                    sd[cs[c]] = {'pct': round(pct, 2), 'change_pct': round(chg, 2)}
                ca[f'{sn}_{ty}'] = sd
                log(f"  {sn} {ty}: " + ", ".join(f"{cs[c]}:{sd[cs[c]]['pct']:.1f}%"
                                                   for c in range(7) if sd[cs[c]]['pct'] > 0.5))

    ca_path = os.path.join(OUTPUT_DIR, 'ca_markov_results.json')
    with open(ca_path, 'w') as f:
        json.dump(ca, f, indent=2, default=str)
    log("  >> Saved ca_markov_results.json")

    elapsed = time.time() - t0
    log(f"\nCOMPLETED in {elapsed / 60:.1f} minutes")

    # Summary
    summary_path = os.path.join(OUTPUT_DIR, 'analysis_summary.json')
    with open(summary_path) as f:
        summary = json.load(f)
    summary['t4_final_fix'] = {
        'date': datetime.now().isoformat(),
        'hansen_version': 'v1.12 (through 2024)',
        'fix': 'lossyear.gte(year_offset) instead of gt - captures pixels losing forest in target year',
        't4_areas': {str(k): v for k, v in t4_metrics['class_areas_ha'].items()},
    }
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)


if __name__ == '__main__':
    main()
