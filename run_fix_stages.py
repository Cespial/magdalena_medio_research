"""
run_fix_stages.py - Fix for Sentinel-2 band type mismatch and re-run failed stages.
Issue: S2 images with inconsistent band ranges when merged with Landsat.
Fix: Cast all bands to Float32 before median, and use .toInt8() for classified outputs.
"""

import ee
import os
import sys
import json
import time
import numpy as np
from datetime import datetime

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_DIR)
from gee_config import PERIODS, LULC_CLASSES, CARBON_POOLS
from scripts.utils import get_study_area, mask_landsat_clouds, mask_sentinel2_clouds

import importlib
training_mod = importlib.import_module('scripts.02_training_samples')
change_mod = importlib.import_module('scripts.05_change_detection')
hotspot_mod = importlib.import_module('scripts.07_hotspot_analysis')
ecosystem_mod = importlib.import_module('scripts.08_ecosystem_services')
gwr_mod = importlib.import_module('scripts.10_gwr_drivers')
ca_markov_mod = importlib.import_module('scripts.11_ca_markov')

OUTPUT_DIR = os.path.join(PROJECT_DIR, 'outputs', 'phase3_stats')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)


def save_json(data, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    # Merge with existing data if file exists
    if os.path.exists(path):
        with open(path) as f:
            existing = json.load(f)
        existing.update(data)
        data = existing
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
# FIXED COMPOSITE - explicit Float32 cast
# ================================================================
def create_fixed_composite(start, end, region, year):
    """Composite with explicit Float32 casting to prevent band type mismatch."""
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

    # Indices
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
# RE-CLASSIFY periods with band type fix
# ================================================================
def reclassify_all(region):
    log("=" * 60)
    log("RE-CLASSIFYING ALL PERIODS (fixed band types)")
    log("=" * 60)

    results = {}
    metrics = {}
    importance_all = {}

    for pk, pi in PERIODS.items():
        year = pi['map_year']
        t0 = time.time()
        log(f"\n  [{year}] {pi['label']}")

        composite, n_img, bands = create_fixed_composite(
            pi['start'], pi['end'], region, year
        )
        n = gi(n_img, f"n_{year}")
        log(f"    {n} images, {len(bands)} bands")

        # Training
        ref = training_mod.get_reference_lulc(year, region)
        samples = training_mod.generate_stratified_samples(ref, region, 300, seed=42+year)
        train, val = training_mod.split_train_validation(samples, 0.7, seed=42)

        train_data = composite.select(bands).sampleRegions(
            collection=train, properties=['lulc_reference'],
            scale=30, tileScale=4, geometries=False
        )
        val_data = composite.select(bands).sampleRegions(
            collection=val, properties=['lulc_reference'],
            scale=30, tileScale=4, geometries=False
        )

        # RF
        clf = ee.Classifier.smileRandomForest(
            numberOfTrees=200, minLeafPopulation=5, bagFraction=0.632, seed=42
        ).train(train_data, 'lulc_reference', bands)

        imp = gi(ee.Dictionary(clf.explain().get('importance')), f"imp_{year}")
        if imp:
            top = sorted(imp.items(), key=lambda x: x[1], reverse=True)[:5]
            log(f"    Top: {', '.join(f'{k}:{v:.0f}' for k,v in top)}")
            importance_all[pk] = imp

        # Classify and simplify
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

        if oa:
            log(f"    OA={oa:.4f} ({oa*100:.1f}%), Kappa={kp:.4f}")
        else:
            log(f"    Accuracy deferred")

        # Areas
        log(f"    Areas...")
        areas = {}
        area_img = ee.Image.pixelArea().divide(10000)
        for cid in range(1, 8):
            a = gi(
                area_img.updateMask(classified.eq(cid)).reduceRegion(
                    ee.Reducer.sum(), region, 100, maxPixels=1e12,
                    tileScale=4, bestEffort=True
                ), f"a{cid}_{year}"
            )
            ha = a.get('area', 0) if a else 0
            nm = LULC_CLASSES.get(cid, {}).get('name', f'C{cid}')
            areas[cid] = {'name': nm, 'area_ha': round(ha, 1)}
            log(f"      {nm}: {ha:,.0f} ha")

        n_tr = gi(train_data.size(), f"nt_{year}")
        n_vl = gi(val_data.size(), f"nv_{year}")

        elapsed = time.time() - t0
        log(f"    Time: {elapsed:.0f}s")

        results[pk] = {'classified': classified}
        metrics[pk] = {
            'year': year, 'label': pi['label'],
            'overall_accuracy': round(oa, 4) if oa else None,
            'kappa': round(kp, 4) if kp else None,
            'n_training': n_tr, 'n_validation': n_vl, 'n_images': n,
            'confusion_matrix': cm,
            'producers_accuracy': pa, 'users_accuracy': ua,
            'class_areas_ha': {str(k): v for k, v in areas.items()},
        }

    save_json(metrics, 'classification_metrics.json')
    save_json(importance_all, 'feature_importance.json')
    return results, metrics


# ================================================================
# RE-RUN CHANGE DETECTION
# ================================================================
def rerun_change_detection(maps, region):
    log("\n" + "=" * 60)
    log("CHANGE DETECTION (with simplified classified images)")
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

        lf = maps[pf]['classified']
        lt = maps[pt]['classified']

        _, stats = change_mod.compute_transition_matrix(lf, lt, region, scale=100)
        si = gi(stats, f"tr_{tk}")

        matrix = {}
        if si and 'groups' in si:
            for g in si['groups']:
                code = int(g['transition'])
                ha = round(g['sum'], 1)
                cf, ct = code // 10, code % 10
                fn = LULC_CLASSES.get(cf, {}).get('name', f'C{cf}')
                tn = LULC_CLASSES.get(ct, {}).get('name', f'C{ct}')
                matrix[f"{cf}->{ct}"] = {'from': fn, 'to': tn, 'area_ha': ha}

            non_p = [(k, v) for k, v in matrix.items() if k.split('->')[0] != k.split('->')[1]]
            non_p.sort(key=lambda x: abs(x[1]['area_ha']), reverse=True)
            for k, v in non_p[:6]:
                log(f"    {v['from']} -> {v['to']}: {v['area_ha']:,.0f} ha")

        af, at = {}, {}
        for k, v in matrix.items():
            c1, c2 = int(k.split('->')[0]), int(k.split('->')[1])
            af[c1] = af.get(c1, 0) + v['area_ha']
            at[c2] = at.get(c2, 0) + v['area_ha']

        rates = change_mod.compute_change_rates(af, at, yd)
        for cid, r in rates.items():
            if abs(r['annual_rate_pct']) > 0.01:
                log(f"    {r['name']}: {r['annual_rate_pct']:+.3f}%/yr")

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

    save_json(all_t, 'change_detection_results.json')
    return all_t


# ================================================================
# RE-RUN ECOSYSTEM SERVICES
# ================================================================
def rerun_ecosystem_services(maps, region):
    log("\n" + "=" * 60)
    log("ECOSYSTEM SERVICES")
    log("=" * 60)

    es = {}
    plist = list(PERIODS.keys())

    for pk, pi in PERIODS.items():
        year = pi['map_year']
        lulc = maps[pk]['classified']
        log(f"\n  [{year}] {pi['label']}")

        # Carbon
        carbon = ecosystem_mod.compute_carbon_storage(lulc, region)
        ct = gi(
            carbon['c_total'].multiply(ee.Image.pixelArea().divide(10000)).reduceRegion(
                ee.Reducer.sum(), region, 100, maxPixels=1e12,
                tileScale=4, bestEffort=True
            ), f"c_{year}"
        )
        c_val = ct.get('c_total', 0) if ct else 0
        log(f"    Carbon: {c_val:,.0f} Mg C")

        # Water
        water = ecosystem_mod.compute_water_yield_proxy(lulc, region, year)
        wy = gi(water['water_yield'].reduceRegion(
            ee.Reducer.mean(), region, 1000, maxPixels=1e12,
            tileScale=4, bestEffort=True
        ), f"w_{year}")
        wy_val = wy.get('water_yield', 0) if wy else 0

        bf = gi(water['baseflow'].reduceRegion(
            ee.Reducer.mean(), region, 1000, maxPixels=1e12,
            tileScale=4, bestEffort=True
        ), f"bf_{year}")
        bf_val = bf.get('baseflow', 0) if bf else 0
        log(f"    Water: {wy_val:.1f} mm/yr, Baseflow: {bf_val:.1f} mm/yr")

        # Habitat
        hab = ecosystem_mod.compute_habitat_quality(lulc, region)
        hq = gi(hab['habitat_quality'].reduceRegion(
            ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True),
            region, 1000, maxPixels=1e12, tileScale=4, bestEffort=True
        ), f"hq_{year}")
        hq_m = hq.get('habitat_quality_mean', 0) if hq else 0
        hq_s = hq.get('habitat_quality_stdDev', 0) if hq else 0
        log(f"    Habitat: {hq_m:.3f} +/- {hq_s:.3f}")

        es[pk] = {
            'year': year, 'carbon_Mg_C': round(c_val, 0),
            'water_yield_mm': round(wy_val, 1), 'baseflow_mm': round(bf_val, 1),
            'habitat_quality_mean': round(hq_m, 4), 'habitat_quality_std': round(hq_s, 4),
        }

    # Carbon change
    log("\n  Carbon change:")
    for i in range(len(plist) - 1):
        pf, pt = plist[i], plist[i+1]
        yf, yt = PERIODS[pf]['map_year'], PERIODS[pt]['map_year']
        cc = ecosystem_mod.compute_carbon_change(
            maps[pf]['classified'], maps[pt]['classified'], region
        )
        net = gi(cc['net_change_Mg_C'], f"cd_{yf}_{yt}")
        nv = net.get('c_change', 0) if net else 0
        log(f"    {yf}->{yt}: {nv:+,.0f} Mg C")
        es[f"carbon_change_{yf}_{yt}"] = {'net_Mg_C': round(nv, 0)}

    save_json(es, 'ecosystem_services_results.json')
    return es


# ================================================================
# RE-RUN HOTSPOT ANALYSIS
# ================================================================
def rerun_hotspot(maps, region):
    log("\n" + "=" * 60)
    log("HOTSPOT ANALYSIS")
    log("=" * 60)

    hs = {}
    lulc_t1 = maps['pre_acuerdo']['classified']
    lulc_t4 = maps['post_acuerdo_2']['classified']
    forest_t1 = lulc_t1.eq(1).Or(lulc_t1.eq(2))
    forest_t4 = lulc_t4.eq(1).Or(lulc_t4.eq(2))
    defor = forest_t1.And(forest_t4.Not()).toInt8().rename('deforestation')

    combined = defor.addBands(ee.Image.pixelLonLat())

    log("  Sampling (1500 pts, 1km)...")
    sample = combined.sample(region=region, scale=1000, numPixels=1500,
                             seed=42, geometries=True)
    info = gi(sample.limit(1500), "hs")

    if info and 'features' in info:
        feats = info['features']
        n = len(feats)
        log(f"  {n} points")

        if n > 30:
            vals = np.array([f['properties'].get('deforestation', 0) for f in feats])
            coords = np.array([
                [f['properties'].get('longitude', 0), f['properties'].get('latitude', 0)]
                for f in feats
            ])

            W = hotspot_mod.create_queen_weights(n, coords)

            morans = hotspot_mod.compute_morans_i(vals, W)
            log(f"  Moran I={morans['I']:.4f}, z={morans['z_score']:.2f}, p={morans['p_value']:.6f}")
            hs['morans_i'] = morans

            gi_star = hotspot_mod.compute_getis_ord_gi_star(vals, W)
            cats = hotspot_mod.classify_hotspots(gi_star)
            counts = {}
            for nm, v in [('hotspot_99',3),('hotspot_95',2),('hotspot_90',1),
                          ('not_significant',0),('coldspot_90',-1),
                          ('coldspot_95',-2),('coldspot_99',-3)]:
                counts[nm] = int(np.sum(cats == v))
            log(f"  HS99:{counts['hotspot_99']} HS95:{counts['hotspot_95']} "
                f"NS:{counts['not_significant']} CS99:{counts['coldspot_99']}")
            hs['gi_star'] = counts
            hs['defor_rate_mean'] = round(float(np.mean(vals)), 4)

            # KDE
            dp = coords[vals > 0]
            if len(dp) > 10:
                density, xg, yg = hotspot_mod.compute_kernel_density(dp, 80, 5000)
                hs['kde'] = {'n_points': len(dp), 'max_density': round(float(np.max(density)), 6)}

    save_json(hs, 'hotspot_analysis_results.json')
    return hs


# ================================================================
# RE-RUN GWR (fix SoilGrids path)
# ================================================================
def rerun_gwr(maps, region):
    log("\n" + "=" * 60)
    log("GWR DRIVERS")
    log("=" * 60)

    res = {}
    lulc_16 = maps['transicion']['classified']
    lulc_24 = maps['post_acuerdo_2']['classified']
    f16 = lulc_16.eq(1).Or(lulc_16.eq(2))
    f24 = lulc_24.eq(1).Or(lulc_24.eq(2))
    defor = f16.And(f24.Not()).toInt8().rename('defor_rate')

    # Build driver stack without SoilGrids clay (asset not available)
    srtm = ee.Image('USGS/SRTMGL1_003').clip(region)
    elevation = srtm.select('elevation').rename('elevation').toFloat()
    slope = ee.Terrain.slope(srtm).rename('slope').toFloat()

    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(region)
    rivers = jrc.select('occurrence').gte(50)
    dist_rivers = rivers.fastDistanceTransform().sqrt().multiply(30).rename('dist_rivers').toFloat()

    ghsl_b = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2020').clip(region)
    built = ghsl_b.select('built_surface').gt(0)
    dist_roads = built.fastDistanceTransform().sqrt().multiply(30).rename('dist_roads').toFloat()

    smod = ee.Image('JRC/GHSL/P2023A/GHS_SMOD/2020').clip(region)
    urban = smod.select('smod_code').gte(20)
    dist_urban = urban.fastDistanceTransform().sqrt().multiply(30).rename('dist_urban').toFloat()

    worldpop = (ee.ImageCollection('WorldPop/GP/100m/pop')
                .filter(ee.Filter.eq('country', 'COL'))
                .filter(ee.Filter.eq('year', 2020)).first()
                .clip(region).rename('pop_density').toFloat())

    chirps = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filterDate('2012-01-01', '2024-12-31').filterBounds(region)
              .mean().multiply(365).clip(region).rename('precip').toFloat())

    lst = (ee.ImageCollection('MODIS/061/MOD11A2')
           .filterDate('2012-01-01', '2024-12-31').filterBounds(region)
           .select('LST_Day_1km').mean().multiply(0.02).subtract(273.15)
           .clip(region).rename('lst').toFloat())

    drivers = (elevation.addBands(slope).addBands(dist_rivers)
               .addBands(dist_roads).addBands(dist_urban)
               .addBands(worldpop).addBands(chirps).addBands(lst))

    combined = drivers.addBands(defor)

    log("  Sampling (1500 pts, 1km)...")
    sample = combined.sample(region=region, scale=1000, numPixels=1500,
                             seed=42, geometries=True)
    info = gi(sample.limit(1500), "gwr")

    var_names = ['elevation', 'slope', 'dist_rivers', 'dist_roads',
                 'dist_urban', 'pop_density', 'precip', 'lst']

    if info and 'features' in info:
        feats = info['features']
        log(f"  {len(feats)} features")

        y_l, X_l, c_l = [], [], []
        for f in feats:
            p = f['properties']
            yv = p.get('defor_rate')
            if yv is None:
                continue
            row = []
            ok = True
            for v in var_names:
                val = p.get(v)
                if val is None:
                    ok = False
                    break
                row.append(val)
            if ok and f.get('geometry'):
                y_l.append(yv)
                X_l.append(row)
                c_l.append(f['geometry']['coordinates'][:2])

        nv = len(y_l)
        log(f"  Valid: {nv} points")

        if nv > 50:
            y = np.array(y_l)
            X = np.array(X_l)
            coords = np.array(c_l)

            Xm = np.mean(X, axis=0)
            Xs = np.std(X, axis=0)
            Xs[Xs == 0] = 1
            Xn = (X - Xm) / Xs

            # VIF
            vifs = gwr_mod.compute_vif(Xn)
            res['vif'] = {var_names[i]: vifs[i] for i in range(len(var_names))}
            log("  VIF: " + ", ".join(f"{v}:{vifs[i]:.1f}" for i, v in enumerate(var_names)))

            # OLS
            ols = gwr_mod.fit_ols(Xn, y)
            log(f"  OLS: R2={ols['r_squared']:.4f}, AIC={ols['aic']:.2f}")

            ols_names = ['intercept'] + var_names
            res['ols'] = {
                'r2': ols['r_squared'], 'adj_r2': ols['adj_r_squared'],
                'aic': ols['aic'], 'n': ols['n'],
                'coefficients': {ols_names[i]: round(ols['coefficients'][i], 6)
                                 for i in range(len(ols_names))},
                't_statistics': {ols_names[i]: round(ols['t_statistics'][i], 4)
                                 for i in range(len(ols_names))},
            }
            for i, nm in enumerate(ols_names):
                log(f"    {nm:<20} b={ols['coefficients'][i]:>8.4f} t={ols['t_statistics'][i]:>7.2f}")

            # GWR
            log("  GWR bandwidth optimization...")
            bw, baic = gwr_mod.optimize_bandwidth(Xn, y, coords, n_steps=10)
            log(f"  GWR: bw={bw}, AICc={baic:.2f}")

            gwr = gwr_mod.compute_gwr(Xn, y, coords, bandwidth=bw)
            log(f"  GWR: mean_R2={gwr['mean_r2']:.4f}, AIC={gwr['aic']:.2f}")

            summary = gwr_mod.summarize_gwr_results(gwr, var_names)
            for nm, st in summary.items():
                log(f"    {nm:<20} mean={st['mean']:>8.4f} %pos={st['pct_positive']:>5.1f}%")

            res['gwr'] = {
                'bandwidth': bw, 'mean_r2': gwr['mean_r2'],
                'median_r2': gwr['median_r2'], 'aic': gwr['aic'],
                'summary': summary,
            }

            comp = gwr_mod.compare_ols_gwr(ols, gwr)
            log(f"  OLS vs GWR: dR2={comp['r2_improvement']:.4f}, dAIC={comp['aic_improvement']:.2f}")
            res['comparison'] = comp

    save_json(res, 'gwr_drivers_results.json')
    return res


# ================================================================
# RE-RUN CA-MARKOV
# ================================================================
def rerun_ca_markov(maps, region):
    log("\n" + "=" * 60)
    log("CA-MARKOV PROJECTIONS")
    log("=" * 60)

    ca = {}
    cs = ['BDen', 'BSec', 'Past', 'Cult', 'Agua', 'Urb', 'Suel']

    log("  Extracting LULC arrays (500m, 8000 pts)...")
    lnp = {}
    for pk in ['transicion', 'post_acuerdo_1', 'post_acuerdo_2']:
        year = PERIODS[pk]['map_year']
        lulc = maps[pk]['classified']
        s = lulc.sample(region=region, scale=500, numPixels=8000,
                        seed=42, geometries=False)
        info = gi(s.limit(8000), f"arr_{year}")
        if info and 'features' in info:
            vals = np.array([f['properties'].get('lulc', 0) for f in info['features']])
            vals = vals[vals > 0]  # Remove nodata
            log(f"    {year}: {len(vals)} valid pts, classes: {np.unique(vals)}")
            lnp[pk] = vals

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
            row = f"  {cs[i]:<8}" + "".join(f"{tp[i,j]:>7.3f}" for j in range(7))
            log(row)
        ca['transition_matrix'] = tp.tolist()

        a24 = np.array([np.sum(v24 == c) for c in range(1, 8)]).astype(float)
        total = a24.sum()
        log(f"\n  Areas 2024 ({int(total)} pts):")
        for c in range(7):
            log(f"    {cs[c]}: {a24[c]:.0f} ({a24[c]/total*100:.1f}%)")
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
    log("FIX & RE-RUN: Stages with band type errors")
    log(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    log("=" * 70)

    region = get_study_area()

    # Stage 1: Re-classify all (with fixed band types)
    maps, metrics = reclassify_all(region)

    # Stage 2: Change detection
    changes = rerun_change_detection(maps, region)

    # Stage 3: Ecosystem services
    es = rerun_ecosystem_services(maps, region)

    # Stage 5: Hotspot
    hs = rerun_hotspot(maps, region)

    # Stage 6: GWR
    gwr = rerun_gwr(maps, region)

    # Stage 7: CA-Markov
    ca = rerun_ca_markov(maps, region)

    elapsed = time.time() - t0
    log(f"\nCOMPLETED in {elapsed/60:.1f} minutes")

    files = sorted(f for f in os.listdir(OUTPUT_DIR) if f.endswith('.json'))
    for f in files:
        sz = os.path.getsize(os.path.join(OUTPUT_DIR, f))
        log(f"  {f} ({sz:,} bytes)")

    summary = {
        'date': datetime.now().isoformat(),
        'time_min': round(elapsed / 60, 1),
        'stages_fixed': ['classification', 'change_detection', 'ecosystem_services',
                         'hotspot', 'gwr', 'ca_markov'],
        'fix': 'Float32 cast for band homogeneity, removed SoilGrids clay',
        'classification': {
            k: {'OA': v['overall_accuracy'], 'kappa': v['kappa']}
            for k, v in metrics.items()
        },
    }
    save_json(summary, 'analysis_summary.json')


if __name__ == '__main__':
    main()
