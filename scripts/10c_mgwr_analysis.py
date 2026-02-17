#!/usr/bin/env python3
"""
10c_mgwr_analysis.py
=====================
Phase 1.5: MGWR Analysis with Real Exported Data.

Two-stage workflow:
  Stage 1 (GEE): Export the 1,470-point sample as CSV with real coordinates
                  and predictor values from Google Earth Engine.
  Stage 2 (Local): Run proper MGWR using the `mgwr` Python package
                    (variable-specific bandwidths), plus Random Forest
                    variable importance and partial dependence plots.

Replaces synthetic bandwidth sensitivity table S12 with real values.

Requirements:
    - Stage 1: earthengine-api, authenticated GEE session
    - Stage 2: mgwr, libpysal, spreg, scikit-learn, matplotlib, pandas, numpy

Usage:
    python scripts/10c_mgwr_analysis.py --stage 1   # Export from GEE
    python scripts/10c_mgwr_analysis.py --stage 2   # Run MGWR locally
    python scripts/10c_mgwr_analysis.py --stage all  # Both stages

Outputs:
    - outputs/phase3_stats/gwr_sample_data.csv       (Stage 1)
    - outputs/phase3_stats/mgwr_results.json          (Stage 2)
    - outputs/figures/fig_mgwr_bandwidths.pdf          (Stage 2)
    - outputs/figures/fig_rf_importance.pdf             (Stage 2)
"""

import argparse
import json
import math
import os
import sys
import warnings
from datetime import datetime

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'outputs', 'phase3_stats')
FIGURE_DIR = os.path.join(PROJECT_ROOT, 'outputs', 'figures')
CSV_PATH = os.path.join(OUTPUT_DIR, 'gwr_sample_data.csv')
MGWR_OUTPUT = os.path.join(OUTPUT_DIR, 'mgwr_results.json')

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Predictor variable names (must match GEE export)
# ---------------------------------------------------------------------------
PREDICTOR_NAMES = [
    'elevation', 'slope', 'dist_rivers', 'dist_roads',
    'dist_urban', 'pop_density', 'precip_annual', 'lst_mean',
]
RESPONSE_NAME = 'deforestation_rate'


# =========================================================================
# STAGE 1: Export sample data from GEE
# =========================================================================

def stage1_export_from_gee():
    """
    Export the 1,470-point stratified sample from GEE with coordinates
    and all 8 predictor variables plus the deforestation rate response.
    """
    print('=' * 65)
    print('STAGE 1: EXPORTING GWR SAMPLE DATA FROM GEE')
    print('=' * 65)

    import ee
    from gee_config import STUDY_AREA_BBOX, COLLECTIONS

    region = STUDY_AREA_BBOX

    # ------------------------------------------------------------------
    # 1. Load classified LULC maps for T2 (2016) and T4 (2024)
    # ------------------------------------------------------------------
    print('\n  Loading LULC maps for deforestation rate computation...')

    # Attempt to load from GEE assets; fall back to classification
    try:
        lulc_t2 = ee.Image(f'projects/{os.getenv("GEE_PROJECT_ID", "ee-maestria-tesis")}'
                           '/assets/lulc_2016').clip(region)
        lulc_t4 = ee.Image(f'projects/{os.getenv("GEE_PROJECT_ID", "ee-maestria-tesis")}'
                           '/assets/lulc_2024').clip(region)
        print('    Loaded from GEE assets.')
    except Exception:
        print('    GEE assets not found; computing from classification pipeline...')
        # Minimal fallback: generate from Hansen + proxy
        hansen = ee.Image(COLLECTIONS['hansen'])
        tc2000 = hansen.select('treecover2000')
        lossyear = hansen.select('lossyear')

        # Forest mask at T2 (2016): treecover > 25% AND no loss before 2017
        forest_t2 = tc2000.gte(25).And(lossyear.eq(0).Or(lossyear.gt(16)))
        forest_t4 = tc2000.gte(25).And(lossyear.eq(0).Or(lossyear.gt(24)))

        # Binary: 1=forest, 3=non-forest (matching class scheme)
        lulc_t2 = forest_t2.where(forest_t2.eq(1), 1).where(forest_t2.eq(0), 3)
        lulc_t4 = forest_t4.where(forest_t4.eq(1), 1).where(forest_t4.eq(0), 3)
        print('    Computed from Hansen GFC proxy.')

    # ------------------------------------------------------------------
    # 2. Compute deforestation rate at 1km grid
    # ------------------------------------------------------------------
    print('  Computing deforestation rate at 1km resolution...')

    # Forest fraction at each period (reproject to 1km)
    forest_frac_t2 = (lulc_t2.eq(1).Or(lulc_t2.eq(2))).reduceResolution(
        reducer=ee.Reducer.mean(), maxPixels=2048
    ).reproject(crs='EPSG:4326', scale=1000)

    forest_frac_t4 = (lulc_t4.eq(1).Or(lulc_t4.eq(2))).reduceResolution(
        reducer=ee.Reducer.mean(), maxPixels=2048
    ).reproject(crs='EPSG:4326', scale=1000)

    # Deforestation rate (%/yr) over post-agreement period (8 years: 2016-2024)
    defor_rate = forest_frac_t2.subtract(forest_frac_t4).divide(8).multiply(100).rename(
        RESPONSE_NAME
    )

    # ------------------------------------------------------------------
    # 3. Prepare predictor layers at 1km
    # ------------------------------------------------------------------
    print('  Preparing predictor layers at 1km...')

    srtm = ee.Image('USGS/SRTMGL1_003')
    elevation = srtm.select('elevation').reproject(crs='EPSG:4326', scale=1000).rename('elevation')
    slope = ee.Terrain.slope(srtm).reproject(crs='EPSG:4326', scale=1000).rename('slope')

    # Distance to rivers (JRC water occurrence > 50%)
    jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence')
    rivers = jrc.gte(50).selfMask()
    dist_rivers = rivers.fastDistanceTransform(1024).sqrt().multiply(30).divide(1000).reproject(
        crs='EPSG:4326', scale=1000
    ).rename('dist_rivers')  # km

    # Distance to roads (GHSL built-up proxy)
    ghsl = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S').mosaic().select('built_surface')
    roads_proxy = ghsl.gt(0).selfMask()
    dist_roads = roads_proxy.fastDistanceTransform(1024).sqrt().multiply(100).divide(1000).reproject(
        crs='EPSG:4326', scale=1000
    ).rename('dist_roads')  # km

    # Distance to urban centres (GHSL SMOD >= 30)
    ghsl_smod = ee.Image('JRC/GHSL/P2023A/GHS_SMOD_V2-0').select('smod_code')
    urban_centers = ghsl_smod.gte(30).selfMask()
    dist_urban = urban_centers.fastDistanceTransform(1024).sqrt().multiply(1000).divide(1000).reproject(
        crs='EPSG:4326', scale=1000
    ).rename('dist_urban')  # km

    # Population density (WorldPop 2020)
    pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(
        ee.Filter.eq('year', 2020)
    ).filter(
        ee.Filter.eq('country', 'COL')
    ).mosaic().reproject(crs='EPSG:4326', scale=1000).rename('pop_density')

    # Mean annual precipitation (CHIRPS 2012-2024)
    chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(
        '2012-01-01', '2024-12-31'
    ).filterBounds(region)
    precip = chirps.select('precipitation').sum().divide(13).reproject(
        crs='EPSG:4326', scale=1000
    ).rename('precip_annual')

    # Mean LST (MODIS 2012-2024)
    modis_lst = ee.ImageCollection('MODIS/061/MOD11A2').filterDate(
        '2012-01-01', '2024-12-31'
    ).filterBounds(region).select('LST_Day_1km')
    lst = modis_lst.mean().multiply(0.02).subtract(273.15).reproject(
        crs='EPSG:4326', scale=1000
    ).rename('lst_mean')

    # ------------------------------------------------------------------
    # 4. Stack all bands and sample
    # ------------------------------------------------------------------
    print('  Stacking bands and sampling...')

    stack = (defor_rate
             .addBands(elevation)
             .addBands(slope)
             .addBands(dist_rivers)
             .addBands(dist_roads)
             .addBands(dist_urban)
             .addBands(pop)
             .addBands(precip)
             .addBands(lst))

    # Stratified random sample: ~1,500 points within study area
    sample = stack.sample(
        region=region,
        scale=1000,
        numPixels=1500,
        seed=42,
        geometries=True,  # Include coordinates
    )

    # Add lat/lon as properties
    def add_coords(feature):
        coords = feature.geometry().coordinates()
        return feature.set('longitude', coords.get(0)).set('latitude', coords.get(1))

    sample = sample.map(add_coords)

    # ------------------------------------------------------------------
    # 5. Export to CSV
    # ------------------------------------------------------------------
    print('  Exporting sample to CSV...')

    # Get info (this triggers the computation)
    try:
        sample_info = sample.getInfo()
        features = sample_info['features']

        rows = []
        for f in features:
            props = f['properties']
            row = {
                'longitude': props.get('longitude'),
                'latitude': props.get('latitude'),
                RESPONSE_NAME: props.get(RESPONSE_NAME),
            }
            for pred in PREDICTOR_NAMES:
                row[pred] = props.get(pred)
            rows.append(row)

        df = pd.DataFrame(rows)
        # Drop rows with NaN
        df = df.dropna()
        df.to_csv(CSV_PATH, index=False)

        print(f'\n  Exported {len(df)} samples to: {CSV_PATH}')
        print(f'  Columns: {list(df.columns)}')
        print(f'  Response range: [{df[RESPONSE_NAME].min():.4f}, {df[RESPONSE_NAME].max():.4f}]')

    except Exception as e:
        print(f'\n  ERROR exporting: {e}')
        print('  Attempting Drive export as fallback...')

        from scripts.utils import export_table_to_drive
        task = export_table_to_drive(sample, 'gwr_sample_data', 'magdalena_medio')
        print(f'  Export task started: {task.status()}')
        print('  Download from Google Drive and place at:', CSV_PATH)
        return False

    return True


# =========================================================================
# STAGE 2: Run MGWR and Random Forest locally
# =========================================================================

def stage2_run_mgwr():
    """
    Run Multiscale GWR (MGWR) using the mgwr package, plus Random Forest
    variable importance and partial dependence as complements.
    """
    print('=' * 65)
    print('STAGE 2: MGWR + RANDOM FOREST ANALYSIS')
    print('=' * 65)

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    if not os.path.exists(CSV_PATH):
        print(f'\n  ERROR: Sample CSV not found at {CSV_PATH}')
        print('  Run stage 1 first: python scripts/10c_mgwr_analysis.py --stage 1')
        return False

    df = pd.read_csv(CSV_PATH)
    print(f'\n  Loaded {len(df)} samples from {CSV_PATH}')

    # Prepare arrays
    coords = df[['longitude', 'latitude']].values
    y = df[RESPONSE_NAME].values.reshape(-1, 1)
    X = df[PREDICTOR_NAMES].values

    n, k = X.shape
    print(f'  n = {n}, k = {k} predictors')
    print(f'  Response mean = {y.mean():.4f}, std = {y.std():.4f}')

    # Standardize predictors (important for GWR/MGWR)
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X)

    # ------------------------------------------------------------------
    # 2. OLS baseline
    # ------------------------------------------------------------------
    print('\n  Running OLS baseline...')
    from scipy import stats as scipy_stats

    # Add intercept
    X_ols = np.column_stack([np.ones(n), X_std])
    beta_ols = np.linalg.lstsq(X_ols, y, rcond=None)[0].flatten()
    y_pred_ols = X_ols @ beta_ols.reshape(-1, 1)
    resid_ols = y - y_pred_ols
    ss_res_ols = float(np.sum(resid_ols ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2_ols = 1 - ss_res_ols / ss_tot
    adj_r2_ols = 1 - (1 - r2_ols) * (n - 1) / (n - k - 1)
    aic_ols = n * np.log(ss_res_ols / n) + 2 * (k + 1)

    print(f'    R² = {r2_ols:.4f}, Adj R² = {adj_r2_ols:.4f}, AIC = {aic_ols:.1f}')

    # Compute t-stats and p-values
    mse_ols = ss_res_ols / (n - k - 1)
    var_beta = mse_ols * np.linalg.inv(X_ols.T @ X_ols).diagonal()
    se_beta = np.sqrt(var_beta)
    t_stats = beta_ols / se_beta
    p_values = 2 * scipy_stats.t.sf(np.abs(t_stats), df=n - k - 1)

    ols_results = {
        'r2': round(r2_ols, 4),
        'adj_r2': round(adj_r2_ols, 4),
        'aic': round(float(aic_ols), 1),
        'sse': round(ss_res_ols, 2),
        'n': n,
        'k': k,
        'coefficients': {},
    }
    var_names_ols = ['intercept'] + PREDICTOR_NAMES
    for i, name in enumerate(var_names_ols):
        ols_results['coefficients'][name] = {
            'beta': round(float(beta_ols[i]), 6),
            'se': round(float(se_beta[i]), 6),
            't': round(float(t_stats[i]), 4),
            'p': round(float(p_values[i]), 6),
            'significant': bool(p_values[i] < 0.05),
        }

    # ------------------------------------------------------------------
    # 3. GWR (standard, single bandwidth)
    # ------------------------------------------------------------------
    print('\n  Running standard GWR...')

    try:
        from mgwr.gwr import GWR, MGWR
        from mgwr.sel_bw import Sel_BW

        # Select bandwidth for GWR
        print('    Selecting GWR bandwidth (AICc)...')
        gwr_selector = Sel_BW(coords, y, X_std, fixed=False, kernel='bisquare')
        gwr_bw = gwr_selector.search(criterion='AICc')
        print(f'    Optimal GWR bandwidth: {gwr_bw}')

        # Fit GWR
        gwr_model = GWR(coords, y, X_std, bw=gwr_bw, fixed=False, kernel='bisquare')
        gwr_results = gwr_model.fit()

        r2_gwr = float(gwr_results.R2)
        adj_r2_gwr = float(gwr_results.adj_R2)
        aic_gwr = float(gwr_results.aicc)
        enp_gwr = float(gwr_results.ENP)

        # Local R² statistics
        local_r2 = gwr_results.localR2.flatten()
        mean_local_r2 = float(np.mean(local_r2))
        median_local_r2 = float(np.median(local_r2))

        print(f'    GWR R² = {r2_gwr:.4f}, AICc = {aic_gwr:.1f}')
        print(f'    Mean local R² = {mean_local_r2:.4f}, Median = {median_local_r2:.4f}')
        print(f'    ENP = {enp_gwr:.1f}, ENP/n = {enp_gwr/n:.3f}')

        gwr_dict = {
            'bandwidth': int(gwr_bw),
            'r2': round(r2_gwr, 4),
            'adj_r2': round(adj_r2_gwr, 4),
            'aicc': round(aic_gwr, 1),
            'enp': round(enp_gwr, 1),
            'enp_over_n': round(enp_gwr / n, 4),
            'mean_local_r2': round(mean_local_r2, 4),
            'median_local_r2': round(median_local_r2, 4),
            'local_r2_q25': round(float(np.percentile(local_r2, 25)), 4),
            'local_r2_q75': round(float(np.percentile(local_r2, 75)), 4),
            'sse': round(float(np.sum(gwr_results.resid_response ** 2)), 2),
        }

        # Coefficient summary per variable
        gwr_coefs = gwr_results.params  # n x (k+1) with intercept
        gwr_dict['coefficient_summary'] = {}
        coef_names = ['intercept'] + PREDICTOR_NAMES
        for j, name in enumerate(coef_names):
            col = gwr_coefs[:, j]
            gwr_dict['coefficient_summary'][name] = {
                'mean': round(float(np.mean(col)), 4),
                'median': round(float(np.median(col)), 4),
                'min': round(float(np.min(col)), 4),
                'max': round(float(np.max(col)), 4),
                'std': round(float(np.std(col)), 4),
            }

        gwr_success = True

    except ImportError:
        print('    WARNING: mgwr package not installed. Skipping GWR/MGWR.')
        print('    Install with: pip install mgwr')
        gwr_dict = {'error': 'mgwr package not installed'}
        gwr_success = False
    except Exception as e:
        print(f'    ERROR in GWR: {e}')
        gwr_dict = {'error': str(e)}
        gwr_success = False

    # ------------------------------------------------------------------
    # 4. MGWR (multiscale, variable-specific bandwidths)
    # ------------------------------------------------------------------
    mgwr_dict = {}
    if gwr_success:
        print('\n  Running MGWR (multiscale)...')
        try:
            mgwr_selector = Sel_BW(coords, y, X_std, fixed=False, kernel='bisquare',
                                    multi=True)
            mgwr_bw = mgwr_selector.search(criterion='AICc')
            print(f'    Variable-specific bandwidths: {mgwr_bw}')

            mgwr_model = MGWR(coords, y, X_std, selector=mgwr_selector,
                              fixed=False, kernel='bisquare')
            mgwr_results_fit = mgwr_model.fit()

            r2_mgwr = float(mgwr_results_fit.R2)
            adj_r2_mgwr = float(mgwr_results_fit.adj_R2)
            aic_mgwr = float(mgwr_results_fit.aicc)
            enp_mgwr = float(mgwr_results_fit.ENP)

            print(f'    MGWR R² = {r2_mgwr:.4f}, AICc = {aic_mgwr:.1f}')
            print(f'    ENP = {enp_mgwr:.1f}, ENP/n = {enp_mgwr/n:.3f}')

            mgwr_dict = {
                'variable_bandwidths': {},
                'r2': round(r2_mgwr, 4),
                'adj_r2': round(adj_r2_mgwr, 4),
                'aicc': round(aic_mgwr, 1),
                'enp': round(enp_mgwr, 1),
                'enp_over_n': round(enp_mgwr / n, 4),
                'sse': round(float(np.sum(mgwr_results_fit.resid_response ** 2)), 2),
            }

            # Variable-specific bandwidths
            bw_names = ['intercept'] + PREDICTOR_NAMES
            bw_array = list(mgwr_bw) if hasattr(mgwr_bw, '__iter__') else [mgwr_bw]
            for i, name in enumerate(bw_names):
                if i < len(bw_array):
                    mgwr_dict['variable_bandwidths'][name] = int(bw_array[i])

            # Coefficient summary
            mgwr_coefs = mgwr_results_fit.params
            mgwr_dict['coefficient_summary'] = {}
            for j, name in enumerate(bw_names):
                if j < mgwr_coefs.shape[1]:
                    col = mgwr_coefs[:, j]
                    mgwr_dict['coefficient_summary'][name] = {
                        'mean': round(float(np.mean(col)), 4),
                        'median': round(float(np.median(col)), 4),
                        'min': round(float(np.min(col)), 4),
                        'max': round(float(np.max(col)), 4),
                        'std': round(float(np.std(col)), 4),
                        'bandwidth': mgwr_dict['variable_bandwidths'].get(name),
                    }

            # Bandwidth sensitivity interpretation
            mgwr_dict['interpretation'] = {
                'small_bandwidth_vars': [],
                'large_bandwidth_vars': [],
            }
            median_bw = np.median(bw_array)
            for i, name in enumerate(bw_names):
                if i < len(bw_array):
                    if bw_array[i] < median_bw:
                        mgwr_dict['interpretation']['small_bandwidth_vars'].append(name)
                    else:
                        mgwr_dict['interpretation']['large_bandwidth_vars'].append(name)

        except Exception as e:
            print(f'    ERROR in MGWR: {e}')
            mgwr_dict = {'error': str(e)}

    # ------------------------------------------------------------------
    # 5. Leung F-test: OLS vs GWR
    # ------------------------------------------------------------------
    print('\n  Computing Leung F-test (OLS vs GWR)...')
    leung_test = {}
    if gwr_success:
        try:
            sse_ols = ss_res_ols
            sse_gwr = gwr_dict['sse']
            df_ols = n - k - 1
            df_gwr = n - enp_gwr

            if df_gwr > 0 and sse_gwr > 0:
                f_stat = ((sse_ols - sse_gwr) / (enp_gwr - k - 1)) / (sse_gwr / df_gwr)
                df1 = enp_gwr - k - 1
                df2 = df_gwr

                if df1 > 0 and df2 > 0:
                    p_value_f = float(scipy_stats.f.sf(f_stat, df1, df2))
                    leung_test = {
                        'f_statistic': round(float(f_stat), 4),
                        'df1': round(float(df1), 1),
                        'df2': round(float(df2), 1),
                        'p_value': round(p_value_f, 6),
                        'significant': bool(p_value_f < 0.05),
                        'interpretation': 'GWR significantly improves over OLS' if p_value_f < 0.05
                                          else 'No significant improvement',
                    }
                    print(f'    F = {f_stat:.4f}, p = {p_value_f:.6f} '
                          f'({"significant" if p_value_f < 0.05 else "not significant"})')
        except Exception as e:
            leung_test = {'error': str(e)}

    # ------------------------------------------------------------------
    # 6. Random Forest variable importance
    # ------------------------------------------------------------------
    print('\n  Running Random Forest variable importance...')

    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.model_selection import cross_val_score
        from sklearn.inspection import permutation_importance

        rf = RandomForestRegressor(
            n_estimators=500, max_depth=15, min_samples_leaf=5,
            random_state=42, n_jobs=-1
        )
        rf.fit(X_std, y.ravel())

        # Cross-validated R²
        cv_scores = cross_val_score(rf, X_std, y.ravel(), cv=5, scoring='r2')
        rf_r2_cv = float(np.mean(cv_scores))
        rf_r2_train = float(rf.score(X_std, y.ravel()))

        print(f'    RF R² (train) = {rf_r2_train:.4f}')
        print(f'    RF R² (5-fold CV) = {rf_r2_cv:.4f} ± {np.std(cv_scores):.4f}')

        # Feature importance (impurity-based)
        imp_mdi = rf.feature_importances_

        # Permutation importance (more reliable)
        perm_imp = permutation_importance(rf, X_std, y.ravel(), n_repeats=10,
                                          random_state=42, n_jobs=-1)
        imp_perm = perm_imp.importances_mean

        rf_dict = {
            'r2_train': round(rf_r2_train, 4),
            'r2_cv_mean': round(rf_r2_cv, 4),
            'r2_cv_std': round(float(np.std(cv_scores)), 4),
            'n_estimators': 500,
            'variable_importance': {},
        }

        for i, name in enumerate(PREDICTOR_NAMES):
            rf_dict['variable_importance'][name] = {
                'mdi_importance': round(float(imp_mdi[i]), 4),
                'permutation_importance': round(float(imp_perm[i]), 4),
                'permutation_std': round(float(perm_imp.importances_std[i]), 4),
                'rank_mdi': int(np.argsort(-imp_mdi).tolist().index(i) + 1),
                'rank_perm': int(np.argsort(-imp_perm).tolist().index(i) + 1),
            }

        # ------------------------------------------------------------------
        # 6b. Partial dependence (top 4 variables by permutation importance)
        # ------------------------------------------------------------------
        print('    Computing partial dependence for top 4 variables...')

        top4_idx = np.argsort(-imp_perm)[:4]
        top4_names = [PREDICTOR_NAMES[i] for i in top4_idx]
        rf_dict['partial_dependence_top4'] = top4_names

        try:
            from sklearn.inspection import partial_dependence
            pd_results = {}
            for idx in top4_idx:
                pd_result = partial_dependence(rf, X_std, features=[idx], grid_resolution=50)
                pd_results[PREDICTOR_NAMES[idx]] = {
                    'grid': pd_result['grid_values'][0].tolist(),
                    'average': pd_result['average'][0].tolist(),
                }
            rf_dict['partial_dependence'] = pd_results
        except Exception as e:
            rf_dict['partial_dependence'] = {'error': str(e)}

    except ImportError:
        print('    WARNING: scikit-learn not installed. Skipping RF analysis.')
        rf_dict = {'error': 'scikit-learn not installed'}
    except Exception as e:
        print(f'    ERROR in RF: {e}')
        rf_dict = {'error': str(e)}

    # ------------------------------------------------------------------
    # 7. Bandwidth sensitivity analysis (real values, replacing Table S12)
    # ------------------------------------------------------------------
    print('\n  Running bandwidth sensitivity analysis...')
    bandwidth_sensitivity = []

    if gwr_success:
        try:
            test_bandwidths = [int(gwr_bw)] + [25, 50, 100, 150, 200, 300, 500]
            test_bandwidths = sorted(set(test_bandwidths))

            for bw in test_bandwidths:
                if bw < k + 2:
                    continue
                try:
                    gwr_test = GWR(coords, y, X_std, bw=bw, fixed=False, kernel='bisquare')
                    gwr_test_fit = gwr_test.fit()

                    local_r2_test = gwr_test_fit.localR2.flatten()

                    entry = {
                        'bandwidth': bw,
                        'mean_r2': round(float(np.mean(local_r2_test)), 4),
                        'median_r2': round(float(np.median(local_r2_test)), 4),
                        'aicc': round(float(gwr_test_fit.aicc), 1),
                        'enp': round(float(gwr_test_fit.ENP), 1),
                        'enp_over_n': round(float(gwr_test_fit.ENP / n), 4),
                        'r2': round(float(gwr_test_fit.R2), 4),
                    }
                    bandwidth_sensitivity.append(entry)
                    print(f'    k={bw}: R²={entry["r2"]:.4f}, '
                          f'mean local R²={entry["mean_r2"]:.4f}, '
                          f'AICc={entry["aicc"]:.1f}, ENP={entry["enp"]:.1f}')
                except Exception:
                    pass

        except Exception as e:
            print(f'    ERROR in bandwidth sensitivity: {e}')

    # ------------------------------------------------------------------
    # 8. Save results
    # ------------------------------------------------------------------
    output = {
        'ols': ols_results,
        'gwr': gwr_dict,
        'mgwr': mgwr_dict,
        'leung_f_test': leung_test,
        'random_forest': rf_dict,
        'bandwidth_sensitivity': bandwidth_sensitivity,
        '_metadata': {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'n_samples': n,
            'n_predictors': k,
            'predictor_names': PREDICTOR_NAMES,
            'response': RESPONSE_NAME,
            'csv_path': CSV_PATH,
            'standardization': 'StandardScaler (zero mean, unit variance)',
            'gwr_kernel': 'adaptive bisquare',
            'mgwr_kernel': 'adaptive bisquare',
        },
    }

    with open(MGWR_OUTPUT, 'w') as f:
        json.dump(output, f, indent=2)

    print(f'\n  Results saved to: {MGWR_OUTPUT}')

    # ------------------------------------------------------------------
    # 9. Generate figures
    # ------------------------------------------------------------------
    print('\n  Generating figures...')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Figure: RF Variable Importance
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # MDI importance
        sorted_idx_mdi = np.argsort(imp_mdi)
        axes[0].barh([PREDICTOR_NAMES[i] for i in sorted_idx_mdi], imp_mdi[sorted_idx_mdi],
                     color='steelblue')
        axes[0].set_xlabel('MDI Importance')
        axes[0].set_title('(A) Mean Decrease Impurity')

        # Permutation importance
        sorted_idx_perm = np.argsort(imp_perm)
        axes[1].barh([PREDICTOR_NAMES[i] for i in sorted_idx_perm], imp_perm[sorted_idx_perm],
                     xerr=perm_imp.importances_std[sorted_idx_perm], color='coral')
        axes[1].set_xlabel('Permutation Importance')
        axes[1].set_title('(B) Permutation Importance')

        plt.tight_layout()
        fig_path = os.path.join(FIGURE_DIR, 'fig_rf_importance.pdf')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'    Saved: {fig_path}')

        # Figure: Bandwidth sensitivity
        if bandwidth_sensitivity:
            fig, ax1 = plt.subplots(figsize=(8, 5))
            bws = [e['bandwidth'] for e in bandwidth_sensitivity]
            r2s = [e['r2'] for e in bandwidth_sensitivity]
            mean_r2s = [e['mean_r2'] for e in bandwidth_sensitivity]

            ax1.plot(bws, r2s, 'o-', color='steelblue', label='Global R²')
            ax1.plot(bws, mean_r2s, 's--', color='coral', label='Mean local R²')
            ax1.axhline(y=r2_ols, color='gray', linestyle=':', label=f'OLS R² = {r2_ols:.3f}')
            ax1.set_xlabel('Bandwidth (k nearest neighbors)')
            ax1.set_ylabel('R²')
            ax1.set_title('GWR Bandwidth Sensitivity')
            ax1.legend()
            ax1.set_xscale('log')

            fig_path = os.path.join(FIGURE_DIR, 'fig_mgwr_bandwidths.pdf')
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f'    Saved: {fig_path}')

    except Exception as e:
        print(f'    Figure generation error: {e}')

    return True


# =========================================================================
# MAIN
# =========================================================================

def main():
    parser = argparse.ArgumentParser(description='MGWR Analysis')
    parser.add_argument('--stage', choices=['1', '2', 'all'], default='all',
                        help='Stage to run: 1=GEE export, 2=MGWR analysis, all=both')
    args = parser.parse_args()

    print('=' * 65)
    print('10c  MGWR ANALYSIS WITH REAL EXPORTED DATA')
    print(f'     {datetime.now().strftime("%Y-%m-%d %H:%M")}')
    print('=' * 65)

    if args.stage in ('1', 'all'):
        success = stage1_export_from_gee()
        if not success and args.stage == 'all':
            print('\nStage 1 failed. Cannot proceed to Stage 2.')
            return

    if args.stage in ('2', 'all'):
        stage2_run_mgwr()

    print('\nDone.')


if __name__ == '__main__':
    main()
