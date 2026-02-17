"""
10b_gwr_diagnostics.py
======================
Fase 3.4b: GWR Bandwidth Sensitivity Analysis and Residual Spatial
Autocorrelation Diagnostics.

Implements two complementary analyses:
  1. Bandwidth sensitivity: run GWR at bandwidths 11, 25, 50, 100, 150, 200
     and report mean R2, median R2, AIC, ENP, ENP/n for each.
  2. Moran's I on OLS vs GWR residuals to quantify the reduction in spatial
     autocorrelation achieved by the local model.

Data strategy:
  - Reads existing results from gwr_drivers_results.json.
  - Checks for a local tabular file (gwr_sample_data.csv) first.
  - If not available, generates *synthetic* sample data whose global OLS
    statistics (R2, coefficient signs, sample size) are consistent with the
    published results.  This is clearly flagged in the output JSON
    ("data_source": "synthetic").

Does NOT require Google Earth Engine credentials.

Outputs:
  - outputs/phase3_stats/gwr_bandwidth_sensitivity.json
"""

import os
import sys
import json
import numpy as np
from datetime import datetime

# ---------------------------------------------------------------------------
# Import compute_gwr and fit_ols from the main GWR script
# ---------------------------------------------------------------------------
# We add the scripts directory to sys.path so we can import directly from the
# module file without triggering GEE initialisation that lives in gee_config.
# Because 10_gwr_drivers.py does `from gee_config import PERIODS, STUDY_AREA_BBOX`
# at module level, we monkey-patch the import so that it does not require an
# active GEE session.  We only need the pure-numpy functions.

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
_project_dir = os.path.dirname(_scripts_dir)

# Create a lightweight stub for gee_config so the import of 10_gwr_drivers
# does not require Earth Engine credentials.
import types
_gee_stub = types.ModuleType("gee_config")
_gee_stub.PERIODS = {}
_gee_stub.STUDY_AREA_BBOX = None
sys.modules["gee_config"] = _gee_stub

# Now we can safely import the helpers.
sys.path.insert(0, _scripts_dir)
from importlib import util as _importlib_util

_spec = _importlib_util.spec_from_file_location(
    "gwr_drivers", os.path.join(_scripts_dir, "10_gwr_drivers.py")
)
_gwr_mod = _importlib_util.module_from_spec(_spec)
_spec.loader.exec_module(_gwr_mod)

compute_gwr = _gwr_mod.compute_gwr
fit_ols = _gwr_mod.fit_ols


# ============================================================
# CONSTANTS
# ============================================================

BANDWIDTHS = [11, 25, 50, 100, 150, 200]

VARIABLE_NAMES = [
    "elevation", "slope", "dist_rivers", "dist_roads",
    "dist_urban", "pop_density", "precip_annual", "lst_mean",
]

# Bounding box for synthetic coordinate generation (lon, lat ranges for
# the Magdalena Medio study area).
LON_RANGE = (-75.0, -73.5)
LAT_RANGE = (6.0, 8.0)

RESULTS_JSON = os.path.join(
    _project_dir, "outputs", "phase3_stats", "gwr_drivers_results.json"
)
SAMPLE_CSV = os.path.join(
    _project_dir, "outputs", "phase3_stats", "gwr_sample_data.csv"
)
OUTPUT_JSON = os.path.join(
    _project_dir, "outputs", "phase3_stats", "gwr_bandwidth_sensitivity.json"
)


# ============================================================
# DATA LOADING / GENERATION
# ============================================================

def load_existing_results(path):
    """Load gwr_drivers_results.json and return the dict."""
    with open(path, "r") as f:
        return json.load(f)


def load_sample_csv(path):
    """
    Load tabular GWR sample data exported as CSV.

    Expected columns: lon, lat, deforestation_rate, elevation, slope,
    dist_rivers, dist_roads, dist_urban, pop_density, precip_annual, lst_mean.

    Returns
    -------
    X : ndarray (n, 8)
    y : ndarray (n,)
    coords : ndarray (n, 2)
    """
    import csv

    with open(path, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    n = len(rows)
    X = np.zeros((n, len(VARIABLE_NAMES)))
    y = np.zeros(n)
    coords = np.zeros((n, 2))

    for i, row in enumerate(rows):
        coords[i, 0] = float(row["lon"])
        coords[i, 1] = float(row["lat"])
        y[i] = float(row["deforestation_rate"])
        for j, v in enumerate(VARIABLE_NAMES):
            X[i, j] = float(row[v])

    return X, y, coords


def generate_synthetic_data(results, seed=42):
    """
    Generate synthetic sample data consistent with the known OLS results.

    The procedure is:
      1. Sample n = results['ols']['n'] points with coordinates uniformly
         distributed across the Magdalena Medio bounding box.
      2. Generate 8 standardised independent variables drawn from a
         multivariate normal with modest inter-variable correlations
         (informed by VIF values in the results).
      3. Construct y = X * beta_ols + noise, scaling noise so that
         R2 ~ results['ols']['r2'].

    The resulting data is *not* the real GEE-extracted data, but has
    matching dimensionality, coefficient structure, and global fit so that
    bandwidth sensitivity diagnostics are meaningful.

    Returns
    -------
    X : ndarray (n, 8)   -- standardised independent variables
    y : ndarray (n,)      -- synthetic deforestation rate
    coords : ndarray (n, 2) -- (lon, lat)
    """
    rng = np.random.RandomState(seed)
    n = results["ols"]["n"]  # 1470
    p = len(VARIABLE_NAMES)   # 8

    # ----- coordinates -----
    lons = rng.uniform(LON_RANGE[0], LON_RANGE[1], size=n)
    lats = rng.uniform(LAT_RANGE[0], LAT_RANGE[1], size=n)
    coords = np.column_stack([lons, lats])

    # ----- X matrix (standardised) -----
    # Build a correlation matrix loosely reflecting the VIF structure.
    # VIF_j = 1/(1 - R2_j);  high VIF => correlated with other variables.
    vif = results["vif"]
    r2_from_vif = {v: 1 - 1 / vif[v] for v in VARIABLE_NAMES}

    # Use the average pairwise |r| implied by each variable's R2_j.
    # For simplicity: generate independent standard normals, then mix
    # elevation and lst_mean (the two highest-VIF variables) with a
    # shared latent factor so their marginal R2_j is realistic.
    Z = rng.randn(n, p)  # independent base

    # Add spatial structure: a smooth trend across lat/lon so that
    # residual Moran's I will be informative.
    lat_norm = (lats - lats.mean()) / lats.std()
    lon_norm = (lons - lons.mean()) / lons.std()
    spatial_trend = 0.4 * lat_norm + 0.3 * lon_norm

    # Inject spatial trend into some variables
    Z[:, 0] += spatial_trend * 0.8   # elevation
    Z[:, 7] += spatial_trend * 0.7   # lst_mean
    Z[:, 2] += spatial_trend * 0.3   # dist_rivers
    Z[:, 5] += -spatial_trend * 0.2  # pop_density

    # Introduce cross-correlations between elevation and lst_mean (both
    # have highest VIFs: 6.17 and 4.95).
    shared_factor = rng.randn(n)
    Z[:, 0] = 0.6 * Z[:, 0] + 0.8 * shared_factor   # elevation
    Z[:, 7] = 0.65 * Z[:, 7] + 0.75 * shared_factor  # lst_mean

    # Standardise each column to mean 0, std 1
    X = (Z - Z.mean(axis=0)) / Z.std(axis=0)

    # ----- y vector -----
    # Use the OLS coefficients (excluding intercept) on the standardised X
    ols_coefs = results["ols"]["coefficients"]
    beta = np.array([ols_coefs[v] for v in VARIABLE_NAMES])

    signal = X @ beta
    target_r2 = results["ols"]["r2"]

    # var(signal) / (var(signal) + var(noise)) = target_r2
    # => var(noise) = var(signal) * (1 - r2) / r2
    var_signal = np.var(signal)
    if target_r2 > 0 and target_r2 < 1:
        var_noise = var_signal * (1 - target_r2) / target_r2
    else:
        var_noise = var_signal

    noise = rng.randn(n) * np.sqrt(var_noise)

    # Add mild spatial autocorrelation to noise so that Moran's I on OLS
    # residuals is positive (simulating omitted spatial variables).
    spatial_noise_component = 0.25 * (
        np.sin(2 * np.pi * lat_norm) * np.cos(np.pi * lon_norm)
    ) * np.sqrt(var_noise)
    noise += spatial_noise_component

    y = float(ols_coefs["intercept"]) + signal + noise

    return X, y, coords


# ============================================================
# MORAN'S I
# ============================================================

def spatial_weights_knn(coords, k=8):
    """
    Build a row-standardised k-nearest-neighbour spatial weights matrix.

    Parameters
    ----------
    coords : ndarray (n, 2)
    k : int  -- number of neighbours

    Returns
    -------
    W : ndarray (n, n) -- row-standardised weight matrix
    """
    from scipy.spatial.distance import cdist

    n = coords.shape[0]
    dists = cdist(coords, coords)

    W = np.zeros((n, n))
    for i in range(n):
        # indices of k nearest (exclude self)
        neighbours = np.argsort(dists[i])[1: k + 1]
        W[i, neighbours] = 1.0

    # Row-standardise
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W = W / row_sums
    return W


def morans_i(residuals, W):
    """
    Compute Moran's I statistic for spatial autocorrelation.

    Parameters
    ----------
    residuals : ndarray (n,)
    W : ndarray (n, n) -- row-standardised weight matrix

    Returns
    -------
    dict with keys I, E_I, Var_I, z_score, p_value
    """
    from scipy import stats as sp_stats

    n = len(residuals)
    z = residuals - residuals.mean()
    S0 = W.sum()

    numerator = n * float(z @ W @ z)
    denominator = S0 * float(z @ z)
    I = numerator / denominator if denominator != 0 else 0.0

    # Expected value under randomisation
    E_I = -1.0 / (n - 1)

    # Variance under randomisation (normal approximation)
    S1 = 0.5 * np.sum((W + W.T) ** 2)
    S2 = np.sum((W.sum(axis=1) + W.sum(axis=0)) ** 2)

    z2 = z ** 2
    b2 = (n * np.sum(z ** 4)) / (np.sum(z2) ** 2)  # kurtosis

    A = n * ((n ** 2 - 3 * n + 3) * S1 - n * S2 + 3 * S0 ** 2)
    B = b2 * ((n ** 2 - n) * S1 - 2 * n * S2 + 6 * S0 ** 2)
    C = (n - 1) * (n - 2) * (n - 3) * S0 ** 2

    Var_I = (A - B) / C - E_I ** 2 if C != 0 else 1e-10
    # Guard against negative variance from numerical issues
    Var_I = max(Var_I, 1e-12)

    z_score = (I - E_I) / np.sqrt(Var_I)
    p_value = 2.0 * (1.0 - sp_stats.norm.cdf(abs(z_score)))

    return {
        "I": round(float(I), 6),
        "E_I": round(float(E_I), 6),
        "Var_I": round(float(Var_I), 8),
        "z_score": round(float(z_score), 4),
        "p_value": round(float(p_value), 6),
    }


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 64)
    print("FASE 3.4b: GWR BANDWIDTH SENSITIVITY & MORAN'S I DIAGNOSTICS")
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 64)

    # Ensure output directory exists
    output_dir = os.path.dirname(OUTPUT_JSON)
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load existing results (mandatory)
    # ------------------------------------------------------------------
    print(f"\nLoading existing GWR results from:\n  {RESULTS_JSON}")
    results = load_existing_results(RESULTS_JSON)
    print(f"  OLS R2 = {results['ols']['r2']:.4f},  n = {results['ols']['n']}")
    print(f"  GWR mean R2 (bw=11) = {results['gwr']['mean_r2']:.4f}")

    # ------------------------------------------------------------------
    # 2. Load or generate sample data
    # ------------------------------------------------------------------
    if os.path.exists(SAMPLE_CSV):
        print(f"\nLocal sample data found:\n  {SAMPLE_CSV}")
        X, y, coords = load_sample_csv(SAMPLE_CSV)
        data_source = "local_csv"
        print(f"  Loaded {X.shape[0]} samples, {X.shape[1]} variables.")
    else:
        print(f"\nNo local sample CSV found at:\n  {SAMPLE_CSV}")
        print("  Generating synthetic data consistent with OLS results ...")
        X, y, coords = generate_synthetic_data(results, seed=42)
        data_source = "synthetic"
        print(f"  Generated {X.shape[0]} samples, {X.shape[1]} variables.")
        print("  NOTE: Results use synthetic data; re-run with real GEE")
        print("        export (gwr_sample_data.csv) for production figures.")

    # Quick OLS check on the generated/loaded data
    ols_check = fit_ols(X, y)
    print(f"\n  OLS sanity check on loaded data: R2 = {ols_check['r_squared']:.4f}")

    # ------------------------------------------------------------------
    # 3. Bandwidth Sensitivity Analysis
    # ------------------------------------------------------------------
    print("\n" + "-" * 50)
    print("BANDWIDTH SENSITIVITY ANALYSIS")
    print("-" * 50)
    print(f"  Bandwidths to evaluate: {BANDWIDTHS}")

    n = X.shape[0]
    sensitivity_rows = []

    for bw in BANDWIDTHS:
        # Ensure bandwidth does not exceed n-1
        bw_eff = min(bw, n - 1)
        print(f"\n  Running GWR with bandwidth = {bw_eff} ...")
        try:
            gwr_out = compute_gwr(X, y, coords, bandwidth=bw_eff, kernel="adaptive")
            row = {
                "bandwidth": bw_eff,
                "mean_r2": round(gwr_out["mean_r2"], 6),
                "median_r2": round(gwr_out["median_r2"], 6),
                "aic": round(gwr_out["aic"], 4),
                "enp": round(gwr_out["hat_matrix_trace"], 4),
                "enp_n_ratio": round(gwr_out["enp_n_ratio"], 6),
            }
            sensitivity_rows.append(row)
            print(f"    mean R2={row['mean_r2']:.4f}  median R2={row['median_r2']:.4f}"
                  f"  AIC={row['aic']:.1f}  ENP={row['enp']:.1f}"
                  f"  ENP/n={row['enp_n_ratio']:.4f}")
        except Exception as e:
            print(f"    ERROR at bw={bw_eff}: {e}")
            sensitivity_rows.append({
                "bandwidth": bw_eff,
                "mean_r2": None,
                "median_r2": None,
                "aic": None,
                "enp": None,
                "enp_n_ratio": None,
                "error": str(e),
            })

    # ------------------------------------------------------------------
    # 4. Moran's I: OLS residuals vs GWR residuals
    # ------------------------------------------------------------------
    print("\n" + "-" * 50)
    print("RESIDUAL SPATIAL AUTOCORRELATION (MORAN'S I)")
    print("-" * 50)

    # Use k=8 nearest neighbours for spatial weights
    k_neighbours = 8
    print(f"  Building spatial weights matrix (k={k_neighbours} nearest neighbours) ...")
    W = spatial_weights_knn(coords, k=k_neighbours)

    # OLS residuals
    ols_residuals = ols_check["residuals"]
    print("  Computing Moran's I on OLS residuals ...")
    moran_ols = morans_i(ols_residuals, W)
    print(f"    I = {moran_ols['I']:.4f},  z = {moran_ols['z_score']:.2f},"
          f"  p = {moran_ols['p_value']:.6f}")

    # GWR residuals (use the bandwidth from the original results: 11)
    original_bw = results["gwr"]["bandwidth"]
    bw_gwr = min(original_bw, n - 1)
    print(f"\n  Running GWR (bw={bw_gwr}) for residual comparison ...")
    gwr_for_moran = compute_gwr(X, y, coords, bandwidth=bw_gwr, kernel="adaptive")
    gwr_residuals = gwr_for_moran["local_residuals"]

    print("  Computing Moran's I on GWR residuals ...")
    moran_gwr = morans_i(gwr_residuals, W)
    print(f"    I = {moran_gwr['I']:.4f},  z = {moran_gwr['z_score']:.2f},"
          f"  p = {moran_gwr['p_value']:.6f}")

    # Reduction
    reduction_pct = (
        (abs(moran_ols["I"]) - abs(moran_gwr["I"])) / abs(moran_ols["I"]) * 100
        if abs(moran_ols["I"]) > 0 else 0.0
    )
    print(f"\n  Reduction in |Moran's I|: {reduction_pct:.1f}%")

    # ------------------------------------------------------------------
    # 5. Assemble and save output
    # ------------------------------------------------------------------
    output = {
        "metadata": {
            "script": "10b_gwr_diagnostics.py",
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data_source": data_source,
            "data_note": (
                "Synthetic data generated to match OLS R2, coefficients, "
                "and sample size from gwr_drivers_results.json. "
                "Replace with real GEE-exported gwr_sample_data.csv for "
                "production-quality diagnostics."
                if data_source == "synthetic"
                else "Local CSV sample data used."
            ),
            "n_samples": int(X.shape[0]),
            "n_variables": int(X.shape[1]),
            "variable_names": VARIABLE_NAMES,
            "ols_r2_check": round(ols_check["r_squared"], 6),
            "spatial_weights": f"k={k_neighbours} nearest neighbours, row-standardised",
        },
        "bandwidth_sensitivity": sensitivity_rows,
        "morans_i": {
            "ols_residuals": moran_ols,
            "gwr_residuals": moran_gwr,
            "gwr_bandwidth_used": bw_gwr,
            "reduction_in_abs_I_pct": round(reduction_pct, 2),
        },
    }

    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults written to:\n  {OUTPUT_JSON}")
    print("\n" + "=" * 64)
    print("DONE")
    print("=" * 64)

    return output


if __name__ == "__main__":
    main()
