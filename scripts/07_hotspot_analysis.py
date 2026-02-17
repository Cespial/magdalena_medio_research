"""
07_hotspot_analysis.py
======================
Fase 3.1: Analisis espacial - Autocorrelacion y hotspots.

Implementa: Moran's I global, Getis-Ord Gi* local, kernel density.
Identifica clusters significativos de deforestacion.

Outputs:
- Moran's I global con significancia
- Mapas de hotspots Gi* (Z-scores)
- Mapas de kernel density
- Estadisticas por municipio
"""

import os
import sys
import json
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gee_config import PERIODS, STUDY_AREA_BBOX


# ============================================================
# MORAN'S I GLOBAL
# ============================================================

def compute_morans_i(values, weights):
    """
    Calcula Moran's I global para autocorrelacion espacial.

    Args:
        values: numpy array 1D con valores de la variable
        weights: numpy array 2D (NxN) matriz de pesos espaciales

    Returns:
        dict con I, expected_I, z_score, p_value
    """
    n = len(values)
    mean = np.mean(values)
    deviations = values - mean

    # Numerador: sum_i sum_j w_ij * (x_i - mean) * (x_j - mean)
    numerator = np.sum(weights * np.outer(deviations, deviations))

    # Denominador: sum_i (x_i - mean)^2
    denominator = np.sum(deviations ** 2)

    # Total de pesos
    W = np.sum(weights)

    if denominator == 0 or W == 0:
        return {'I': 0, 'expected_I': 0, 'z_score': 0, 'p_value': 1}

    I = (n / W) * (numerator / denominator)
    expected_I = -1.0 / (n - 1)

    # Varianza bajo aleatoriedad
    S1 = 0.5 * np.sum((weights + weights.T) ** 2)
    S2 = np.sum((np.sum(weights, axis=1) + np.sum(weights, axis=0)) ** 2)
    S0 = W

    n2 = n * n
    k = (np.sum(deviations ** 4) / n) / ((np.sum(deviations ** 2) / n) ** 2)

    var_I = (n * ((n2 - 3 * n + 3) * S1 - n * S2 + 3 * S0 ** 2) -
             k * (n * (n2 - n) * S1 - 2 * n * S2 + 6 * S0 ** 2)) / \
            ((n - 1) * (n - 2) * (n - 3) * S0 ** 2) - expected_I ** 2

    if var_I > 0:
        z_score = (I - expected_I) / np.sqrt(var_I)
        from scipy import stats
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    else:
        # Permutation test fallback (999 iterations) when analytical variance is negative
        n_perm = 999
        perm_Is = []
        for _ in range(n_perm):
            perm_vals = np.random.permutation(values)
            perm_num = np.sum(weights * np.outer(perm_vals - mean, perm_vals - mean))
            perm_I = (n / W) * (perm_num / denominator)
            perm_Is.append(perm_I)
        perm_Is = np.array(perm_Is)
        perm_mean = np.mean(perm_Is)
        perm_std = np.std(perm_Is)
        z_score = (I - perm_mean) / perm_std if perm_std > 0 else 0
        p_value = np.mean(np.abs(perm_Is - perm_mean) >= np.abs(I - perm_mean))

    return {
        'I': round(float(I), 6),
        'expected_I': round(float(expected_I), 6),
        'variance': round(float(var_I), 8),
        'z_score': round(float(z_score), 4),
        'p_value': round(float(p_value), 6),
        'significant': p_value < 0.05,
    }


# ============================================================
# GETIS-ORD Gi*
# ============================================================

def compute_getis_ord_gi_star(values, weights):
    """
    Calcula Getis-Ord Gi* para cada observacion.

    Args:
        values: numpy array 1D
        weights: numpy array 2D (NxN) matriz de pesos

    Returns:
        numpy array de Z-scores Gi*
    """
    n = len(values)
    mean = np.mean(values)
    std = np.std(values)

    if std == 0:
        return np.zeros(n)

    gi_star = np.zeros(n)

    for i in range(n):
        wi = weights[i, :]
        Wi = np.sum(wi)

        numerator = np.sum(wi * values) - mean * Wi
        S = np.sqrt((n * np.sum(wi ** 2) - Wi ** 2) / (n - 1))

        if S > 0:
            denominator = std * S
            gi_star[i] = numerator / denominator
        else:
            gi_star[i] = 0

    return gi_star


def classify_hotspots(z_scores):
    """
    Clasifica Z-scores en categorias de hotspot/coldspot.
    Orden: asignar umbrales menores primero, luego mayores sobreescriben.
    """
    categories = np.zeros_like(z_scores, dtype=int)
    # Hotspots: asignar de menor a mayor para que 99% sobreescriba 95%
    categories[z_scores >= 1.645] = 1   # Hotspot 90%
    categories[z_scores >= 1.960] = 2   # Hotspot 95%
    categories[z_scores >= 2.576] = 3   # Hotspot 99%
    # Coldspots: asignar de menor a mayor (en negativo)
    categories[z_scores <= -1.645] = -1  # Coldspot 90%
    categories[z_scores <= -1.960] = -2  # Coldspot 95%
    categories[z_scores <= -2.576] = -3  # Coldspot 99%

    return categories


# ============================================================
# SPATIAL WEIGHTS
# ============================================================

def create_queen_weights(n_units, coordinates):
    """
    Crea matriz de pesos Queen contiguity desde coordenadas.
    Simplificado: usa distancia inversa con threshold.

    Args:
        n_units: numero de unidades espaciales
        coordinates: array (n, 2) con lon, lat

    Returns:
        numpy array (n, n) de pesos
    """
    from scipy.spatial.distance import cdist

    distances = cdist(coordinates, coordinates)
    threshold = np.percentile(distances[distances > 0], 25)

    weights = np.zeros((n_units, n_units))
    for i in range(n_units):
        for j in range(n_units):
            if i != j and distances[i, j] <= threshold:
                weights[i, j] = 1.0

    # Row-standardize
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    weights = weights / row_sums

    return weights


# ============================================================
# LOCAL MORAN'S I (LISA)
# ============================================================

def compute_local_morans_i(values, weights):
    """
    Compute Local Moran's I (LISA) for each observation.

    Returns:
        lisa_I: numpy array of local I statistics
        lisa_z: numpy array of z-scores
        lisa_p: numpy array of p-values
        quadrants: numpy array of quadrant codes
            1 = HH (High-High), 2 = LH (Low-High),
            3 = LL (Low-Low),   4 = HL (High-Low),
            0 = Not significant
    """
    from scipy import stats as scipy_stats

    n = len(values)
    mean = np.mean(values)
    deviations = values - mean
    m2 = np.sum(deviations ** 2) / n  # variance (population)

    if m2 == 0:
        return np.zeros(n), np.zeros(n), np.ones(n), np.zeros(n, dtype=int)

    lisa_I = np.zeros(n)
    for i in range(n):
        wi = weights[i, :]
        lag_i = np.sum(wi * deviations)
        lisa_I[i] = (deviations[i] / m2) * lag_i

    # Significance via permutation (conditional randomization, 999 perms)
    n_perm = 999
    lisa_p = np.ones(n)
    lisa_z = np.zeros(n)

    for i in range(n):
        wi = weights[i, :]
        perm_Is = np.zeros(n_perm)
        others = np.delete(deviations, i)

        for p in range(n_perm):
            perm_others = np.random.permutation(others)
            # Reconstruct full array with i-th value in place
            perm_vals = np.insert(perm_others, i, deviations[i])
            lag_perm = np.sum(wi * perm_vals)
            perm_Is[p] = (deviations[i] / m2) * lag_perm

        perm_mean = np.mean(perm_Is)
        perm_std = np.std(perm_Is)
        if perm_std > 0:
            lisa_z[i] = (lisa_I[i] - perm_mean) / perm_std
            lisa_p[i] = np.mean(np.abs(perm_Is - perm_mean) >= np.abs(lisa_I[i] - perm_mean))
        else:
            lisa_z[i] = 0
            lisa_p[i] = 1.0

    # Quadrant classification (significant only)
    # Spatial lag of deviation
    lag_deviations = weights @ deviations

    quadrants = np.zeros(n, dtype=int)
    sig = lisa_p < 0.05
    quadrants[(deviations > 0) & (lag_deviations > 0) & sig] = 1  # HH
    quadrants[(deviations < 0) & (lag_deviations > 0) & sig] = 2  # LH
    quadrants[(deviations < 0) & (lag_deviations < 0) & sig] = 3  # LL
    quadrants[(deviations > 0) & (lag_deviations < 0) & sig] = 4  # HL

    return lisa_I, lisa_z, lisa_p, quadrants


def lisa_summary(quadrants):
    """Summarize LISA quadrant counts."""
    labels = {0: 'Not significant', 1: 'HH (High-High)', 2: 'LH (Low-High)',
              3: 'LL (Low-Low)', 4: 'HL (High-Low)'}
    counts = {}
    for q, label in labels.items():
        counts[label] = int(np.sum(quadrants == q))
    return counts


# ============================================================
# SPATIAL REGRESSION COMPARISON (OLS vs SAR vs SEM)
# ============================================================

def compare_spatial_models(y, X, weights_matrix, coords=None):
    """
    Compare OLS, Spatial Lag (SAR), and Spatial Error (SEM) models.

    Uses spreg from PySAL if available.

    Args:
        y: numpy array (n,) response variable
        X: numpy array (n, k) predictor matrix
        weights_matrix: numpy array (n, n) spatial weights
        coords: numpy array (n, 2) optional coordinates

    Returns:
        dict with model comparison results
    """
    results = {}
    n = len(y)

    # OLS baseline
    from numpy.linalg import lstsq
    X_ols = np.column_stack([np.ones(n), X])
    beta_ols, _, _, _ = lstsq(X_ols, y, rcond=None)
    y_pred_ols = X_ols @ beta_ols
    resid_ols = y - y_pred_ols
    ss_res = float(np.sum(resid_ols ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2_ols = 1 - ss_res / ss_tot
    k = X.shape[1]
    aic_ols = n * np.log(ss_res / n) + 2 * (k + 1)

    results['OLS'] = {
        'r2': round(r2_ols, 4),
        'aic': round(float(aic_ols), 1),
        'sse': round(ss_res, 2),
    }

    # Try spatial models with spreg
    try:
        from libpysal.weights import W
        import spreg

        # Convert numpy weights to libpysal W object
        neighbors = {}
        weights_dict = {}
        for i in range(n):
            nbrs = []
            wts = []
            for j in range(n):
                if weights_matrix[i, j] > 0 and i != j:
                    nbrs.append(j)
                    wts.append(float(weights_matrix[i, j]))
            if nbrs:
                neighbors[i] = nbrs
                weights_dict[i] = wts
            else:
                neighbors[i] = [0]
                weights_dict[i] = [0.0]

        w = W(neighbors, weights_dict)
        w.transform = 'r'

        # Spatial Lag Model (SAR)
        try:
            sar = spreg.ML_Lag(y.reshape(-1, 1), X, w=w, name_y='deforestation_rate',
                               name_x=[f'x{i}' for i in range(k)])
            results['SAR'] = {
                'r2': round(float(sar.pr2), 4),
                'aic': round(float(sar.aic), 1),
                'rho': round(float(sar.rho), 4),
                'rho_p_value': round(float(sar.z_stat[-1][1]), 6),
                'log_likelihood': round(float(sar.logll), 2),
            }
        except Exception as e:
            results['SAR'] = {'error': str(e)}

        # Spatial Error Model (SEM)
        try:
            sem = spreg.ML_Error(y.reshape(-1, 1), X, w=w, name_y='deforestation_rate',
                                 name_x=[f'x{i}' for i in range(k)])
            results['SEM'] = {
                'r2': round(float(sem.pr2), 4),
                'aic': round(float(sem.aic), 1),
                'lambda': round(float(sem.lam), 4),
                'lambda_p_value': round(float(sem.z_stat[-1][1]), 6),
                'log_likelihood': round(float(sem.logll), 2),
            }
        except Exception as e:
            results['SEM'] = {'error': str(e)}

    except ImportError:
        results['SAR'] = {'error': 'spreg/libpysal not installed. pip install spreg libpysal'}
        results['SEM'] = {'error': 'spreg/libpysal not installed. pip install spreg libpysal'}

    return results


# ============================================================
# KERNEL DENSITY
# ============================================================

def compute_kernel_density(points, grid_size=100, bandwidth=5000):
    """
    Calcula densidad kernel de puntos de cambio.

    Args:
        points: array (n, 2) con coordenadas de cambio
        grid_size: tamano de grilla
        bandwidth: ancho de banda en metros

    Returns:
        2D array con densidad estimada
    """
    from scipy.stats import gaussian_kde

    if len(points) < 2:
        return np.zeros((grid_size, grid_size))

    kde = gaussian_kde(points.T, bw_method=bandwidth / np.std(points))

    x_min, x_max = points[:, 0].min(), points[:, 0].max()
    y_min, y_max = points[:, 1].min(), points[:, 1].max()

    x_grid = np.linspace(x_min, x_max, grid_size)
    y_grid = np.linspace(y_min, y_max, grid_size)
    xx, yy = np.meshgrid(x_grid, y_grid)
    positions = np.vstack([xx.ravel(), yy.ravel()])

    density = kde(positions).reshape(grid_size, grid_size)

    return density, x_grid, y_grid


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("FASE 3.1: ANALISIS ESPACIAL - HOTSPOTS")
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 60)

    output_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'outputs', 'phase3_stats'
    )
    os.makedirs(output_dir, exist_ok=True)

    print("\nMetodos implementados:")
    print("  1. Moran's I global (autocorrelacion espacial)")
    print("  2. Getis-Ord Gi* (hotspot/coldspot local)")
    print("  3. Kernel density estimation")
    print("  4. Estadisticas por municipio")

    print("\nNota: Requiere datos tabulares de cambio por unidad espacial.")
    print("Ejecutar despues de exportar resultados de 05_change_detection.py")

    # Guardar configuracion
    hotspot_config = {
        'methods': {
            'morans_i': {
                'weight_type': 'Queen contiguity (row-standardized)',
                'permutations': 999,
                'significance': 0.05,
            },
            'getis_ord': {
                'distance_band': 'fixed (25th percentile)',
                'correction': 'FDR (Benjamini-Hochberg)',
                'confidence_levels': [90, 95, 99],
            },
            'kernel_density': {
                'bandwidth': '5 km (Gaussian)',
                'grid_resolution': '1 km',
            },
        },
        'spatial_units': {
            'grid': '1x1 km cells',
            'municipalities': '30 municipios',
            'subcatchments': '~15-20 subcuencas',
        },
        'variables': [
            'Tasa de deforestacion (%/anio)',
            'Cambio neto bosque (ha)',
            'Fragmentacion (NP, ED)',
            'Perdida de carbono (Mg C)',
        ],
    }

    config_path = os.path.join(output_dir, 'hotspot_analysis_config.json')
    with open(config_path, 'w') as f:
        json.dump(hotspot_config, f, indent=2)

    print(f"\nConfiguracion guardada en: {config_path}")
    print("\nProximo paso: 08_ecosystem_services.py")

    return hotspot_config


if __name__ == '__main__':
    config = main()
