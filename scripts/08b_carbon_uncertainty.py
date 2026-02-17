#!/usr/bin/env python3
"""
08b_carbon_uncertainty.py
=========================
Phase 2 enhancement: Carbon stock estimation with propagated uncertainty.

Combines:
  - Olofsson et al. (2014) stratified area estimates (with area SE)
  - Tier 2 Colombia carbon pools (Alvarez et al. 2012) with pool-level SE

Uncertainty propagation follows the error-propagation formula for products
of independent random variables:
  Var(C_total) = sum_i [ c_i^2 * Var(A_i) + A_i^2 * Var(c_i)
                        + Var(A_i) * Var(c_i) ]

where:
  c_i   = total carbon density for class i (Mg C/ha)
  A_i   = Olofsson-adjusted area for class i (ha)
  Var(A_i) = (area_se_ha)^2
  Var(c_i) = sum of squared pool standard errors (c_above_se^2 + c_below_se^2
             + c_soil_se^2 + c_dead_se^2)

Net carbon change between periods propagates uncertainty as:
  Var(delta_C) = Var(C_t1) + Var(C_t2)

Outputs:
  outputs/phase3_stats/ecosystem_services_results.json  (overwritten)
"""

import json
import math
import os
import sys
from datetime import datetime

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

OLOFSSON_PATH = os.path.join(
    PROJECT_ROOT, 'outputs', 'phase3_stats', 'olofsson_area_estimates.json'
)
OUTPUT_PATH = os.path.join(
    PROJECT_ROOT, 'outputs', 'phase3_stats', 'ecosystem_services_results.json'
)

# ---------------------------------------------------------------------------
# Import Tier 2 carbon pools from gee_config
# We avoid initialising Earth Engine by importing only the dict we need.
# ---------------------------------------------------------------------------
# gee_config.py runs ee.Initialize() at import time, which would fail
# without credentials.  Instead we parse the CARBON_POOLS dict directly
# from the config file to keep this script self-contained and runnable
# without GEE credentials.
# ---------------------------------------------------------------------------

def _load_carbon_pools():
    """
    Parse CARBON_POOLS from gee_config.py without triggering ee.Initialize().
    Returns a dict identical in structure to gee_config.CARBON_POOLS.
    """
    config_path = os.path.join(PROJECT_ROOT, 'gee_config.py')
    with open(config_path, 'r') as f:
        source = f.read()

    # Extract the CARBON_POOLS block
    start_marker = 'CARBON_POOLS = {'
    start_idx = source.index(start_marker)
    # Find the matching closing brace
    brace_depth = 0
    end_idx = start_idx
    for i in range(start_idx, len(source)):
        if source[i] == '{':
            brace_depth += 1
        elif source[i] == '}':
            brace_depth -= 1
            if brace_depth == 0:
                end_idx = i + 1
                break

    block = source[start_idx:end_idx]
    # Remove inline comments (everything after # on each line)
    import re
    block_clean = re.sub(r'#[^\n]*', '', block)

    local_ns = {}
    exec(block_clean, {}, local_ns)  # noqa: S102
    return local_ns['CARBON_POOLS']


CARBON_POOLS = _load_carbon_pools()

# ---------------------------------------------------------------------------
# LULC class names (lightweight; avoids importing gee_config)
# ---------------------------------------------------------------------------
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
# Helper functions
# ---------------------------------------------------------------------------

def total_carbon_density(pools: dict) -> float:
    """Sum of all pool mean densities (Mg C/ha)."""
    return pools['c_above'] + pools['c_below'] + pools['c_soil'] + pools['c_dead']


def total_carbon_density_variance(pools: dict) -> float:
    """Variance of total carbon density = sum of squared pool SEs."""
    return (
        pools['c_above_se'] ** 2
        + pools['c_below_se'] ** 2
        + pools['c_soil_se'] ** 2
        + pools['c_dead_se'] ** 2
    )


def compute_period_carbon(period_data: dict) -> dict:
    """
    For a single period, compute total carbon stock and its uncertainty.

    Returns a dict with:
      - carbon_Mg_C: total carbon (Mg C)
      - carbon_se: standard error (Mg C)
      - carbon_ci95: half-width of 95 % CI (Mg C)
      - per_class: breakdown by class
    """
    per_class = period_data.get('per_class', {})

    total_carbon = 0.0
    total_variance = 0.0
    class_details = {}

    # Iterate over ALL 7 LULC classes (some may be absent from Olofsson)
    for cls_id in range(1, 8):
        cls_str = str(cls_id)
        pools = CARBON_POOLS[cls_id]
        c_i = total_carbon_density(pools)          # Mg C/ha
        var_c_i = total_carbon_density_variance(pools)  # (Mg C/ha)^2

        if cls_str in per_class:
            A_i = per_class[cls_str]['adjusted_area_ha']  # ha
            var_A_i = per_class[cls_str]['area_se_ha'] ** 2  # ha^2
        else:
            # Class not reported in Olofsson (e.g. Cultivos, Suelo desnudo
            # had zero mapped area and were merged).  Treat area = 0.
            A_i = 0.0
            var_A_i = 0.0

        # Carbon stock for this class
        C_i = c_i * A_i  # Mg C

        # Variance of C_i using error propagation for product of
        # two independent random variables:
        # Var(X*Y) = E[X]^2 Var(Y) + E[Y]^2 Var(X) + Var(X) Var(Y)
        var_C_i = (c_i ** 2) * var_A_i + (A_i ** 2) * var_c_i + var_A_i * var_c_i

        total_carbon += C_i
        total_variance += var_C_i

        class_details[cls_str] = {
            'name': CLASS_NAMES[cls_id],
            'c_density_MgC_ha': round(c_i, 2),
            'c_density_se_MgC_ha': round(math.sqrt(var_c_i), 2),
            'area_ha': round(A_i, 1),
            'area_se_ha': round(math.sqrt(var_A_i), 1),
            'carbon_Mg_C': round(C_i, 0),
            'carbon_se_Mg_C': round(math.sqrt(var_C_i), 0),
        }

    total_se = math.sqrt(total_variance)
    ci95 = 1.96 * total_se

    return {
        'carbon_Mg_C': round(total_carbon, 0),
        'carbon_se': round(total_se, 0),
        'carbon_ci95': round(ci95, 0),
        'carbon_lower_ci95': round(total_carbon - ci95, 0),
        'carbon_upper_ci95': round(total_carbon + ci95, 0),
        'per_class': class_details,
    }


def compute_carbon_change(result_t1: dict, result_t2: dict, label: str) -> dict:
    """
    Net carbon change (t2 - t1) with propagated uncertainty.
    Var(delta) = Var(C_t1) + Var(C_t2)  (independent periods).
    """
    net = result_t2['carbon_Mg_C'] - result_t1['carbon_Mg_C']
    var_net = result_t1['carbon_se'] ** 2 + result_t2['carbon_se'] ** 2
    se_net = math.sqrt(var_net)
    ci95 = 1.96 * se_net

    return {
        'label': label,
        'net_Mg_C': round(net, 0),
        'net_se': round(se_net, 0),
        'net_ci95': round(ci95, 0),
        'net_lower_ci95': round(net - ci95, 0),
        'net_upper_ci95': round(net + ci95, 0),
    }


def compute_carbon_change_correlated(period_data_t1: dict, period_data_t2: dict,
                                     label: str) -> dict:
    """
    Net carbon change with CORRELATED propagation.

    Key insight: the SAME Tier-2 carbon densities c_i are used for both
    periods.  In the difference ΔC = C_t2 - C_t1 = Σ c_i(A_i,t2 - A_i,t1),
    the c_i uncertainty partially cancels:

      Var(ΔC) = Σ_i [ c_i² × (Var(A_i,t1) + Var(A_i,t2))
                     + (A_i,t2 - A_i,t1)² × Var(c_i)
                     + (Var(A_i,t1) + Var(A_i,t2)) × Var(c_i) ]

    This is SMALLER than the independent formula because the dominant
    c_i² × Var(c_i) terms vanish when differencing.
    """
    per_class_t1 = period_data_t1.get('per_class', {})
    per_class_t2 = period_data_t2.get('per_class', {})

    net_carbon = 0.0
    total_variance = 0.0

    for cls_id in range(1, 8):
        cls_str = str(cls_id)
        pools = CARBON_POOLS[cls_id]
        c_i = total_carbon_density(pools)
        var_c_i = total_carbon_density_variance(pools)

        A_t1 = per_class_t1.get(cls_str, {}).get('adjusted_area_ha', 0.0)
        var_A_t1 = per_class_t1.get(cls_str, {}).get('area_se_ha', 0.0) ** 2
        A_t2 = per_class_t2.get(cls_str, {}).get('adjusted_area_ha', 0.0)
        var_A_t2 = per_class_t2.get(cls_str, {}).get('area_se_ha', 0.0) ** 2

        delta_A = A_t2 - A_t1
        delta_C_i = c_i * delta_A
        net_carbon += delta_C_i

        # Correlated variance: c_i terms partially cancel in the difference
        var_delta_i = (
            c_i ** 2 * (var_A_t1 + var_A_t2)       # area uncertainty (dominant)
            + delta_A ** 2 * var_c_i                  # density uncertainty on NET change only
            + (var_A_t1 + var_A_t2) * var_c_i         # cross term
        )
        total_variance += var_delta_i

    se_net = math.sqrt(total_variance)
    ci95 = 1.96 * se_net

    return {
        'label': label,
        'method': 'correlated (same c_i for both periods)',
        'net_Mg_C': round(net_carbon, 0),
        'net_se': round(se_net, 0),
        'net_ci95': round(ci95, 0),
        'net_lower_ci95': round(net_carbon - ci95, 0),
        'net_upper_ci95': round(net_carbon + ci95, 0),
    }


def monte_carlo_carbon(periods_data: dict, n_draws: int = 10000, seed: int = 42) -> dict:
    """
    Monte Carlo simulation for cumulative carbon change with correlated
    carbon densities (same c_i draw for all periods).

    Returns P(net loss), correlated credible intervals, and distribution stats.
    """
    np.random.seed(seed)

    period_keys = list(periods_data.keys())
    n_periods = len(period_keys)

    # Draw carbon densities ONCE per simulation (correlated across periods)
    carbon_draws = np.zeros((n_draws, n_periods))

    for draw in range(n_draws):
        # Same c_i for all periods in this draw
        c_draws = {}
        for cls_id in range(1, 8):
            pools = CARBON_POOLS[cls_id]
            c_i = total_carbon_density(pools)
            se_i = math.sqrt(total_carbon_density_variance(pools))
            c_draws[cls_id] = max(0, np.random.normal(c_i, se_i))

        # Independent area draws per period
        for p_idx, pk in enumerate(period_keys):
            per_class = periods_data[pk].get('per_class', {})
            total_c = 0.0
            for cls_str, cls_data in per_class.items():
                cls_id = int(cls_str)
                if cls_id not in c_draws:
                    continue
                area_mean = cls_data.get('adjusted_area_ha', 0.0)
                area_se = cls_data.get('area_se_ha', 0.0)
                area_draw = max(0, np.random.normal(area_mean, area_se))
                total_c += c_draws[cls_id] * area_draw
            carbon_draws[draw, p_idx] = total_c

    # Cumulative change: T4 - T1
    cumulative_change = carbon_draws[:, -1] - carbon_draws[:, 0]

    # Period-to-period changes
    period_changes = {}
    change_pairs = [
        (0, 1, f'{periods_data[period_keys[0]]["year"]}-{periods_data[period_keys[1]]["year"]}'),
        (1, 2, f'{periods_data[period_keys[1]]["year"]}-{periods_data[period_keys[2]]["year"]}'),
        (2, 3, f'{periods_data[period_keys[2]]["year"]}-{periods_data[period_keys[3]]["year"]}'),
    ]

    for i1, i2, label in change_pairs:
        delta = carbon_draws[:, i2] - carbon_draws[:, i1]
        period_changes[label] = {
            'mean_Mg_C': round(float(np.mean(delta)), 0),
            'se_Mg_C': round(float(np.std(delta)), 0),
            'ci95_lower': round(float(np.percentile(delta, 2.5)), 0),
            'ci95_upper': round(float(np.percentile(delta, 97.5)), 0),
            'p_loss': round(float(np.mean(delta < 0)), 4),
        }

    return {
        'n_draws': n_draws,
        'seed': seed,
        'method': 'Correlated Monte Carlo: same c_i draw for all periods, '
                  'independent area draws per period',
        'cumulative_change': {
            'mean_Mg_C': round(float(np.mean(cumulative_change)), 0),
            'se_Mg_C': round(float(np.std(cumulative_change)), 0),
            'ci95_lower': round(float(np.percentile(cumulative_change, 2.5)), 0),
            'ci95_upper': round(float(np.percentile(cumulative_change, 97.5)), 0),
            'p_net_loss': round(float(np.mean(cumulative_change < 0)), 4),
            'p_loss_gt_10pct': round(float(np.mean(
                cumulative_change < -0.10 * carbon_draws[:, 0]
            )), 4),
        },
        'period_changes': period_changes,
        'period_carbon_means': {
            pk: round(float(np.mean(carbon_draws[:, i])), 0)
            for i, pk in enumerate(period_keys)
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print('=' * 65)
    print('08b  CARBON STOCK ESTIMATION WITH PROPAGATED UNCERTAINTY')
    print(f'     {datetime.now().strftime("%Y-%m-%d %H:%M")}')
    print('=' * 65)

    # ------------------------------------------------------------------
    # 1. Load Olofsson area estimates
    # ------------------------------------------------------------------
    with open(OLOFSSON_PATH, 'r') as f:
        olofsson = json.load(f)

    periods_data = olofsson['periods']
    period_keys = list(periods_data.keys())

    print(f'\nLoaded Olofsson estimates for {len(period_keys)} periods')
    print(f'Carbon pool source: Tier 2 Colombia (Alvarez et al. 2012)')

    # Print carbon density table
    print('\nCarbon density table (Mg C/ha):')
    print(f'  {"Class":<22s} {"Above":>6s} {"Below":>6s} {"Soil":>6s} {"Dead":>6s} {"Total":>6s} {"SE_tot":>7s}')
    for cls_id in range(1, 8):
        p = CARBON_POOLS[cls_id]
        c_tot = total_carbon_density(p)
        se_tot = math.sqrt(total_carbon_density_variance(p))
        print(f'  {CLASS_NAMES[cls_id]:<22s} {p["c_above"]:6.1f} {p["c_below"]:6.1f} '
              f'{p["c_soil"]:6.1f} {p["c_dead"]:6.1f} {c_tot:6.1f} {se_tot:7.2f}')

    # ------------------------------------------------------------------
    # 2. Compute per-period carbon with uncertainty
    # ------------------------------------------------------------------
    period_results = {}

    for pk in period_keys:
        pd = periods_data[pk]
        result = compute_period_carbon(pd)
        year = pd['year']
        label = pd['label']

        period_results[pk] = result
        period_results[pk]['year'] = year

        # Preserve non-carbon ecosystem service fields from any existing
        # results (water_yield_mm, baseflow_mm, habitat_quality) if we
        # later want to merge.  For now we store only carbon fields.

        total_Tg = result['carbon_Mg_C'] / 1e6
        se_Tg = result['carbon_se'] / 1e6
        ci95_Tg = result['carbon_ci95'] / 1e6

        print(f'\n  {label} ({year})')
        print(f'    Total C: {total_Tg:,.2f} Tg C  (SE: {se_Tg:,.2f}, 95%CI: +/-{ci95_Tg:,.2f})')

        # Per-class summary
        for cls_str in sorted(result['per_class'].keys(), key=int):
            cd = result['per_class'][cls_str]
            if cd['carbon_Mg_C'] == 0:
                continue
            print(f'      {cd["name"]:<22s}  {cd["carbon_Mg_C"]/1e6:8.2f} Tg C  '
                  f'(SE: {cd["carbon_se_Mg_C"]/1e6:6.2f})')

    # ------------------------------------------------------------------
    # 3. Carbon change between consecutive periods
    # ------------------------------------------------------------------
    print('\n' + '-' * 65)
    print('  NET CARBON CHANGE (with uncertainty)')
    print('-' * 65)

    change_pairs = [
        ('pre_acuerdo', 'transicion', '2013-2016'),
        ('transicion', 'post_acuerdo_1', '2016-2020'),
        ('post_acuerdo_1', 'post_acuerdo_2', '2020-2024'),
    ]

    change_results = {}

    for pk1, pk2, label_years in change_pairs:
        change_key = f'carbon_change_{label_years.replace("-", "_")}'
        change = compute_carbon_change(period_results[pk1], period_results[pk2], label_years)
        change_results[change_key] = change

        net_Tg = change['net_Mg_C'] / 1e6
        se_Tg = change['net_se'] / 1e6
        ci95_Tg = change['net_ci95'] / 1e6
        sign = '+' if net_Tg >= 0 else ''

        print(f'  {label_years}: {sign}{net_Tg:,.2f} Tg C  '
              f'(SE: {se_Tg:,.2f}, 95%CI: [{(net_Tg - ci95_Tg):,.2f}, {(net_Tg + ci95_Tg):,.2f}])')

    # ------------------------------------------------------------------
    # 4. Correlated propagation for carbon CHANGE
    # ------------------------------------------------------------------
    print('\n' + '-' * 65)
    print('  CORRELATED PROPAGATION (same c_i for both periods)')
    print('-' * 65)

    correlated_change_results = {}

    for pk1, pk2, label_years in change_pairs:
        change_key = f'carbon_change_correlated_{label_years.replace("-", "_")}'
        change = compute_carbon_change_correlated(
            periods_data[pk1], periods_data[pk2], label_years
        )
        correlated_change_results[change_key] = change

        net_Tg = change['net_Mg_C'] / 1e6
        se_Tg = change['net_se'] / 1e6
        ci95_Tg = change['net_ci95'] / 1e6
        sign = '+' if net_Tg >= 0 else ''

        # Compare with independent
        indep_key = f'carbon_change_{label_years.replace("-", "_")}'
        indep_se = change_results[indep_key]['net_se'] / 1e6

        print(f'  {label_years}: {sign}{net_Tg:,.2f} Tg C  '
              f'(SE: {se_Tg:,.2f} vs independent: {indep_se:,.2f})  '
              f'95%CI: [{(net_Tg - ci95_Tg):,.2f}, {(net_Tg + ci95_Tg):,.2f}]')

    # Cumulative correlated
    cum_change = compute_carbon_change_correlated(
        periods_data[period_keys[0]], periods_data[period_keys[-1]],
        f'{periods_data[period_keys[0]]["year"]}-{periods_data[period_keys[-1]]["year"]}'
    )
    correlated_change_results['cumulative_correlated'] = cum_change
    cum_Tg = cum_change['net_Mg_C'] / 1e6
    cum_se_Tg = cum_change['net_se'] / 1e6
    cum_ci_Tg = cum_change['net_ci95'] / 1e6
    print(f'\n  Cumulative (correlated): {cum_Tg:,.2f} Tg C  '
          f'(SE: {cum_se_Tg:,.2f}, 95%CI: [{(cum_Tg - cum_ci_Tg):,.2f}, {(cum_Tg + cum_ci_Tg):,.2f}])')

    # ------------------------------------------------------------------
    # 5. Monte Carlo simulation (correlated)
    # ------------------------------------------------------------------
    print('\n' + '-' * 65)
    print('  MONTE CARLO SIMULATION (10,000 draws, correlated c_i)')
    print('-' * 65)

    mc_results = monte_carlo_carbon(periods_data, n_draws=10000)

    cum_mc = mc_results['cumulative_change']
    print(f'  Cumulative ΔC (MC): {cum_mc["mean_Mg_C"]/1e6:,.2f} Tg C  '
          f'(SE: {cum_mc["se_Mg_C"]/1e6:,.2f})')
    print(f'  MC 95% CI: [{cum_mc["ci95_lower"]/1e6:,.2f}, {cum_mc["ci95_upper"]/1e6:,.2f}] Tg C')
    print(f'  P(net loss): {cum_mc["p_net_loss"]:.4f}')
    print(f'  P(loss > 10% of baseline): {cum_mc["p_loss_gt_10pct"]:.4f}')

    for label, pc in mc_results['period_changes'].items():
        print(f'  {label}: mean={pc["mean_Mg_C"]/1e6:,.2f} Tg C, P(loss)={pc["p_loss"]:.4f}')

    # ------------------------------------------------------------------
    # 6. Build output JSON
    # ------------------------------------------------------------------
    # Read existing results to preserve non-carbon fields
    existing = {}
    if os.path.exists(OUTPUT_PATH):
        with open(OUTPUT_PATH, 'r') as f:
            existing = json.load(f)

    output = {}

    for pk in period_keys:
        pr = period_results[pk]
        # Merge with existing ecosystem service fields if present
        entry = {}
        if pk in existing:
            entry.update(existing[pk])

        entry['year'] = pr['year']
        entry['carbon_Mg_C'] = pr['carbon_Mg_C']
        entry['carbon_se'] = pr['carbon_se']
        entry['carbon_ci95'] = pr['carbon_ci95']
        entry['carbon_lower_ci95'] = pr['carbon_lower_ci95']
        entry['carbon_upper_ci95'] = pr['carbon_upper_ci95']
        entry['carbon_per_class'] = pr['per_class']

        output[pk] = entry

    # Independent propagation changes
    for ck, cv in change_results.items():
        output[ck] = cv

    # Correlated propagation changes
    for ck, cv in correlated_change_results.items():
        output[ck] = cv

    # Monte Carlo results
    output['monte_carlo'] = mc_results

    # Metadata
    output['_metadata'] = {
        'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        'carbon_pool_source': 'Tier 2 Colombia (Alvarez et al. 2012)',
        'area_source': 'Olofsson et al. (2014) stratified estimators',
        'uncertainty_methods': {
            'independent': 'Var(C_i) = c_i^2*Var(A_i) + A_i^2*Var(c_i) + Var(A_i)*Var(c_i); '
                           'Var(delta) = Var(C_t1) + Var(C_t2)',
            'correlated': 'Var(delta) = sum_i [c_i^2*(Var(A_t1)+Var(A_t2)) + '
                          'dA_i^2*Var(c_i) + (Var(A_t1)+Var(A_t2))*Var(c_i)]',
            'monte_carlo': 'Same c_i draw for all periods, independent area draws, '
                           '10000 iterations',
        },
        'ci_level': '95% (z=1.96 for analytical, 2.5/97.5 percentiles for MC)',
    }

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f'\nResults written to: {OUTPUT_PATH}')
    print('Done.')


if __name__ == '__main__':
    main()
