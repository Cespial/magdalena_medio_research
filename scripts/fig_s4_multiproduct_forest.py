#!/usr/bin/env python3
"""
fig_s4_multiproduct_forest.py
=============================
Generates Supplementary Figure S4: Multi-product forest area time series
for the Magdalena Medio manuscript.

Compares Olofsson-adjusted forest area estimates (with 95% CI) against
MapBiomas Colombia Collection 1 and ESA WorldCover v200 (single-year).

Output:
    outputs/figures/fig_s4_multiproduct_forest.pdf
    outputs/figures/fig_s4_multiproduct_forest.png

Usage:
    python scripts/fig_s4_multiproduct_forest.py
"""

import os
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
FIG_DIR = os.path.join(PROJECT_DIR, "outputs", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

OUTPUT_BASE = os.path.join(FIG_DIR, "fig_s4_multiproduct_forest")

# ---------------------------------------------------------------------------
# Import project figure style
# ---------------------------------------------------------------------------
sys.path.insert(0, SCRIPT_DIR)
from figure_style import setup_journal_style, save_figure, DOUBLE_COL_WIDTH

plt = setup_journal_style()

# Override a few rcParams for this specific figure (slightly larger text
# since this is a single-panel supplementary figure at full page width)
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'lines.linewidth': 1.5,
    'lines.markersize': 7,
})

# ---------------------------------------------------------------------------
# DATA  (all values in thousands of hectares)
# ---------------------------------------------------------------------------

# Product 1: This study (Olofsson-adjusted) with 95% CI
study_years = np.array([2013, 2016, 2020, 2024])
study_forest = np.array([2021, 2029, 1868, 1591])  # x10^3 ha
study_ci = np.array([249, 246, 258, 262])           # +/- x10^3 ha
study_lo = study_forest - study_ci
study_hi = study_forest + study_ci

# Product 2: MapBiomas Colombia Collection 1
mb_years = np.array([2013, 2016, 2020, 2022])
mb_forest = np.array([1483, 1421, 1400, 1313])  # x10^3 ha

# Product 3: ESA WorldCover v200 (single epoch)
esa_year = 2021
esa_forest = 1790  # x10^3 ha

# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(DOUBLE_COL_WIDTH, 4.2))

# --- Shaded 95% CI band for this study ---
ax.fill_between(
    study_years, study_lo, study_hi,
    color='#7fbf7b', alpha=0.30, linewidth=0,
    label='95% CI (this study)',
)

# --- This study line ---
ax.plot(
    study_years, study_forest,
    color='#1b7837', linestyle='-', marker='o',
    markersize=7, markeredgecolor='white', markeredgewidth=0.8,
    zorder=5, label='This study (Olofsson-adjusted)',
)

# --- MapBiomas line ---
ax.plot(
    mb_years, mb_forest,
    color='#e6550d', linestyle='--', marker='^',
    markersize=7, markeredgecolor='white', markeredgewidth=0.8,
    zorder=5, label='MapBiomas Colombia Coll. 1',
)

# --- ESA WorldCover single point ---
ax.plot(
    esa_year, esa_forest,
    color='#3182bd', marker='D', markersize=8,
    markeredgecolor='white', markeredgewidth=0.8,
    linestyle='none', zorder=6,
    label='ESA WorldCover v200 (2021)',
)

# ---------------------------------------------------------------------------
# Annotations: data labels
# ---------------------------------------------------------------------------

# This study — labels above points
offsets_study = [(0, 12), (0, 12), (0, -18), (0, -18)]
for i, (x, y) in enumerate(zip(study_years, study_forest)):
    va = 'bottom' if offsets_study[i][1] > 0 else 'top'
    ax.annotate(
        f'{y:,}',
        xy=(x, y), xytext=offsets_study[i],
        textcoords='offset points', fontsize=8,
        color='#1b7837', ha='center', va=va,
        fontweight='bold',
    )

# MapBiomas — labels below points
offsets_mb = [(0, -16), (0, -16), (0, -16), (0, 10)]
for i, (x, y) in enumerate(zip(mb_years, mb_forest)):
    va = 'top' if offsets_mb[i][1] < 0 else 'bottom'
    ax.annotate(
        f'{y:,}',
        xy=(x, y), xytext=offsets_mb[i],
        textcoords='offset points', fontsize=8,
        color='#e6550d', ha='center', va=va,
    )

# ESA WorldCover
ax.annotate(
    f'{esa_forest:,}',
    xy=(esa_year, esa_forest), xytext=(10, 8),
    textcoords='offset points', fontsize=8,
    color='#3182bd', ha='left', va='bottom',
)

# ---------------------------------------------------------------------------
# Peace Agreement vertical line
# ---------------------------------------------------------------------------
ax.axvline(
    x=2016, color='#888888', linestyle=':', linewidth=1.0, zorder=1,
)
ax.text(
    2016.15, 2340, 'Peace Agreement\n(Nov 2016)',
    fontsize=8, color='#666666', ha='left', va='top',
    style='italic',
)

# ---------------------------------------------------------------------------
# Axes, grid, limits
# ---------------------------------------------------------------------------
ax.set_xlim(2012, 2025)
ax.set_ylim(1000, 2400)

ax.set_xlabel('Year')
ax.set_ylabel(r'Forest area ($\times$10$^{3}$ ha)')

# X-ticks: every year
ax.set_xticks(np.arange(2012, 2026, 1))
ax.set_xticklabels(
    [str(y) for y in range(2012, 2026)],
    rotation=0,
)

# Y-ticks
ax.set_yticks(np.arange(1000, 2500, 200))

# Horizontal gridlines only
ax.yaxis.grid(True, color='#d9d9d9', linewidth=0.5, linestyle='-')
ax.xaxis.grid(False)
ax.set_axisbelow(True)

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
leg = ax.legend(
    loc='upper right',
    frameon=True, framealpha=0.95,
    edgecolor='#cccccc',
    borderpad=0.6,
    handlelength=2.5,
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig.tight_layout()

# PDF (vector)
fig.savefig(OUTPUT_BASE + '.pdf', format='pdf')
print(f"  [OK] {OUTPUT_BASE}.pdf")

# PNG (raster, 300 dpi as requested)
fig.savefig(OUTPUT_BASE + '.png', format='png', dpi=300)
print(f"  [OK] {OUTPUT_BASE}.png")

plt.close(fig)
print("\nDone — fig_s4_multiproduct_forest generated.")
