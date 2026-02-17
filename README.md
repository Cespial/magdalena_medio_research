# Post-conflict LULCC in Colombia's Magdalena Medio (2012-2024)

Multi-temporal land use/land cover change analysis using Google Earth Engine, Random Forest classification, and spatial statistical methods.

## Citation

Espinal, C. (2026). Post-conflict land use transitions and ecosystem service loss in Colombia's Magdalena Medio: A multi-temporal remote sensing analysis (2012-2024). *[Manuscript in preparation]*.

## Project structure

```
magdalena_medio_research/
├── gee_config.py                  # Central configuration (study area, periods, classes, params)
├── requirements.txt               # Python dependencies
├── scripts/                       # Google Earth Engine analysis scripts
│   ├── 01_preprocessing.py        # Image compositing, cloud masking, band harmonization
│   ├── 02_training_samples.py     # Reference data generation (Hansen GFC + JRC + GHSL)
│   ├── 03_classification.py       # Random Forest classification (7 classes, 4 periods)
│   ├── 04_accuracy_assessment.py  # Error matrices, OA, Kappa, PA/UA
│   ├── 05_change_detection.py     # Transition matrices, annual rates
│   ├── 06_fragmentation.py        # Landscape fragmentation metrics
│   ├── 07_hotspot_analysis.py     # Moran's I, Getis-Ord Gi*, kernel density
│   ├── 08_ecosystem_services.py   # Carbon (IPCC Tier 1), water yield, habitat quality
│   ├── 09_climate_analysis.py     # CHIRPS precip, MODIS LST, Mann-Kendall, SPI
│   ├── 10_gwr_drivers.py         # OLS + GWR driver analysis (8 variables)
│   ├── 11_ca_markov.py           # CA-Markov scenarios (BAU, Conservation, PDET)
│   ├── 12_visualization.py       # Figure generation templates
│   └── utils.py                  # Shared utility functions
├── run_analysis.py               # Main execution pipeline
├── run_phase4_figures.py         # Figure and table generation (11 figs, 6 tables)
├── run_phase6_qc.py             # Quality control validation (247 checks)
├── outputs/
│   ├── phase3_stats/             # JSON results from all analyses
│   │   ├── classification_metrics.json
│   │   ├── change_detection_results.json
│   │   ├── ecosystem_services_results.json
│   │   ├── climate_analysis_results.json
│   │   ├── hotspot_analysis_results.json
│   │   ├── gwr_drivers_results.json
│   │   ├── feature_importance.json
│   │   └── ca_markov_results.json
│   ├── figures/                  # Publication-ready figures (300 DPI PNG)
│   ├── tables/                   # CSV tables for manuscript
│   └── phase6_qc/              # QC reports
└── paper/
    ├── manuscript_v2.md          # Full manuscript (43 references, ~7,500 words)
    ├── supplementary_materials.md # Tables S1-S11, Figure S1
    ├── cover_letter.md
    └── submission_checklist.md
```

## Requirements

- Python 3.11+
- Google Earth Engine account with project access
- Conda environment recommended

### Setup

```bash
# Create conda environment
conda create -n magdalena_medio python=3.11
conda activate magdalena_medio

# Install dependencies
pip install -r requirements.txt

# Configure GEE authentication
earthengine authenticate
```

### Environment variables

Create a `.env` file in the project root:

```
GEE_PROJECT_ID=ee-maestria-tesis
```

## Reproduction

### Step 1: Run GEE analyses (Phase 3)

The analysis scripts in `scripts/` are designed to run sequentially through Google Earth Engine. The main pipeline is:

```bash
python run_analysis.py
```

This executes scripts 01-11 and generates JSON outputs in `outputs/phase3_stats/`. Computation time: ~25 minutes (depends on GEE server load).

**Key parameters** (configured in `gee_config.py`):
- Study area: Magdalena Medio bounding box [-75.0, 6.0, -73.5, 8.0]
- Periods: T1 (2013), T2 (2016), T3 (2020), T4 (2024)
- Classification: Random Forest, 200 trees, 12 features, tileScale=4
- Hansen GFC: v1.12 (UMD/hansen/global_forest_change_2024_v1_12)

### Step 2: Generate figures and tables (Phase 4)

```bash
python run_phase4_figures.py
```

Generates 11 figures (300 DPI PNG) and 6 CSV tables from JSON outputs.

### Step 3: Quality control (Phase 6)

```bash
python run_phase6_qc.py
```

Validates 247 cross-consistency checks across all outputs, figures, tables, and manuscript.

## Data sources

| Dataset | Source | Resolution | Access |
|---------|--------|-----------|--------|
| Landsat 8/9 C2 L2 SR | USGS via GEE | 30 m | `LANDSAT/LC08/C02/T1_L2` |
| Sentinel-2 SR Harmonized | ESA via GEE | 10-20 m | `COPERNICUS/S2_SR_HARMONIZED` |
| SRTM Elevation | USGS via GEE | 30 m | `USGS/SRTMGL1_003` |
| Hansen GFC v1.12 | UMD via GEE | 30 m | `UMD/hansen/global_forest_change_2024_v1_12` |
| CHIRPS Daily v2.0 | UCSB via GEE | ~5.5 km | `UCSB-CHG/CHIRPS/DAILY` |
| MODIS LST (MOD11A2) | NASA via GEE | 1 km | `MODIS/061/MOD11A2` |
| WorldPop | WorldPop via GEE | 100 m | `WorldPop/GP/100m/pop` |
| JRC Global Surface Water | JRC via GEE | 30 m | `JRC/GSW1_4/GlobalSurfaceWater` |
| GHSL Settlement Model | JRC via GEE | 1 km | `JRC/GHSL/P2023A/GHS_SMOD` |

## Key results

- **Forest loss**: Total forest cover declined 43.9% (2,599,275 to 1,457,337 ha) between 2013-2024
- **Dominant transition**: Dense forest to pastures (986,100 ha in 2020-2024 alone)
- **Carbon loss**: Net 190 Mt C (-33.5%), equivalent to ~697 Mt CO2
- **Deforestation hotspots**: 454 significant clusters at 99% confidence (Getis-Ord Gi*)
- **Driver analysis**: GWR R^2 = 0.609 vs OLS R^2 = 0.143; elevation, LST, distance to rivers as key drivers
- **Scenarios**: Conservation pathway projects 72.5% dense forest by 2040 vs 60.9% under BAU

## Classification scheme

| ID | Class | Description |
|----|-------|-------------|
| 1 | Dense forest | Canopy cover >60% |
| 2 | Secondary forest | Canopy cover 30-60%, successional |
| 3 | Pastures/grasslands | Cattle ranching, natural and introduced grasses |
| 4 | Croplands | Oil palm, rice, coca, other crops |
| 5 | Water bodies | Rivers, wetlands, reservoirs |
| 6 | Urban/built-up | Settlements, infrastructure |
| 7 | Bare soil/mining | Exposed areas, artisanal mining |

## License

Code: MIT License
Data: All input data are publicly available through Google Earth Engine.

## Contact

Cristian Espinal
Camara de Comercio de Medellin, Medellin, Colombia
