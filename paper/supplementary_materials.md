# Supplementary Materials

## Post-conflict land use transitions and ecosystem service loss in Colombia's Magdalena Medio: A multi-temporal remote sensing analysis (2012-2024)

Cristian Espinal

---

## Table S1. Complete list of municipalities in the study area

| # | Municipality | Department | Area (km2) | Population (est.) | PDET | Economic activities |
|---|-------------|------------|-----------|-------------------|------|---------------------|
| 1 | Barrancabermeja | Santander | 1,154 | 191,000 | Yes | Petroleum, services |
| 2 | Puerto Wilches | Santander | 1,529 | 32,000 | Yes | Oil palm, fisheries |
| 3 | Sabana de Torres | Santander | 1,386 | 20,000 | Yes | Cattle, oil palm |
| 4 | Cimitarra | Santander | 3,189 | 42,000 | Yes | Cattle, mining |
| 5 | Puerto Parra | Santander | 671 | 7,000 | Yes | Cattle, agriculture |
| 6 | Rionegro | Santander | 1,244 | 28,000 | Yes | Cattle, agriculture |
| 7 | San Vicente de Chucuri | Santander | 1,151 | 34,000 | Yes | Cacao, cattle |
| 8 | El Carmen de Chucuri | Santander | 905 | 19,000 | Yes | Agriculture, cattle |
| 9 | Betulia | Santander | 577 | 5,000 | Yes | Agriculture |
| 10 | Bolivar | Santander | 1,349 | 13,000 | Yes | Agriculture, cattle |
| 11 | Landazuri | Santander | 651 | 16,000 | Yes | Cacao, agriculture |
| 12 | Simacota | Santander | 1,093 | 9,000 | Yes | Agriculture |
| 13 | Santa Helena del Opon | Santander | 740 | 4,000 | Yes | Agriculture |
| 14 | El Penon | Santander | 381 | 5,500 | Yes | Agriculture |
| 15 | Puerto Berrio | Antioquia | 1,184 | 48,000 | Yes | Cattle, mining |
| 16 | Yondo | Antioquia | 1,881 | 19,000 | Yes | Petroleum, cattle, coca |
| 17 | Puerto Nare | Antioquia | 554 | 18,000 | Yes | Cattle, fisheries |
| 18 | Puerto Triunfo | Antioquia | 365 | 20,000 | Yes | Tourism, cattle |
| 19 | Caracoli | Antioquia | 307 | 4,500 | Yes | Agriculture |
| 20 | Maceo | Antioquia | 451 | 7,000 | Yes | Agriculture, mining |
| 21 | San Pablo | Bolivar | 1,970 | 33,000 | Yes | Mining, coca, cattle |
| 22 | Simiti | Bolivar | 1,241 | 19,000 | Yes | Mining, cattle |
| 23 | Santa Rosa del Sur | Bolivar | 2,816 | 42,000 | Yes | Mining, cattle, coca |
| 24 | Cantagallo | Bolivar | 886 | 8,500 | Yes | Petroleum, cattle |
| 25 | Morales | Bolivar | 1,351 | 21,000 | Yes | Agriculture, cattle |
| 26 | Arenal | Bolivar | 465 | 18,000 | Yes | Agriculture, fisheries |
| 27 | Aguachica | Cesar | 876 | 97,000 | Yes | Commerce, agriculture |
| 28 | Gamarra | Cesar | 336 | 16,000 | Yes | Fisheries, agriculture |
| 29 | San Martin | Cesar | 853 | 18,000 | Yes | Cattle, oil palm |
| 30 | San Alberto | Cesar | 608 | 25,000 | Yes | Oil palm, cattle |

---

## Table S2. Google Earth Engine collections and parameters

| Dataset | GEE Asset ID | Resolution | Variables used |
|---------|-------------|-----------|----------------|
| Landsat 8 SR C2 | LANDSAT/LC08/C02/T1_L2 | 30 m | SR_B2-B7, QA_PIXEL |
| Landsat 9 SR C2 | LANDSAT/LC09/C02/T1_L2 | 30 m | SR_B2-B7, QA_PIXEL |
| Sentinel-2 SR | COPERNICUS/S2_SR_HARMONIZED | 10-20 m | B2-B12, SCL |
| MODIS LST | MODIS/061/MOD11A2 | 1 km | LST_Day_1km |
| CHIRPS Daily | UCSB-CHG/CHIRPS/DAILY | ~5.5 km | precipitation |
| SRTM DEM | USGS/SRTMGL1_003 | 30 m | elevation |
| Hansen GFC v1.12 | UMD/hansen/global_forest_change_2024_v1_12 | 30 m | treecover2000, lossyear, gain |
| JRC Water | JRC/GSW1_4/GlobalSurfaceWater | 30 m | occurrence |
| WorldPop | WorldPop/GP/100m/pop | 100 m | population |
| GHSL SMOD | JRC/GHSL/P2023A/GHS_SMOD | 1 km | smod_code |

Note: Hansen GFC v1.12 (through 2024) was used instead of v1.11 (through 2023) to enable accurate forest classification for the T4 (2024) period. The lossyear band encoding in v1.12 includes year 24 (2024).

---

## Table S3. Random Forest classification parameters

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Algorithm | ee.Classifier.smileRandomForest | Standard GEE implementation |
| Number of trees (ntree) | 200 | Balanced accuracy vs. computation |
| Min leaf population | 5 | Prevent overfitting |
| Bag fraction | 0.632 | Default bootstrap proportion |
| Seed | 42 | Reproducibility |
| Training/validation split | 70/30 | Stratified random |
| Samples per class | ~300 (varies by class availability) | Balanced stratified |
| Total samples per period | ~1,000-1,100 | Sum of stratified samples |
| tileScale | 4 | Manage GEE memory for large areas |
| bestEffort | True | Automatic scale adjustment |
| Image casting | Float32 | Ensure homogeneous collections |
| Output type | Int8 | Reduce computation chain complexity |

---

## Table S4. Spectral features used in classification (12 total)

| # | Feature | Type | Formula / Source |
|---|---------|------|-----------------|
| 1 | blue | Reflectance | SR_B2 (Landsat) / B2 (S2) |
| 2 | green | Reflectance | SR_B3 (Landsat) / B3 (S2) |
| 3 | red | Reflectance | SR_B4 (Landsat) / B4 (S2) |
| 4 | nir | Reflectance | SR_B5 (Landsat) / B8 (S2) |
| 5 | swir1 | Reflectance | SR_B6 (Landsat) / B11 (S2) |
| 6 | swir2 | Reflectance | SR_B7 (Landsat) / B12 (S2) |
| 7 | NDVI | Index | (NIR - Red) / (NIR + Red) |
| 8 | NDWI | Index | (Green - NIR) / (Green + NIR) |
| 9 | NDBI | Index | (SWIR1 - NIR) / (SWIR1 + NIR) |
| 10 | NBR | Index | (NIR - SWIR2) / (NIR + SWIR2) |
| 11 | elevation | Topographic | SRTM 30 m |
| 12 | slope | Topographic | ee.Terrain.slope(SRTM) |

Note: Surface reflectance scaling applied before index calculation: Landsat multiply(0.0000275).add(-0.2); Sentinel-2 multiply(0.0001).

---

## Table S5. Carbon pool values by LULC class (Mg C ha^-1)

| LULC Class | C_above | C_below | C_soil | C_dead | Total | Source |
|-----------|---------|---------|--------|--------|-------|--------|
| Dense forest | 120 | 30 | 80 | 12 | 242 | IPCC 2006, tropical humid |
| Secondary forest | 60 | 15 | 65 | 6 | 146 | IPCC 2006, secondary |
| Pastures | 5 | 3 | 40 | 0.5 | 48.5 | IPCC 2006, grassland |
| Croplands | 8 | 2 | 35 | 0.5 | 45.5 | IPCC 2006, annual crop |
| Water | 0 | 0 | 0 | 0 | 0 | -- |
| Urban | 2 | 0 | 20 | 0 | 22 | IPCC 2006, settlements |
| Bare soil | 0 | 0 | 15 | 0 | 15 | IPCC 2006, other land |

Carbon stock change per hectare for key transitions:
- Dense forest -> Pastures: -193.5 Mg C ha^-1
- Dense forest -> Croplands: -196.5 Mg C ha^-1
- Dense forest -> Secondary forest: -96 Mg C ha^-1
- Secondary forest -> Pastures: -97.5 Mg C ha^-1
- Pastures -> Secondary forest: +97.5 Mg C ha^-1 (recovery)

---

## Table S6. Habitat quality model parameters

| Parameter | Value |
|-----------|-------|
| **Habitat suitability** | |
| Dense forest | 1.0 |
| Secondary forest | 0.7 |
| Pastures | 0.2 |
| Croplands | 0.1 |
| Water | 0.5 |
| Urban | 0.0 |
| Bare soil | 0.0 |
| **Threats** | |
| Agriculture/pastures: max distance | 5 km |
| Agriculture/pastures: weight | 0.6 |
| Urban: max distance | 10 km |
| Urban: weight | 0.3 |
| Bare soil: max distance | 3 km |
| Bare soil: weight | 0.1 |
| **Model parameters** | |
| Decay function | Exponential |
| Half-saturation constant (k) | 0.5 |
| Scaling parameter (z) | 2.5 |

---

## Table S7. Water yield model parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Precipitation (P) | CHIRPS Daily, annual sum | ~5.5 km |
| Evapotranspiration (ET) | 60% of precipitation (simplified proxy) | Approximation |
| **Baseflow recharge coefficients** | | |
| Dense forest | 0.80 | High infiltration |
| Secondary forest | 0.65 | Moderate infiltration |
| Pastures | 0.40 | Compacted soil |
| Croplands | 0.35 | Moderate |
| Water | 0.00 | Direct surface |
| Urban | 0.10 | Mostly impervious |
| Bare soil | 0.15 | Low organic matter |

Note: The simplified ET approximation was used due to band compatibility issues with MODIS MOD16A2 over the study area. This represents a limitation and should be improved in future studies.

---

## Table S8. GWR model specifications and results

| Parameter | Value |
|-----------|-------|
| Dependent variable | Deforestation rate (% yr^-1), post-agreement |
| Spatial resolution | 1 km grid |
| Sample size (n) | 1,470 |
| Kernel type | Adaptive bisquare |
| Bandwidth | 11 nearest neighbors |
| Bandwidth selection | AICc minimization |
| VIF threshold | 10 (max observed: 5.94 for elevation) |
| **OLS results** | |
| R^2 | 0.143 |
| Adjusted R^2 | 0.138 |
| AIC | -2,249 |
| **GWR results** | |
| Mean R^2 | 0.609 |
| Median R^2 | 0.920 |
| AIC | -7,190 |
| **Improvement** | |
| R^2 gain | +0.466 |
| AIC improvement | 4,942 |
| GWR preferred | Yes |

### VIF diagnostics (all < 10):

| Variable | VIF | OLS coefficient | t-statistic | Significant |
|----------|-----|----------------|-------------|-------------|
| elevation | 5.94 | -0.325 | -11.03 | Yes |
| slope | 1.76 | 0.023 | 1.44 | No |
| dist_rivers | 1.78 | 0.140 | 8.69 | Yes |
| dist_roads | 1.83 | 0.062 | 3.79 | Yes |
| dist_urban | 1.32 | -0.010 | -0.69 | No |
| pop_density | 1.08 | -0.028 | -2.19 | Yes |
| precip | 1.17 | 0.004 | 0.28 | No |
| lst | 5.00 | -0.259 | -9.57 | Yes |

---

## Table S9. CA-Markov scenario specifications

| Parameter | BAU | Conservation | PDET |
|-----------|-----|-------------|------|
| Deforestation rate modifier | 1.0x | 0.5x | 0.7x |
| Forest recovery modifier | 1.0x | 1.3x | 1.15x |
| Agricultural expansion | Current | Reduced | Diversified |
| Projection years | 2030, 2040 | 2030, 2040 | 2030, 2040 |
| Calibration points | 2,938 | 2,938 | 2,938 |
| Active classes | 5 (BDen, BSec, Past, Agua, Urb) | 5 | 5 |

### Transition probability matrix (calibrated from T2-T3-T4):

| From \ To | BDen | BSec | Past | Cult | Agua | Urb | Suel |
|-----------|------|------|------|------|------|-----|------|
| BDen | 0.811 | 0.004 | 0.185 | 0.000 | 0.000 | 0.000 | 0.000 |
| BSec | 0.624 | 0.356 | 0.007 | 0.000 | 0.012 | 0.000 | 0.000 |
| Past | 0.226 | 0.222 | 0.516 | 0.000 | 0.005 | 0.032 | 0.000 |
| Cult | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Agua | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 |
| Urb | 0.255 | 0.346 | 0.109 | 0.000 | 0.151 | 0.139 | 0.000 |
| Suel | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

### Projected compositions (%):

| Scenario | Year | BDen | BSec | Past | Agua | Urb |
|----------|------|------|------|------|------|-----|
| Current | 2024 | 46.9 | 25.1 | 20.1 | 4.3 | 3.6 |
| BAU | 2030 | 59.2 | 14.8 | 19.6 | 5.2 | 1.1 |
| BAU | 2040 | 60.9 | 8.7 | 23.2 | 6.4 | 0.8 |
| Conservation | 2030 | 63.0 | 15.9 | 15.0 | 5.1 | 1.1 |
| Conservation | 2040 | 72.5 | 7.1 | 14.1 | 5.8 | 0.5 |
| PDET | 2030 | 61.5 | 14.9 | 17.3 | 5.1 | 1.1 |
| PDET | 2040 | 67.1 | 7.4 | 18.8 | 6.1 | 0.7 |

---

## Table S10. Hansen GFC forest loss by period

| Period | Label | Hansen loss (ha) | Annual rate (ha/yr) |
|--------|-------|-----------------|---------------------|
| T1 (2010-2013) | Pre-agreement | 64,850 | ~16,213 |
| T2 (2014-2016) | Transition | 94,201 | ~31,400 |
| T3 (2017-2020) | Early post-agreement | 81,106 | ~20,277 |
| T4 (2021-2024) | Recent post-agreement | 26,765 | ~6,691 |

Hansen treecover2000 mean for study area: 69.7%

---

## Table S11. Climate analysis summary (2012-2024)

| Year | Precip (mm) | LST (C) | SPI mean |
|------|-------------|---------|----------|
| 2012 | 2,846 | 27.93 | -- |
| 2013 | 2,985 | 27.72 | -0.14 |
| 2014 | 2,826 | 28.09 | -- |
| 2015 | 2,705 | 28.65 | -- |
| 2016 | 2,706 | 28.36 | -0.86 |
| 2017 | 2,944 | 27.85 | -- |
| 2018 | 2,945 | 27.97 | -- |
| 2019 | 2,739 | 28.24 | -- |
| 2020 | 2,822 | 28.21 | -0.48 |
| 2021 | 2,973 | 27.80 | -- |
| 2022 | 3,464 | 27.01 | -- |
| 2023 | 2,597 | 27.78 | -- |
| 2024 | 2,726 | 27.69 | -0.89 |

Trend analysis (Mann-Kendall):
- Precipitation: tau = -0.051, p = 0.858, Sen's slope = -2.04 mm/yr (not significant)
- LST: tau = -0.333, p = 0.129, Sen's slope = -0.037 C/yr (not significant)
- Mean drought frequency: 0.263 (max: 0.692)

---

## Figure S1. Random Forest feature importance by period

See `outputs/figures/fig_s1_feature_importance.png`

Top 5 features per period:
- T1 (2013): green (89), NDWI (86), swir1 (84), NDVI (79), swir2 (78)
- T2 (2016): swir1 (83), elevation (81), swir2 (78), NDWI (77), blue (71)
- T3 (2020): swir1 (93), slope (89), swir2 (86), NDBI (85), blue (85)
- T4 (2024): swir1 (91), elevation (91), green (88), blue (87), swir2 (83)

---

## Code availability

All analysis scripts are available at: `magdalena_medio_research/scripts/`

| Script | Description |
|--------|-------------|
| `gee_config.py` | Central configuration (study area, periods, classes, parameters) |
| `scripts/utils.py` | Utility functions (compositing, cloud masking, indices) |
| `scripts/01_preprocessing.py` | Image compositing and auxiliary data |
| `scripts/02_training_samples.py` | Training sample generation (Hansen + JRC + GHSL) |
| `scripts/03_classification.py` | Random Forest classification |
| `scripts/04_accuracy_assessment.py` | Accuracy evaluation |
| `scripts/05_change_detection.py` | Transition matrices |
| `scripts/06_fragmentation.py` | Landscape fragmentation metrics |
| `scripts/07_hotspot_analysis.py` | Spatial autocorrelation (Moran's I) and Gi* hotspots |
| `scripts/08_ecosystem_services.py` | Carbon, water yield, habitat quality |
| `scripts/09_climate_analysis.py` | Climate trends and drought indices |
| `scripts/10_gwr_drivers.py` | OLS and GWR driver analysis |
| `scripts/11_ca_markov.py` | CA-Markov future scenario modeling |
| `scripts/12_visualization.py` | Publication-ready figures |
| `run_phase4_figures.py` | Generate all figures and tables from JSON results |

GEE project: `ee-maestria-tesis`
Python environment: conda `magdalena_medio` (Python 3.11, earthengine-api, numpy, scipy)
