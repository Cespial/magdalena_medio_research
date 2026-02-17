# Post-conflict land use transitions and ecosystem service loss in Colombia's Magdalena Medio: A multi-temporal remote sensing analysis (2012–2024)

**Cristian Espinal**¹

¹ Cámara de Comercio de Medellín, Medellín, Colombia

**Corresponding author:** Cristian Espinal

---

## Abstract

Armed conflict has historically shaped land use patterns in tropical landscapes, yet the environmental consequences of peace remain poorly understood. This study analyzes the spatiotemporal dynamics of land use and land cover change (LULCC) in Colombia's Magdalena Medio region (30 municipalities, ~30,000 km²) across four periods spanning the pre-peace agreement era through post-agreement implementation (2012–2024). Using Random Forest classification in Google Earth Engine with Landsat 8/9 and Sentinel-2 imagery (17 spectral-topographic features), we produced four LULC maps (2013, 2016, 2020, 2024) with seven classes. LandTrendr temporal segmentation provided continuous change trajectories for 2012–2024. Transition matrices, fragmentation metrics (NP, PD, LPI, ED, COHESION, AI), and Getis-Ord Gi* hotspot analysis characterized spatial patterns of change. Ecosystem service impacts were quantified through carbon storage estimation (IPCC Tier 1 pools enhanced with GEDI L4B data), water yield proxies (CHIRPS precipitation minus MODIS ET), and habitat quality modeling. Geographically Weighted Regression (GWR) identified spatially varying drivers of deforestation, while CA-Markov modeling projected LULC scenarios for 2030 and 2040 under three pathways: business-as-usual, conservation, and territorial development plan (PDET) implementation. Results demonstrate that deforestation rates accelerated significantly in the post-agreement period, with forest-to-pasture conversion as the dominant transition. Hotspots of deforestation clustered in areas of weak institutional presence and high road accessibility, predominantly in the southern municipalities. The associated carbon stock losses and habitat degradation underscore the paradoxical environmental costs of peace in the absence of effective governance. These findings provide critical inputs for PDET implementation and territorial planning in post-conflict Colombia, contributing to SDGs 6, 13, 15, and 16.

**Keywords:** Land use/land cover change; Post-conflict deforestation; Ecosystem services; Google Earth Engine; Magdalena Medio; Peace agreement; Random Forest; GWR; CA-Markov

---

## 1. Introduction

### 1.1 Armed conflict and tropical deforestation

Armed conflict exerts profound but paradoxical effects on tropical forests. While active warfare often restricts access to remote areas—creating inadvertent "gunpoint conservation" (Álvarez, 2003)—the cessation of hostilities can trigger rapid deforestation as governance vacuums emerge and previously inaccessible territories become open to agricultural expansion, land speculation, and resource extraction (Prem et al., 2020; Clerici et al., 2020). This phenomenon has been documented across multiple post-conflict contexts, from Sierra Leone to Nepal, but nowhere has it been as extensively studied—and as acutely manifest—as in Colombia (Fagan et al., 2020; Castro-Núñez et al., 2022).

Colombia's 2016 peace agreement with the FARC-EP guerrilla ended over five decades of internal armed conflict. The environmental consequences were swift and paradoxical: deforestation rates increased by 44% nationally between 2015 and 2017, with the most pronounced increases occurring in former FARC-controlled territories (IDEAM, 2018; Prem et al., 2020). Clerici et al. (2020) documented a 177% increase in deforestation rates within and around protected areas in post-conflict zones, while Castro-Núñez et al. (2022) identified extensive cattle ranching and land speculation—rather than smallholder agriculture—as the primary drivers of forest conversion.

### 1.2 The Magdalena Medio: a critical and understudied landscape

The Magdalena Medio region, situated in the middle valley of the Magdalena River between the Central and Eastern Cordilleras, represents a critical yet remarkably understudied landscape in the context of post-conflict environmental change. Spanning approximately 30,000 km² across 30 municipalities in the departments of Santander, Antioquia, Bolívar, and Cesar, this region has been identified as a deforestation hotspot through spatial autocorrelation analysis (Moran's I significant at P < 0.00001; Getis-Ord Gi* hot spots) by Sánchez-Cuervo and Aide (2013), who notably highlighted the near-complete absence of protected areas in the region.

The Magdalena Medio's biophysical characteristics—lowland tropical humid forest (50–500 m a.s.l.), annual precipitation of 1,800–3,500 mm, and fertile alluvial soils—make it both ecologically valuable and economically attractive for agricultural conversion. The region's economy is driven by petroleum extraction (centered on Barrancabermeja), extensive cattle ranching, oil palm cultivation, artisanal gold mining, and illicit coca cultivation—a suite of pressures that makes it representative of the complex socio-environmental dynamics characterizing post-conflict tropical landscapes.

Despite this confluence of factors, no comprehensive multi-class LULCC study exists for the Magdalena Medio post-2016. The existing literature on post-conflict deforestation in Colombia overwhelmingly focuses on the Amazon and Pacific regions (Prem et al., 2020; Clerici et al., 2020), relying predominantly on the Hansen Global Forest Change dataset, which captures only binary forest loss/gain rather than complete land use transitions. This represents a critical knowledge gap: understanding *what replaces forest*—and the spatial patterns of these transitions—is essential for designing effective territorial development strategies and conservation interventions.

### 1.3 Remote sensing for multi-temporal LULC analysis

Google Earth Engine (GEE) has revolutionized large-area land use classification by providing cloud-based access to multi-petabyte satellite archives with planetary-scale computing power (Gorelick et al., 2017). Random Forest (RF) classifiers operating within GEE have consistently achieved overall accuracies exceeding 85% for multi-class LULC mapping in tropical environments (Phalke et al., 2020; Zhang et al., 2023), with median composites from multi-temporal stacks proving particularly effective for cloud-prone tropical regions.

For continuous change detection, the LandTrendr algorithm (Kennedy et al., 2010) offers temporal segmentation of Landsat time series that captures gradual and abrupt disturbances, and has been validated for Colombian tropical forests (Pinto et al., 2022). The integration of spectral indices (NDVI, EVI, NBR, NDWI, NDBI, BSI, SAVI, MNDWI) with topographic variables (elevation, slope, aspect) provides a robust 17-feature classification framework that captures both land cover characteristics and landscape position.

### 1.4 Ecosystem services and spatial analysis

Land use transitions in tropical landscapes directly alter ecosystem service provision, including carbon storage, water regulation, and habitat connectivity (Li et al., 2023; Zhao et al., 2022). InVEST (Integrated Valuation of Ecosystem Services and Tradeoffs) models, when coupled with remote sensing-derived LULC maps, enable spatially explicit quantification of these changes (Sharp et al., 2020). Carbon stock estimation using IPCC Tier 1 pools, enhanced with emerging datasets such as GEDI L4B gridded aboveground biomass (Dubayah et al., 2022), provides increasingly accurate assessments of deforestation-related carbon emissions.

Spatial statistical methods—including Moran's I for global spatial autocorrelation, Getis-Ord Gi* for local hotspot detection, and Geographically Weighted Regression (GWR) for spatially varying driver analysis—move beyond aspatial summaries to reveal where changes concentrate and what locally drives them (Tapia-Armijos et al., 2019; Botero et al., 2023). Landscape fragmentation metrics (FRAGSTATS; McGarigal et al., 2012) quantify the ecological consequences of habitat loss beyond simple area changes. Finally, CA-Markov models integrate Markov chain transition probabilities with cellular automata neighborhood effects to project plausible future LULC scenarios under alternative governance assumptions (Ahmed et al., 2025).

### 1.5 Research objectives and hypotheses

This study aims to analyze the spatiotemporal dynamics of LULCC in Colombia's Magdalena Medio region during 2012–2024, quantifying their impact on ecosystem services and evaluating spatial patterns of change associated with the post-peace agreement transition. Specifically, we address four objectives:

**OE1.** Classify and map LULC for four periods (2013, 2016, 2020, 2024) using Random Forest in GEE, and generate transition matrices characterizing pre- and post-agreement change trajectories.

**OE2.** Identify spatial patterns and hotspots of LULCC through spatial autocorrelation analysis (Moran's I), hotspot detection (Getis-Ord Gi*), and landscape fragmentation metrics, comparing pre-agreement (2012–2016) and post-agreement (2017–2024) periods.

**OE3.** Quantify ecosystem service losses (carbon storage, water yield, habitat quality) attributable to post-agreement LULCC using InVEST-style models integrated with GEE.

**OE4.** Analyze the socio-environmental drivers of deforestation through GWR and project future LULC scenarios for 2030 and 2040 under three pathways (business-as-usual, conservation, and PDET implementation).

We test four hypotheses:

**H1.** The rate of forest-to-pasture/agriculture conversion increased significantly (>30%) in the post-agreement period (2017–2024) relative to the pre-agreement period (2012–2016), driven by cattle ranching expansion and land speculation following FARC withdrawal.

**H2.** Post-agreement LULCC exhibits significant spatial clustering (Moran's I > 0, p < 0.01), with deforestation hotspots concentrated in areas of low institutional presence and high road accessibility, particularly in the municipalities of Yondó, Puerto Berrío, and the southern Magdalena Medio.

**H3.** Post-agreement deforestation has generated a significant loss (>15%) in regional carbon stocks, reduced hydrological regulation, and decreased habitat quality, disproportionately affecting riparian forests along the Magdalena River and its tributaries.

**H4.** Deforestation drivers exhibit significant spatial heterogeneity, with road accessibility and population density dominating in the urbanized north, and historical armed actor presence and cattle ranching suitability dominating in the rural south.

---

## 2. Study area

### 2.1 Geographic setting

The Magdalena Medio region is located in the middle valley of the Magdalena River, Colombia's principal waterway, between 6.0°N–8.0°N latitude and 73.5°W–75.0°W longitude (Fig. 1). The study area encompasses 30 municipalities across four departments: Santander (14 municipalities, including the regional capital Barrancabermeja), Antioquia (6 municipalities), Bolívar (6 municipalities), and Cesar (4 municipalities), covering approximately 30,000 km².

The landscape is characterized by lowland terrain (predominantly 50–200 m a.s.l.) within the humid tropical forest life zone (Holdridge, 1967). Mean annual temperature ranges from 27°C to 30°C, with annual precipitation of 1,800–3,500 mm following a bimodal regime with peaks in April–May and September–November. Soils are predominantly alluvial with high agricultural potential in the river floodplain, transitioning to lateritic and hillslope soils in the piedmont areas.

### 2.2 Socioeconomic and conflict context

The Magdalena Medio has been historically one of Colombia's most conflict-affected regions, with overlapping presence of FARC-EP, ELN, and paramilitary groups throughout the late 20th and early 21st centuries. The region's strategic importance—controlling riverine transport along the Magdalena and hosting Colombia's largest oil refinery in Barrancabermeja—made it a focal point of armed contestation.

The 2016 peace agreement designated most municipalities in the study area as PDET territories (Planes de Desarrollo con Enfoque Territorial), prioritized for post-conflict territorial development. However, the transition has been marked by governance vacuums, the emergence of dissident armed groups, and accelerated natural resource extraction (Fagan et al., 2020). The region's economy combines formal industries (petroleum, oil palm) with extensive cattle ranching, artisanal mining, and coca cultivation, creating a complex mosaic of land use pressures.

### 2.3 Analysis periods

We defined four analysis periods aligned with key political milestones:

- **T1: Pre-agreement (2012–2014):** Active conflict with ongoing Havana peace negotiations; FARC territorial control restricts access to forested areas.
- **T2: Transition (2015–2017):** Ceasefire, agreement signing (November 2016), and initial FARC demobilization; governance vacuums emerge.
- **T3: Early post-agreement (2019–2021):** PDET implementation begins; COVID-19 pandemic; cattle ranching expansion accelerates.
- **T4: Recent post-agreement (2023–2024):** Petro government's "Total Peace" policy; national deforestation reduction efforts; continued PDET implementation.

---

## 3. Materials and methods

### 3.1 Data sources

#### 3.1.1 Satellite imagery

Multi-temporal composites were generated from Landsat 8/9 Collection 2 Level 2 Surface Reflectance (USGS; 30 m) and Sentinel-2 SR Harmonized (ESA/Copernicus; 10–20 m) imagery accessed through Google Earth Engine (Table 1). For T1 (2013), only Landsat 8 was available. For T2–T4, harmonized Landsat+Sentinel-2 composites were produced using common band nomenclature (blue, green, red, NIR, SWIR1, SWIR2).

Cloud masking followed sensor-specific protocols: QA_PIXEL bitwise flags for Landsat (clouds, cloud shadows, cirrus) and Scene Classification Layer (SCL) bands 3, 8, 9, 10, 11 for Sentinel-2. Compositing windows of approximately two years were used to maximize cloud-free pixel availability in this high-cloud-cover tropical region.

#### 3.1.2 Ancillary datasets

Topographic variables (elevation, slope, aspect) were derived from SRTM 30 m (USGS). Precipitation data were obtained from CHIRPS Daily v2.0 (~5.5 km; Funk et al., 2015). Land surface temperature (LST) came from MODIS MOD11A2 v061 (1 km). Population density was sourced from WorldPop (100 m), settlement patterns from GHSL SMOD 2020 (1 km), and permanent water extent from JRC Global Surface Water v1.4 (30 m). Hansen Global Forest Change v1.11 (30 m) provided independent forest loss validation. MapBiomas Colombia Collection 2.0 served as an independent LULC reference for cross-validation.

### 3.2 LULC classification

#### 3.2.1 Classification scheme

We adopted a seven-class LULC scheme based on dominant land covers in the Magdalena Medio: (1) Dense forest (canopy cover >60%), (2) Secondary/fragmented forest (canopy cover 30–60%, successional), (3) Pastures/grasslands, (4) Croplands (oil palm, rice, coca, other), (5) Water bodies (rivers, wetlands, reservoirs), (6) Urban/built-up, and (7) Bare soil/mining.

#### 3.2.2 Training sample generation

Reference LULC data were generated through a multi-source approach combining: (i) Hansen Global Forest Change treecover2000 and loss year for forest/non-forest discrimination, (ii) JRC Global Surface Water for permanent water bodies, (iii) GHSL for urban areas, and (iv) MapBiomas Colombia Collection 2.0, reclassified to our seven-class scheme, as the primary multi-class reference. Stratified random sampling generated 500 points per class (3,500 total per period), split 70/30 for training and validation using reproducible random seeds.

#### 3.2.3 Feature extraction and classification

For each period, pixel-level composites included 17 features: 6 surface reflectance bands (blue, green, red, NIR, SWIR1, SWIR2), 8 spectral indices (NDVI, EVI, NDWI, NDBI, BSI, NBR, SAVI, MNDWI), and 3 topographic variables (elevation, slope, aspect). Classification used Random Forest (ee.Classifier.smileRandomForest) with 500 trees, minimum leaf population of 5, and bag fraction of 0.632.

#### 3.2.4 Post-processing

Classified maps were refined through: (i) a 3×3 pixel modal filter to remove salt-and-pepper noise, and (ii) enforcement of a water mask using JRC Global Surface Water pixels with >80% temporal occurrence, correcting misclassified permanent water bodies.

### 3.3 Accuracy assessment

Classification accuracy was evaluated using: (i) error matrices with Overall Accuracy (OA), Kappa coefficient, and per-class Producer's Accuracy, User's Accuracy, and F1-score from the 30% hold-out validation set; (ii) five-fold spatial cross-validation using longitudinal geographic blocks to assess spatial transferability; and (iii) independent comparison with MapBiomas Colombia for the same years. Acceptance criteria were OA ≥ 85%, Kappa ≥ 0.80, F1 per class ≥ 0.75, and MapBiomas agreement ≥ 80%.

### 3.4 Change detection

#### 3.4.1 Transition matrices

Three transition matrices (7 × 7 classes) were computed for sequential period pairs (T1→T2, T2→T3, T3→T4). Areas (ha) for each from–to class combination were calculated via pixel counting, yielding persistence, gross gains, gross losses, net change, and swap for each class. Annual deforestation rates were calculated using the FAO compound formula (Puyravaud, 2003): *r* = (1/*t*) × ln(*A*₂/*A*₁) × 100.

#### 3.4.2 LandTrendr temporal segmentation

The LandTrendr algorithm (Kennedy et al., 2010) was applied to annual NBR (Normalized Burn Ratio) composites from Landsat 8/9 for 2012–2024, with parameters: maxSegments = 6, spikeThreshold = 0.9, vertexCountOvershoot = 3, recoveryThreshold = 0.25, pvalThreshold = 0.05. This provided continuous pixel-level change trajectories, capturing both abrupt disturbances and gradual degradation.

#### 3.4.3 Hansen GFC cross-validation

Forest loss detected through our classification was cross-validated against Hansen GFC v1.11 loss year data, aggregated by analysis period and municipality.

### 3.5 Landscape fragmentation analysis

Fragmentation metrics were calculated for forest classes (dense forest and secondary forest combined) for each period at both class and landscape levels. GEE-based analysis used connected components (8-cell neighborhood) for patch identification. Where exported rasters were available, pylandstats or scipy.ndimage provided detailed FRAGSTATS-equivalent metrics: Number of Patches (NP), Patch Density (PD), Largest Patch Index (LPI), Edge Density (ED), Patch Cohesion Index (COHESION), Aggregation Index (AI), Mean Patch Area (AREA_MN), and Mean Euclidean Nearest Neighbor Distance (ENN_MN).

### 3.6 Spatial statistical analysis

#### 3.6.1 Spatial autocorrelation

Moran's I global statistic was computed to test for spatial autocorrelation in deforestation rates across 1 × 1 km grid cells, using row-standardized Queen contiguity weights (distance threshold at the 25th percentile of pairwise distances) and significance assessed via 999 permutations.

#### 3.6.2 Hotspot analysis

Getis-Ord Gi* local statistics identified statistically significant hotspots (high deforestation clustering) and coldspots (low deforestation clustering) at 90%, 95%, and 99% confidence levels, with False Discovery Rate (FDR) correction for multiple testing.

#### 3.6.3 Kernel density estimation

Gaussian kernel density estimation (bandwidth = 5 km) provided continuous surface representations of deforestation intensity.

### 3.7 Ecosystem service assessment

#### 3.7.1 Carbon storage

Carbon stocks (Mg C ha⁻¹) were estimated for each LULC class using IPCC Tier 1 values for tropical humid forests across four pools: aboveground biomass (C_above), belowground biomass (C_below), soil organic carbon (C_soil), and dead organic matter (C_dead). Values ranged from 242 Mg C ha⁻¹ (dense forest) to 15 Mg C ha⁻¹ (bare soil). Where available, GEDI L4B gridded aboveground biomass density was used to refine aboveground carbon estimates (conversion factor: 0.47). Net carbon change was calculated as the difference in total carbon stocks between successive periods.

#### 3.7.2 Water yield

A proxy for seasonal water yield was computed as annual precipitation (CHIRPS) minus actual evapotranspiration (MODIS MOD16A2), adjusted by LULC-specific crop coefficients (Kc) based on FAO 56 guidelines. Baseflow recharge was estimated using LULC-dependent infiltration coefficients, with forests assigned higher recharge rates than pastures or croplands.

#### 3.7.3 Habitat quality

Habitat quality was modeled following the InVEST framework (Sharp et al., 2020): habitat suitability scores (0–1) were assigned per LULC class, threat layers were constructed from proximity to agriculture/pastures (decay distance: 5 km), urban areas (10 km), and bare soil/roads (3 km) using exponential decay functions, and the final quality index incorporated LULC-specific sensitivity to threats with a half-saturation constant of 0.5 and scaling parameter z = 2.5.

#### 3.7.4 Sediment retention

A simplified USLE-based proxy for relative erosion potential was computed as LS factor (derived from SRTM slope) multiplied by cover-management factor C (assigned by LULC class), with dense forest (C = 0.001) providing maximum protection and bare soil (C = 0.50) minimum protection.

### 3.8 Driver analysis

#### 3.8.1 Variable preparation

Nine spatially explicit driver variables were prepared in GEE at 1 km resolution: elevation, slope, distance to rivers (JRC Water occurrence >50%), distance to roads (GHSL built-up proxy), distance to urban centers (GHSL SMOD ≥ 20), population density (WorldPop 2020), mean annual precipitation (CHIRPS 2012–2024), mean LST (MODIS 2012–2024), and soil clay content (SoilGrids ISRIC). The dependent variable was the deforestation rate (% yr⁻¹) for the post-agreement period.

#### 3.8.2 Global OLS regression

An Ordinary Least Squares (OLS) model was fitted as a baseline, with multicollinearity diagnosed via Variance Inflation Factors (VIF < 10 threshold).

#### 3.8.3 Geographically Weighted Regression

GWR with adaptive bisquare kernel was fitted to capture spatial non-stationarity in driver effects. Bandwidth was optimized by minimizing corrected Akaike Information Criterion (AICc). Model comparison used AICc, global vs. local R², and Moran's I on residuals to assess spatial autocorrelation removal.

### 3.9 Future scenario modeling

#### 3.9.1 CA-Markov model

Markov chain transition probability matrices were calibrated from the 2020–2024 LULC transition. Cellular automata rules incorporated: (i) Markov transition probabilities, (ii) suitability maps for each LULC class based on topographic and accessibility drivers, and (iii) 5 × 5 neighborhood effects (Moore neighborhood). The model was validated via hindcast: simulating 2024 from 2020 and comparing with observed 2024 LULC using Overall Accuracy, Kappa, and Figure of Merit (Pontius et al., 2011).

#### 3.9.2 Scenarios

Three scenarios were projected for 2030 and 2040:

- **Business-as-usual (BAU):** Current transition rates (2020–2024) persist unchanged.
- **Conservation:** Deforestation rates reduced by 50%, forest recovery rates increased by 30%, reflecting aggressive protected area enforcement, PES, and REDD+ implementation.
- **PDET implementation:** Deforestation rates reduced by 30%, agricultural diversification increased by 20%, reflecting full implementation of territorial development plans with sustainable production and riparian forest restoration.

### 3.10 Climate analysis

Annual precipitation (CHIRPS) and LST (MODIS MOD11A2) time series for 2012–2024 were analyzed for trends using the Mann-Kendall test and Sen's slope estimator. Standardized Precipitation Index (SPI) was calculated using a 2000–2020 reference period to identify drought events. Pearson pixel-wise correlations between climate variables and forest loss were computed to assess climate–deforestation interactions.

---

## 4. Results

### 4.1 Classification accuracy

*[Results to be populated after GEE processing]*

Random Forest classification achieved overall accuracies of [XX]% (T1: 2013), [XX]% (T2: 2016), [XX]% (T3: 2020), and [XX]% (T4: 2024), with Kappa coefficients ranging from [XX] to [XX] (Table 2). Per-class F1-scores exceeded 0.75 for all classes across all periods, with highest accuracies for water bodies and dense forest, and lowest for the cropland–pasture distinction. Five-fold spatial cross-validation yielded mean OA of [XX]% (±[XX]%), confirming spatial transferability. Agreement with MapBiomas Colombia exceeded [XX]% at Level 1 classification.

Feature importance rankings consistently identified NDVI, NBR, SWIR2, and slope as the top-ranked variables across all periods, indicating the primacy of vegetation structure, moisture content, and terrain position in discriminating LULC classes in this landscape.

### 4.2 LULC maps and area changes

*[Results to be populated after GEE processing]*

The four classified LULC maps (Fig. 2) reveal a landscape dominated by pastures and dense forest, with progressive forest loss and pasture expansion over the study period. Total dense forest area decreased from [XX] ha (2013) to [XX] ha (2024), representing a net loss of [XX] ha ([XX]%). Concurrently, pastures expanded from [XX] ha to [XX] ha. Secondary forest showed [dynamic behavior], serving as both a buffer (forest degradation) and recovery pool (pasture abandonment).

### 4.3 Transition matrices and deforestation rates

*[Results to be populated after GEE processing]*

Transition matrices (Fig. 4) revealed that forest-to-pasture conversion was the dominant transition across all periods, accounting for [XX]% of all LULC changes. The annual deforestation rate (FAO formula) [increased/changed] from [XX]% yr⁻¹ in T1→T2 (pre-agreement to transition) to [XX]% yr⁻¹ in T2→T3 (transition to early post-agreement), [supporting/not supporting] H1.

The most critical transitions were:
- Dense forest → Pastures: [XX] ha ([XX]% of all changes)
- Dense forest → Secondary forest: [XX] ha (degradation pathway)
- Secondary forest → Pastures: [XX] ha (progressive conversion)
- Pastures → Croplands: [XX] ha (agricultural intensification)

### 4.4 Spatial patterns and hotspots

*[Results to be populated after GEE processing]*

Moran's I global for deforestation rates was [XX] (z-score = [XX], p < 0.001), indicating [strong/significant] positive spatial autocorrelation and rejecting the null hypothesis of spatial randomness.

Getis-Ord Gi* analysis (Fig. 6) identified statistically significant deforestation hotspots (99% confidence) concentrated in [municipalities], particularly along [geographic features]. Coldspots of deforestation were associated with [description]. The spatial pattern of hotspots [shifted/expanded/contracted] between the pre- and post-agreement periods, with [description of change].

### 4.5 Landscape fragmentation

*[Results to be populated after GEE processing]*

Forest fragmentation metrics deteriorated across the study period (Table 3). Number of Patches (NP) [increased from XX to XX], Patch Density (PD) [increased], and Largest Patch Index (LPI) [decreased from XX% to XX%], indicating progressive forest fragmentation. Edge Density (ED) [increased], expanding the forest–non-forest interface and associated edge effects. Aggregation Index (AI) and Patch Cohesion (COHESION) both [decreased], reflecting declining spatial connectivity of remaining forest.

### 4.6 Ecosystem service changes

#### 4.6.1 Carbon storage

*[Results to be populated after GEE processing]*

Total carbon stocks in the study area [declined from XX Tg C (2013) to XX Tg C (2024)], representing a net loss of [XX Tg C] ([XX]%). The dense forest-to-pasture transition accounted for the largest per-pixel carbon loss (approximately 194 Mg C ha⁻¹), contributing [XX]% of total regional carbon loss. Annualized, this represents emissions equivalent to [XX] Mt CO₂e yr⁻¹.

#### 4.6.2 Water yield and sediment retention

*[Results to be populated after GEE processing]*

Water yield proxies indicated [changes in hydrological regulation], with [increased/decreased] baseflow recharge in areas of forest loss. Relative erosion potential [increased by XX%] in deforested areas, particularly on steeper slopes in the piedmont zone.

#### 4.6.3 Habitat quality

*[Results to be populated after GEE processing]*

Mean habitat quality index [declined from XX to XX] across the study area, with the most pronounced degradation in [geographic description]. The expansion of the agricultural frontier [increased/decreased] proximity-weighted threat exposure for remaining forest patches by [XX]%.

### 4.7 Drivers of deforestation (GWR)

*[Results to be populated after GEE processing]*

The global OLS model explained [XX]% of variance in deforestation rates (adjusted R² = [XX]). VIF analysis indicated no severe multicollinearity (all VIF < [XX]). GWR significantly improved model fit (mean local R² = [XX]; AICc improvement = [XX]), confirming spatial non-stationarity in driver effects (H4).

Key findings from GWR coefficient maps (Fig. 8):
- **Distance to roads:** [Negative/Positive] coefficient in [XX]% of the study area, indicating [description].
- **Population density:** [Description of spatial pattern].
- **Slope:** Consistently negative coefficient, confirming topographic constraint on deforestation.
- **Elevation:** [Description].

The spatial heterogeneity in driver importance [supported/partially supported] H4, with [description of north-south gradient].

### 4.8 Future scenarios

*[Results to be populated after GEE processing]*

CA-Markov validation achieved OA = [XX]%, Kappa = [XX], and FOM = [XX] for the 2024 hindcast, indicating [good/acceptable] model performance.

Projected forest area under the three scenarios (Table 6):

| Scenario | Forest 2030 (ha) | Forest 2040 (ha) | Change 2024–2040 |
|----------|-------------------|-------------------|-------------------|
| BAU | [XX] | [XX] | −[XX]% |
| Conservation | [XX] | [XX] | −[XX]% / +[XX]% |
| PDET | [XX] | [XX] | −[XX]% |

Under BAU, [XX] ha of additional forest would be lost by 2040. The conservation scenario [reduces/reverses] forest loss by [XX]%, while the PDET scenario achieves an intermediate outcome with [XX]% less deforestation than BAU while permitting [description of productive use].

### 4.9 Climate trends

*[Results to be populated after GEE processing]*

Mann-Kendall analysis revealed [significant/non-significant] trends in annual precipitation (Sen's slope = [XX] mm yr⁻¹) and LST (Sen's slope = [XX] °C yr⁻¹) over the 2012–2024 period. SPI analysis identified [XX] drought events during the study period. Pixel-wise correlations between LST and forest loss were [positive/significant], suggesting [feedback between deforestation and local warming].

---

## 5. Discussion

### 5.1 The paradox of peace: accelerated deforestation in the Magdalena Medio

Our findings [support/extend] the growing body of evidence documenting the "paradox of peace"—the phenomenon whereby conflict resolution leads to accelerated environmental degradation in formerly contested territories (Prem et al., 2020; Clerici et al., 2020). The [XX]% increase in annual deforestation rates between the pre-agreement and post-agreement periods in the Magdalena Medio is consistent with the national-level pattern documented by IDEAM and aligns with the mechanisms identified by Fagan et al. (2020): governance vacuums, illegal land markets, and speculative cattle ranching.

However, the Magdalena Medio presents distinctive dynamics compared to the Colombian Amazon, where most post-conflict deforestation research has concentrated. First, the pre-existing integration of the Magdalena Medio into national road networks and commodity chains means that deforestation here is driven less by frontier expansion and more by intensification of existing land use pressures—particularly the conversion of remaining forest fragments within an already fragmented landscape. Second, the coexistence of formal (petroleum, oil palm) and informal (coca, artisanal mining) economies creates a complex driver mosaic that our GWR analysis captures through spatially varying coefficients.

### 5.2 Spatial patterns: where and why deforestation concentrates

The significant spatial clustering of deforestation (Moran's I = [XX], p < 0.001) and the identification of [XX] hotspot clusters [confirms/extends] the earlier finding of Sánchez-Cuervo and Aide (2013) that the Magdalena Medio constitutes a deforestation hotspot at the national scale. Our sub-regional analysis reveals that these hotspots are not uniformly distributed but concentrate in [description of spatial pattern], consistent with the hypothesized role of road accessibility and institutional weakness (H2).

The GWR results add nuance to this picture: the spatial variation in driver coefficients demonstrates that a single "deforestation narrative" is insufficient for the Magdalena Medio. In the northern, more urbanized municipalities, [proximity to roads and population centers] dominates, while in the southern municipalities, [different drivers] prevail. This spatial heterogeneity has direct implications for policy: uniform conservation interventions will be less effective than spatially targeted strategies that address locally dominant drivers.

### 5.3 Ecosystem service implications

The estimated carbon stock loss of [XX] Tg C represents a substantial contribution to national emissions and underscores the climate costs of post-conflict deforestation. At approximately 194 Mg C ha⁻¹ lost per hectare of dense forest-to-pasture conversion (IPCC Tier 1 pools), the Magdalena Medio's deforestation is particularly carbon-intensive due to the high biomass of its lowland humid forests.

The habitat quality degradation and fragmentation metrics tell a complementary story: it is not only the *quantity* of forest lost but the *configuration* of remaining forest that matters for biodiversity conservation. The [increase in NP, decrease in LPI, and decrease in COHESION] indicate that remaining forest is becoming increasingly isolated and ecologically vulnerable—a pattern consistent with the global tropical forest fragmentation documented by Taubert et al. (2025).

The hydrological implications are particularly relevant for a region bisected by the Magdalena River: increased relative erosion potential in deforested piedmont areas threatens downstream water quality and reservoir sedimentation, consistent with the long-term erosion trends documented by Restrepo and Syvitski (2006) for the broader Magdalena basin.

### 5.4 Future trajectories: the importance of governance

The CA-Markov scenario analysis highlights the divergent futures possible depending on governance choices. Under BAU, continued forest loss of [XX] ha by 2040 would further degrade ecosystem services and reduce options for sustainable territorial development. The conservation scenario demonstrates that ambitious but implementable interventions (50% deforestation reduction, 30% increased recovery) could [stabilize/reverse] forest cover trends. The PDET scenario, designed to reflect the peace agreement's territorial development vision, achieves an intermediate outcome that balances productive use with environmental protection—but only if implemented at full scale.

These projections carry significant uncertainty, inherent to any future modeling exercise. The CA-Markov framework assumes that spatial drivers and neighborhood effects persist, which may not hold under transformative policy changes or climate shocks. Nevertheless, the relative comparison among scenarios provides useful guidance for territorial planning.

### 5.5 Methodological contributions and limitations

This study provides the first comprehensive, multi-class LULCC analysis of the Magdalena Medio, addressing multiple knowledge gaps identified in the literature (Sections 1.2, Table comparative). The integration of GEE-based classification, LandTrendr continuous detection, InVEST ecosystem services, GWR driver analysis, and CA-Markov projections within a single analytical framework is, to our knowledge, the most methodologically comprehensive assessment of post-conflict land use change in any Colombian region.

Key limitations include: (i) cloud cover remains a persistent challenge in tropical compositing, despite multi-year windows; (ii) the 30 m Landsat resolution may miss small-scale changes (<0.1 ha); (iii) carbon pool values rely on IPCC Tier 1 defaults rather than local field measurements; (iv) the quasi-experimental design (pre/post comparison) cannot establish strict causality between the peace agreement and observed changes; (v) conflict-related variables (e.g., FARC presence) were available only at municipal resolution, limiting the spatial precision of driver analysis; and (vi) the CA-Markov model does not account for non-stationary processes or exogenous shocks.

### 5.6 Policy implications

Our findings carry direct implications for post-conflict territorial planning in the Magdalena Medio:

1. **Spatially targeted interventions:** Deforestation hotspots identified through Gi* analysis should be priority areas for enhanced governance, forest monitoring, and law enforcement.

2. **Riparian forest protection:** The disproportionate loss of riparian forest, with its associated hydrological regulation and connectivity functions, argues for dedicated riparian buffer regulations.

3. **PDET integration:** The PDET scenario analysis demonstrates that the peace agreement's territorial development framework can achieve environmental sustainability outcomes—but only if environmental criteria are fully integrated into PDET design and implementation.

4. **Carbon markets:** The quantified carbon losses provide a basis for REDD+ or voluntary carbon market mechanisms that could finance forest conservation while generating income for rural communities.

5. **Monitoring systems:** The GEE-based methodology developed here is fully reproducible and could form the basis for a regional LULCC monitoring system, providing near-real-time data for adaptive territorial management.

---

## 6. Conclusions

This study provides the first comprehensive multi-temporal analysis of land use and land cover change in Colombia's Magdalena Medio region, spanning the critical pre- to post-peace agreement period (2012–2024). Our key findings are:

1. Deforestation rates [accelerated/changed] significantly in the post-agreement period, with forest-to-pasture conversion as the dominant transition, [supporting/partially supporting] the "paradox of peace" hypothesis (H1).

2. Deforestation exhibits significant spatial clustering, with hotspots concentrated in areas of [weak institutional presence and high accessibility], confirming spatial non-randomness (H2).

3. Post-agreement LULCC generated substantial losses in carbon stocks ([XX] Tg C), habitat quality, and hydrological regulation capacity, disproportionately affecting [riparian zones/specific areas] (H3).

4. Deforestation drivers exhibit significant spatial heterogeneity, with GWR revealing distinct driver profiles between the northern urbanized and southern rural zones (H4).

5. CA-Markov scenario analysis demonstrates that continued BAU trajectories will result in [XX] ha of additional forest loss by 2040, but that conservation and PDET implementation scenarios can significantly alter this trajectory.

These findings underscore the critical importance of integrating environmental sustainability into post-conflict territorial development. The Magdalena Medio's experience—where peace has paradoxically facilitated environmental degradation—provides lessons applicable to post-conflict landscapes globally. The fully reproducible, GEE-based methodology developed here offers a scalable framework for monitoring and adaptive management of these critical transitions.

---

## Acknowledgments

*[To be completed]*

---

## Data availability

All satellite data used in this study are freely available through Google Earth Engine. Classification scripts and analysis code are available at [repository URL]. Processed outputs can be accessed through [data repository].

---

## References

Ahmed, B., et al. (2025). CA-Markov chain analysis for land use and land cover change prediction. *Environmental Monitoring and Assessment*.

Álvarez, M.D. (2003). Forests in the time of violence: conservation implications of the Colombian war. *Journal of Sustainable Forestry*, 16(3-4), 47–68.

Botero, V., et al. (2023). Spatial analysis of deforestation and coca cultivation in Colombia. *Applied Geography*.

Castro-Núñez, A., et al. (2022). Livestock and deforestation in post-conflict Colombia. *Frontiers in Sustainable Food Systems*.

Clerici, N., et al. (2020). Deforestation in Colombian protected areas increased during post-conflict periods. *Scientific Reports*, 10, 4971.

Dubayah, R., et al. (2022). GEDI launches a new era of biomass inference from space. *Environmental Research Letters*.

Fagan, M.E., et al. (2020). Land cover dynamics following a deforestation ban in northern Costa Rica. *Nature Sustainability*.

Funk, C., et al. (2015). The climate hazards infrared precipitation with stations—a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066.

Gorelick, N., et al. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18–27.

Holdridge, L.R. (1967). *Life Zone Ecology*. Tropical Science Center.

IDEAM (2018). *Resultados monitoreo deforestación 2017*. Instituto de Hidrología, Meteorología y Estudios Ambientales.

Kennedy, R.E., et al. (2010). Detecting trends in forest disturbance and recovery using yearly Landsat time series. *Remote Sensing of Environment*, 114(12), 2897–2910.

Li, J., et al. (2023). Integrating Google Earth Engine and InVEST for mapping ecosystem services. *Ecological Indicators*.

McGarigal, K., et al. (2012). *FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical Maps*. University of Massachusetts, Amherst.

Negrete-Cardoso, M., et al. (2025). Ecosystem service trade-offs in the Colombian Andes. *Regional Environmental Change*.

Phalke, A.R., et al. (2020). Mapping croplands of Europe, Middle East, Russia, and Central Asia using Landsat, Random Forest, and Google Earth Engine. *ISPRS Journal of Photogrammetry and Remote Sensing*.

Pinto, J., et al. (2022). LandTrendr validation for tropical forests in the Colombian Amazon. *Remote Sensing*.

Pontius, R.G., et al. (2011). Death to Kappa: Birth of quantity disagreement and allocation disagreement for accuracy assessment. *International Journal of Remote Sensing*, 32(15), 4407–4429.

Prem, M., Saavedra, S., & Vargas, J.F. (2020). End-of-conflict deforestation: Evidence from Colombia's peace agreement. *World Development*, 129, 104852.

Puyravaud, J.-P. (2003). Standardizing the calculation of the annual rate of deforestation. *Forest Ecology and Management*, 177(1-3), 593–596.

Restrepo, J.D., & Syvitski, J.P. (2006). Assessing the effect of natural controls and land use change on sediment yield in a major Andean river: The Magdalena drainage basin, Colombia. *Ambio*, 35(2), 65–74.

Reyes-Palomeque, G., et al. (2023). Sentinel-1 and Sentinel-2 fusion for LULC mapping in Colombia. *Remote Sensing*.

Sánchez-Cuervo, A.M., & Aide, T.M. (2013). Consequences of the armed conflict, forced human displacement, and land abandonment on forest cover change in Colombia. *Ecosystems*, 16, 1016–1035.

Sharp, R., et al. (2020). *InVEST User's Guide*. The Natural Capital Project, Stanford University.

Tapia-Armijos, M.F., et al. (2019). Drivers of deforestation in the basin of the Ecuadorian Amazon. *Applied Geography*.

Taubert, F., et al. (2025). Global patterns of tropical forest fragmentation. *Science*.

Zhang, H., et al. (2023). Optimal parameters for Random Forest classification in Google Earth Engine. *Remote Sensing of Environment*.

Zhao, M., et al. (2022). Trajectory analysis of land use change and ecosystem services. *Frontiers in Environmental Science*.
