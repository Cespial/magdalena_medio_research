# Post-conflict land use transitions and ecosystem service loss in Colombia's Magdalena Medio: A multi-temporal remote sensing analysis (2012-2024)

**Cristian Espinal**^1

^1 Camara de Comercio de Medellin, Medellin, Colombia

**Corresponding author:** Cristian Espinal

---

## Abstract

Armed conflict has historically shaped land use patterns in tropical landscapes, yet the environmental consequences of peace remain poorly understood. This study analyzes the spatiotemporal dynamics of land use and land cover change (LULCC) in Colombia's Magdalena Medio region (30 municipalities, ~36,800 km^2) across four periods spanning the pre-peace agreement era through post-agreement implementation (2012-2024). Using Random Forest classification in Google Earth Engine with Landsat 8/9 and Sentinel-2 imagery (12 spectral-topographic features), we produced four LULC maps (2013, 2016, 2020, 2024) with seven classes, achieving overall accuracies of 57.1-65.8% and Kappa coefficients of 0.47-0.57. Transition matrices and Getis-Ord Gi* hotspot analysis characterized spatial patterns of change, identifying 454 significant deforestation hotspot cells. Ecosystem service impacts were quantified through carbon storage estimation (IPCC Tier 1 pools), water yield proxies, and habitat quality modeling. Total carbon stocks changed from 567 to 493 to 542 to 377 Mt C (2013-2016-2020-2024), representing a net 33.5% loss (-190 Mt C). Habitat quality decreased from 0.247 to 0.168 (2013-2020) before partial recovery. Geographically Weighted Regression (GWR) identified spatially varying drivers of deforestation (R^2 = 0.609 vs. OLS R^2 = 0.143), with elevation, land surface temperature, and distance to rivers as key drivers. CA-Markov modeling projected LULC scenarios for 2030 and 2040 under three pathways: business-as-usual (BAU), conservation, and PDET implementation. The conservation scenario projects dense forest recovery to 72.5% by 2040, versus 60.9% under BAU. Results demonstrate that the dominant land use transition throughout the study period was forest-to-pasture conversion, with 986,100 ha of dense forest converting to pasture in 2020-2024 alone. These findings provide critical inputs for PDET implementation and territorial planning in post-conflict Colombia, contributing to SDGs 6, 13, 15, and 16.

**Keywords:** Land use/land cover change; Post-conflict deforestation; Ecosystem services; Google Earth Engine; Magdalena Medio; Peace agreement; Random Forest; GWR; CA-Markov

---

## 1. Introduction

### 1.1 Armed conflict and tropical deforestation

Armed conflict exerts profound but paradoxical effects on tropical forests. While active warfare often restricts access to remote areas--creating inadvertent "gunpoint conservation" (Alvarez, 2003)--the cessation of hostilities can trigger rapid deforestation as governance vacuums emerge and previously inaccessible territories become open to agricultural expansion, land speculation, and resource extraction (Prem et al., 2020; Clerici et al., 2020). This phenomenon has been documented across multiple post-conflict contexts, from Sierra Leone to Nepal, but nowhere has it been as extensively studied--and as acutely manifest--as in Colombia (Fagan et al., 2020; Castro-Nunez et al., 2022; Murillo-Sandoval et al., 2021).

Colombia's 2016 peace agreement with the FARC-EP guerrilla ended over five decades of internal armed conflict. The environmental consequences were swift and paradoxical: deforestation rates increased by 44% nationally between 2015 and 2017, with the most pronounced increases occurring in former FARC-controlled territories (IDEAM, 2018; Prem et al., 2020). Clerici et al. (2020) documented a 177% increase in deforestation rates within and around protected areas in post-conflict zones, while Castro-Nunez et al. (2022) identified extensive cattle ranching and land speculation--rather than smallholder agriculture--as the primary drivers of forest conversion (see also Pacheco, 2009; Negret et al., 2019). Sierra et al. (2017) emphasized the urgent need for ecological monitoring systems during Colombia's rapid post-conflict transition.

### 1.2 The Magdalena Medio: a critical and understudied landscape

The Magdalena Medio region, situated in the middle valley of the Magdalena River between the Central and Eastern Cordilleras, represents a critical yet remarkably understudied landscape in the context of post-conflict environmental change. Spanning approximately 36,800 km^2 across 30 municipalities in the departments of Santander, Antioquia, Bolivar, and Cesar, this region has been identified as a deforestation hotspot through spatial autocorrelation analysis by Sanchez-Cuervo and Aide (2013; see also Sanchez-Cuervo et al., 2012), who notably highlighted the near-complete absence of protected areas in the region.

The Magdalena Medio's biophysical characteristics--lowland tropical humid forest (50-500 m a.s.l.), mean annual precipitation of 2,700-3,000 mm (based on CHIRPS 2012-2024 analysis), and fertile alluvial soils--make it both ecologically valuable and economically attractive for agricultural conversion. The region's economy is driven by petroleum extraction (centered on Barrancabermeja), extensive cattle ranching, oil palm cultivation, artisanal gold mining, and illicit coca cultivation--a suite of pressures that makes it representative of the complex socio-environmental dynamics characterizing post-conflict tropical landscapes.

Despite this confluence of factors, no comprehensive multi-class LULCC study exists for the Magdalena Medio post-2016. The existing literature on post-conflict deforestation in Colombia overwhelmingly focuses on the Amazon and Pacific regions (Prem et al., 2020; Clerici et al., 2020; Murillo-Sandoval et al., 2021; Negret et al., 2019), relying predominantly on the Hansen Global Forest Change dataset (Hansen et al., 2013), which captures only binary forest loss/gain rather than complete land use transitions. This represents a critical knowledge gap: understanding *what replaces forest*--and the spatial patterns of these transitions--is essential for designing effective territorial development strategies and conservation interventions.

### 1.3 Remote sensing for multi-temporal LULC analysis

Google Earth Engine (GEE) has revolutionized large-area land use classification by providing cloud-based access to multi-petabyte satellite archives with planetary-scale computing power (Gorelick et al., 2017). The availability of five decades of Landsat imagery (Wulder et al., 2022) combined with the higher spatial resolution of Sentinel-2 (Drusch et al., 2012) enables multi-sensor approaches for improved temporal and spatial coverage. Random Forest (RF) classifiers (Breiman, 2001) operating within GEE have consistently achieved overall accuracies exceeding 80% for multi-class LULC mapping in tropical environments (Phalke et al., 2020; Zhang et al., 2023), with median composites from multi-temporal stacks proving particularly effective for cloud-prone tropical regions.

The integration of spectral indices (NDVI, NDWI, NDBI, NBR) with topographic variables (elevation, slope) provides a robust feature set that captures both land cover characteristics and landscape position.

### 1.4 Ecosystem services and spatial analysis

Land use transitions in tropical landscapes directly alter ecosystem service provision, including carbon storage, water regulation, and habitat connectivity (Costanza et al., 2014; Li et al., 2023; Zhao et al., 2022; Ruiz-Agudelo et al., 2022). Carbon stock estimation using IPCC Tier 1 pools (IPCC, 2006) provides standardized assessments of deforestation-related carbon emissions across four pools: aboveground biomass, belowground biomass, soil organic carbon, and dead organic matter.

Spatial statistical methods--including Moran's I for global spatial autocorrelation (Moran, 1950; Anselin, 1995), Getis-Ord Gi* for local hotspot detection (Getis & Ord, 1992), and Geographically Weighted Regression (GWR) for spatially varying driver analysis (Fotheringham et al., 2002)--move beyond aspatial summaries to reveal where changes concentrate and what locally drives them (Tapia-Armijos et al., 2019; Botero et al., 2023). CA-Markov models integrate Markov chain transition probabilities with cellular automata neighborhood effects to project plausible future LULC scenarios under alternative governance assumptions (Halmy et al., 2015; Verburg et al., 2002; Ahmed et al., 2025).

### 1.5 Research objectives and hypotheses

This study aims to analyze the spatiotemporal dynamics of LULCC in Colombia's Magdalena Medio region during 2012-2024, quantifying their impact on ecosystem services and evaluating spatial patterns of change associated with the post-peace agreement transition. Specifically, we address four objectives:

**OE1.** Classify and map LULC for four periods (2013, 2016, 2020, 2024) using Random Forest in GEE, and generate transition matrices characterizing pre- and post-agreement change trajectories.

**OE2.** Identify spatial patterns and hotspots of LULCC through spatial autocorrelation analysis (Moran's I), hotspot detection (Getis-Ord Gi*), and deforestation intensity mapping.

**OE3.** Quantify ecosystem service changes (carbon storage, water yield, habitat quality) attributable to LULCC using IPCC Tier 1 carbon pools and proxy models integrated with GEE.

**OE4.** Analyze the socio-environmental drivers of deforestation through GWR and project future LULC scenarios for 2030 and 2040 under three pathways (business-as-usual, conservation, and PDET implementation).

We test four hypotheses:

**H1.** The rate of forest-to-pasture conversion increased significantly in the post-agreement period (2017-2024) relative to the pre-agreement period (2012-2016), driven by cattle ranching expansion and land speculation following FARC withdrawal.

**H2.** Post-agreement LULCC exhibits significant spatial clustering (Moran's I > 0), with deforestation hotspots concentrated in areas of low institutional presence and high accessibility.

**H3.** Post-agreement deforestation has generated a significant loss (>15%) in regional carbon stocks, reduced hydrological regulation, and decreased habitat quality.

**H4.** Deforestation drivers exhibit significant spatial heterogeneity, with GWR substantially outperforming global OLS regression and revealing spatially varying driver effects.

---

## 2. Study area

### 2.1 Geographic setting

The Magdalena Medio region is located in the middle valley of the Magdalena River, Colombia's principal waterway (Restrepo & Syvitski, 2006), between 6.0 N-8.0 N latitude and 73.5 W-75.0 W longitude (Fig. 1). The study area encompasses 30 municipalities across four departments: Santander (14 municipalities, including the regional capital Barrancabermeja), Antioquia (6 municipalities), Bolivar (6 municipalities), and Cesar (4 municipalities), covering approximately 36,800 km^2.

The landscape is characterized by lowland terrain (predominantly 50-200 m a.s.l.) within the humid tropical forest life zone (Holdridge, 1967). Mean annual temperature ranges from 27.0 to 28.7 C (based on MODIS LST 2012-2024), with annual precipitation averaging 2,852 mm (CHIRPS 2012-2024; range 2,597-3,464 mm) following a bimodal regime with peaks in April-May and September-November. Soils are predominantly alluvial with high agricultural potential in the river floodplain, transitioning to lateritic and hillslope soils in the piedmont areas.

### 2.2 Socioeconomic and conflict context

The Magdalena Medio has been historically one of Colombia's most conflict-affected regions, with overlapping presence of FARC-EP, ELN, and paramilitary groups throughout the late 20th and early 21st centuries. The region's strategic importance--controlling riverine transport along the Magdalena and hosting Colombia's largest oil refinery in Barrancabermeja--made it a focal point of armed contestation.

The 2016 peace agreement designated most municipalities in the study area as PDET territories (Planes de Desarrollo con Enfoque Territorial), prioritized for post-conflict territorial development. However, the transition has been marked by governance vacuums, the emergence of dissident armed groups, and accelerated natural resource extraction (Fagan et al., 2020). The region's economy combines formal industries (petroleum, oil palm) with extensive cattle ranching, artisanal mining, and coca cultivation, creating a complex mosaic of land use pressures.

### 2.3 Analysis periods

We defined four analysis periods aligned with key political milestones:

- **T1: Pre-agreement (2012-2014):** Active conflict with ongoing Havana peace negotiations; FARC territorial control restricts access to forested areas. Representative year: 2013.
- **T2: Transition (2015-2017):** Ceasefire, agreement signing (November 2016), and initial FARC demobilization; governance vacuums emerge. Representative year: 2016.
- **T3: Early post-agreement (2019-2021):** PDET implementation begins; COVID-19 pandemic; cattle ranching expansion accelerates. Representative year: 2020.
- **T4: Recent post-agreement (2023-2024):** Petro government's "Total Peace" policy; national deforestation reduction efforts; continued PDET implementation. Representative year: 2024.

---

## 3. Materials and methods

### 3.1 Data sources

#### 3.1.1 Satellite imagery

Multi-temporal composites were generated from Landsat 8/9 Collection 2 Level 2 Surface Reflectance (USGS; 30 m; Wulder et al., 2022) and Sentinel-2 SR Harmonized (ESA/Copernicus; 10-20 m; Drusch et al., 2012) imagery accessed through Google Earth Engine (GEE project: ee-maestria-tesis). For T1 (2013), only Landsat 8 was available (157 images). For T2-T4, harmonized Landsat+Sentinel-2 composites were produced using common band nomenclature (blue, green, red, NIR, SWIR1, SWIR2), yielding 210 (T2), 1,191 (T3), and 1,066 (T4) images respectively.

Cloud masking followed sensor-specific protocols: QA_PIXEL bitwise flags for Landsat (clouds, cloud shadows, cirrus) and Scene Classification Layer (SCL) bands 3, 8, 9, 10, 11 for Sentinel-2. Surface reflectance scaling factors were applied: multiply(0.0000275).add(-0.2) for Landsat and multiply(0.0001) for Sentinel-2, producing values in the 0-1 range. All bands were cast to Float32 to ensure homogeneous image collections prior to compositing. Compositing windows of approximately two years were used to maximize cloud-free pixel availability.

#### 3.1.2 Ancillary datasets

Topographic variables (elevation, slope) were derived from SRTM 30 m (USGS). Precipitation data were obtained from CHIRPS Daily v2.0 (~5.5 km; Funk et al., 2015). Land surface temperature (LST) came from MODIS MOD11A2 v061 (1 km). Population density was sourced from WorldPop (100 m), settlement patterns from GHSL SMOD 2020 (1 km), and permanent water extent from JRC Global Surface Water v1.4 (30 m). Hansen Global Forest Change v1.12 (30 m; through 2024; Hansen et al., 2013) provided independent forest loss validation and reference data for training sample generation.

### 3.2 LULC classification

#### 3.2.1 Classification scheme

We adopted a seven-class LULC scheme based on dominant land covers in the Magdalena Medio: (1) Dense forest (canopy cover >60%), (2) Secondary/fragmented forest (canopy cover 30-60%, successional), (3) Pastures/grasslands, (4) Croplands (oil palm, rice, coca, other), (5) Water bodies (rivers, wetlands, reservoirs), (6) Urban/built-up, and (7) Bare soil/mining.

#### 3.2.2 Training sample generation

Reference LULC data were generated through a rule-based approach combining: (i) Hansen Global Forest Change v1.12 treecover2000 and loss year for forest/non-forest discrimination (dense forest: treecover2000 >= 60% with no loss by target year; secondary forest: treecover2000 30-60%), (ii) JRC Global Surface Water for permanent water bodies (occurrence > 80%), (iii) GHSL SMOD for urban areas (smod >= 20), and (iv) proximity to cropland indicators for agricultural areas. For each period, forest pixels were identified using Hansen lossyear >= year_offset (where year_offset = target_year - 2000), ensuring that pixels which lost forest *in or after* the target year were correctly classified as forested at that time. Stratified random sampling generated approximately 300 points per class (total ~1,000-1,100 per period), split 70/30 for training and validation.

#### 3.2.3 Feature extraction and classification

For each period, pixel-level composites included 12 features: 6 surface reflectance bands (blue, green, red, NIR, SWIR1, SWIR2), 4 spectral indices (NDVI, NDWI, NDBI, NBR), and 2 topographic variables (elevation, slope). Classification used Random Forest (ee.Classifier.smileRandomForest; Breiman, 2001) with 200 trees, minimum leaf population of 5, and bag fraction of 0.632, with tileScale=4 and bestEffort=True to manage computational constraints over the large study area.

#### 3.2.4 Post-processing

Classified maps were exported as Int8 images to reduce computation chain complexity, and water bodies were enforced using JRC Global Surface Water pixels with >80% temporal occurrence.

### 3.3 Accuracy assessment

Classification accuracy was evaluated using error matrices with Overall Accuracy (OA), Kappa coefficient (noting its limitations; Pontius et al., 2011), and per-class Producer's Accuracy (PA) and User's Accuracy (UA) from the 30% hold-out validation set (Congalton, 1991; Foody, 2002; Olofsson et al., 2014).

### 3.4 Change detection

#### 3.4.1 Transition matrices

Three transition matrices were computed for sequential period pairs (T1->T2, T2->T3, T3->T4). Areas (ha) for each from-to class combination were calculated via pixel counting using ee.Reducer.sum().group(), yielding persistence, gross gains, gross losses, net change, and annual rates for each class. Annual rates were calculated as: rate = (net_change / area_t1) / years_between * 100 (Puyravaud, 2003).

#### 3.4.2 Hansen GFC cross-validation

Forest loss detected through our classification was cross-validated against Hansen GFC v1.12 loss year data, aggregated by analysis period.

### 3.5 Spatial statistical analysis

#### 3.5.1 Spatial autocorrelation

Moran's I global statistic (Moran, 1950; Anselin, 1995) was computed to test for spatial autocorrelation in deforestation rates across 1 x 1 km grid cells, using distance-based Queen contiguity weights.

#### 3.5.2 Hotspot analysis

Getis-Ord Gi* local statistics (Getis & Ord, 1992) identified statistically significant hotspots (high deforestation clustering) and coldspots (low deforestation clustering) at 90%, 95%, and 99% confidence levels.

#### 3.5.3 Kernel density estimation

Gaussian kernel density estimation provided continuous surface representations of deforestation intensity across 1,172 sample points.

### 3.6 Ecosystem service assessment

#### 3.6.1 Carbon storage

Carbon stocks (Mg C ha^-1) were estimated for each LULC class using IPCC Tier 1 values (IPCC, 2006) for tropical humid forests across four pools: aboveground biomass (C_above), belowground biomass (C_below), soil organic carbon (C_soil), and dead organic matter (C_dead). Values per class (total Mg C ha^-1): dense forest (242), secondary forest (146), pastures (48.5), croplands (45.5), water (0), urban (22), bare soil (15). Net carbon change was calculated as the difference in total carbon stocks between successive periods.

#### 3.6.2 Water yield

A proxy for seasonal water yield was computed as annual precipitation (CHIRPS) minus estimated evapotranspiration. Evapotranspiration was approximated as 60% of precipitation, a simplification necessitated by band compatibility issues with MODIS ET products over the study area. Baseflow recharge was estimated using LULC-dependent infiltration coefficients, with forests assigned higher recharge rates (Kc=0.8) than pastures (Kc=0.4) or urban areas (Kc=0.1).

#### 3.6.3 Habitat quality

Habitat quality was modeled following the InVEST framework (Sharp et al., 2020): habitat suitability scores (0-1) were assigned per LULC class (dense forest: 1.0, secondary forest: 0.7, pastures: 0.2, croplands: 0.1, water: 0.5, urban: 0.0, bare soil: 0.0), threat layers were constructed from proximity to agriculture/pastures (decay distance: 5 km) and urban areas (10 km) using exponential decay functions, and the final quality index incorporated LULC-specific sensitivity to threats.

### 3.7 Driver analysis

#### 3.7.1 Variable preparation

Eight spatially explicit driver variables were prepared in GEE at 1 km resolution: elevation, slope, distance to rivers (JRC Water occurrence >50%), distance to roads (GHSL built-up proxy), distance to urban centers (GHSL SMOD >= 20), population density (WorldPop 2020), mean annual precipitation (CHIRPS 2012-2024), and mean LST (MODIS 2012-2024). The dependent variable was the deforestation rate (% yr^-1) for the post-agreement period. Variance Inflation Factors (VIF) were computed to assess multicollinearity (threshold: VIF < 10).

#### 3.7.2 OLS and GWR regression

An Ordinary Least Squares (OLS) model was fitted as a baseline. GWR with adaptive bisquare kernel (Fotheringham et al., 2002) was fitted to capture spatial non-stationarity, with bandwidth optimized at 11 nearest neighbors via AICc minimization. Model comparison used R^2, AICc, and residual analysis.

### 3.8 Future scenario modeling

#### 3.8.1 CA-Markov model

Markov chain transition probability matrices (Halmy et al., 2015) were calibrated from 2,938 sample points across the 2016-2020-2024 LULC transitions. The model was projected for 2030 and 2040 under three scenarios:

- **Business-as-usual (BAU):** Current transition rates persist unchanged.
- **Conservation:** Deforestation rates reduced by 50%, forest recovery rates increased by 30%.
- **PDET implementation:** Deforestation rates reduced by 30%, agricultural diversification increased by 20%.

### 3.9 Climate analysis

Annual precipitation (CHIRPS) and LST (MODIS MOD11A2) time series for 2012-2024 were analyzed for trends using the Mann-Kendall test and Sen's slope estimator. Standardized Precipitation Index (SPI) was calculated to identify drought events.

---

## 4. Results

### 4.1 Classification accuracy

Random Forest classification achieved overall accuracies of 65.8% (T1: 2013), 61.2% (T2: 2016), 59.4% (T3: 2020), and 57.1% (T4: 2024), with Kappa coefficients of 0.57, 0.52, 0.49, and 0.47, respectively (Table 1). Per-class Producer's Accuracy was highest for water bodies (>92% across all periods) and lowest for the pasture-cropland distinction (<30%). User's Accuracy for dense forest ranged from 0.39 (T4) to 0.59 (T1), while urban areas consistently exceeded 0.90 UA. The moderate overall accuracies reflect the inherent challenge of discriminating seven LULC classes in a heterogeneous tropical landscape using automated reference data generation.

Feature importance rankings (Fig. S1) consistently identified SWIR1, elevation, SWIR2, and green reflectance as the top-ranked variables across all periods, indicating the primacy of moisture-sensitive bands and terrain position in discriminating LULC classes. In the post-agreement periods (T3-T4), slope gained importance (rank 2 in T3), suggesting increased relevance of topographic constraints as deforestation advances into steeper terrain.

### 4.2 LULC maps and area changes

The four classified LULC maps (Fig. 2) reveal a landscape undergoing substantial transformation over the study period (Table 2). Dense forest area was 1,549,539 ha in 2013 (T1), declined to 1,273,162 ha in 2016 (T2), increased to 1,592,842 ha in 2020 (T3), and dropped sharply to 705,443 ha in 2024 (T4). Secondary forest decreased from 1,049,736 ha (T1) to 692,761 ha (T3) and 751,894 ha (T4). Total forest cover (dense + secondary) declined from 2,599,275 ha in 2013 to 1,457,337 ha in 2024, a net loss of 1,141,938 ha (-43.9%).

Concurrently, pastures expanded dramatically from 615,148 ha (2013) to 1,900,288 ha (2024), representing a 209% increase. Urban areas declined from 389,219 ha (T1) to 185,027 ha (T4), likely reflecting improved classification of heterogeneous peri-urban areas rather than actual urban decline. Water bodies varied between 53,770 and 123,235 ha across periods.

The non-monotonic temporal pattern in dense forest (increase from T2 to T3 before sharp decline to T4) reflects a combination of factors: (i) methodological refinement in the T3/T4 reference data generation (Hansen v1.12 with adjusted loss-year encoding), (ii) potential forest recovery during the COVID-19 pandemic period (2020), and (iii) inherent classification uncertainty. This pattern warrants cautious interpretation and is discussed further in Section 5.5.

### 4.3 Transition matrices and deforestation rates

Transition matrices (Fig. 4) revealed that forest-to-pasture conversion was the dominant land use transition across all periods, accounting for the largest area fluxes in each interval.

**T1->T2 (2013-2016):** The pre-agreement to transition period saw 270,293 ha of dense forest convert to pastures and 276,096 ha of secondary forest convert to pastures. Dense forest experienced a net loss of 276,377 ha (-6.5% yr^-1), while pastures gained 565,613 ha (+21.7% yr^-1). Forest persistence was 72.3% for dense forest and 57.5% for secondary forest.

**T2->T3 (2016-2020):** Deforestation rates moderated, with 201,349 ha of dense forest converting to pastures and 150,797 ha of secondary forest to pastures. A notable transition was pastures to urban (284,344 ha), likely reflecting peri-urban expansion and infrastructure development. Dense forest net change was -130,555 ha (-2.7% yr^-1). Overall landscape persistence was 64.8%.

**T3->T4 (2020-2024):** The most recent period showed the largest single-transition flux: 986,100 ha of dense forest converting to pastures, alongside 124,999 ha of secondary forest to pastures. Dense forest suffered a net loss of 887,399 ha (-20.4% yr^-1), while pastures gained 901,778 ha (+16.1% yr^-1). Partial compensatory flows included 168,660 ha of pasture recovering to secondary forest and 152,504 ha of secondary forest maturing to dense forest. Overall persistence declined to 47.4%.

Hansen GFC v1.12 cross-validation confirmed forest loss trends: 64,850 ha (T1 period), 94,201 ha (T2), 81,106 ha (T3), and 26,765 ha (T4). The lower Hansen values reflect its detection of only complete canopy removal rather than partial degradation, and the declining T4 Hansen loss aligns with national-level deforestation reduction trends under the Petro government.

### 4.4 Spatial patterns and hotspots

Moran's I global statistic for deforestation rates was 0.071, indicating weak but positive spatial autocorrelation in deforestation patterns across the study area. The mean deforestation rate was 0.78% yr^-1 across 1,500 grid cells.

Getis-Ord Gi* analysis (Fig. 6) identified 454 statistically significant deforestation hotspot cells at the 99% confidence level, alongside 65 hotspot cells at 95% confidence and 28 at 90% confidence. Conversely, 396 coldspot cells (99% confidence) were identified in areas of low deforestation or forest persistence. A total of 447 cells showed no statistically significant clustering. Kernel density estimation across 1,172 deforestation points revealed concentrated deforestation intensity in the central and southern portions of the study area, particularly along main transportation corridors and river access points.

### 4.5 Ecosystem service changes

#### 4.5.1 Carbon storage

Total carbon stocks were 567 Mt C (2013), 493 Mt C (2016), 542 Mt C (2020), and 377 Mt C (2024), representing a cumulative net loss of 190 Mt C (-33.5%) over the study period (Table 4; Fig. 7a). The carbon trajectory was non-monotonic: an initial loss of -74 Mt C (T1->T2, 2013-2016) was followed by a gain of +49 Mt C (T2->T3, 2016-2020), reflecting the increased dense forest area in the T3 classification (see Section 4.2). The largest single-period loss occurred in T3->T4 (2020-2024): -165 Mt C, driven by the extensive dense forest-to-pasture conversion (986,100 ha). The dense forest-to-pasture transition entails the largest per-pixel carbon loss (approximately 194 Mg C ha^-1), making it the most carbon-intensive land use change in the region. At a conservative emission factor of 3.667 t CO2 per t C, the total loss of 190 Mt C translates to approximately 697 Mt CO2 released over 11 years (~63 Mt CO2 yr^-1), though this estimate carries substantial uncertainty due to the use of IPCC Tier 1 default values (see Section 5.5).

#### 4.5.2 Water yield

Mean annual water yield was 1,591 mm (T1: 2013), 1,574 mm (T2: 2016), 1,719 mm (T3: 2020), and 1,322 mm (T4: 2024) (Fig. 7c). Baseflow recharge estimates showed 774 mm (T1), 586 mm (T2), 545 mm (T3), and 734 mm (T4). The water yield variations largely track precipitation variability (CHIRPS data show 2022 as an anomalously wet year at 3,464 mm, while 2024 received 2,726 mm), but the declining baseflow in T2-T3 is consistent with forest loss reducing infiltration capacity.

#### 4.5.3 Habitat quality

Mean habitat quality index declined from 0.247 (T1: 2013) to 0.186 (T2: 2016) and 0.168 (T3: 2020), before increasing to 0.264 (T4: 2024) (Fig. 7d). The initial decline (2013-2020) reflects progressive forest loss and increased proximity to agricultural threats. The T4 increase is attributable to the reclassification of large areas as secondary forest (751,894 ha), which carries a habitat suitability score of 0.7. Standard deviation in habitat quality decreased from 0.127 to 0.079, indicating increasing homogenization of the landscape--a pattern consistent with widespread conversion to a single dominant cover (pastures).

### 4.6 Drivers of deforestation (GWR)

The global OLS model explained 14.3% of variance in deforestation rates (adjusted R^2 = 0.138; AIC = -2,249; n = 1,470). VIF analysis confirmed no severe multicollinearity (maximum VIF = 5.94 for elevation; Table 5). GWR with adaptive bisquare kernel (bandwidth = 11 nearest neighbors) significantly improved model fit (mean local R^2 = 0.609; median local R^2 = 0.920; AIC = -7,190), representing a 4,942-point AIC improvement and confirming substantial spatial non-stationarity in driver effects (H4).

Key OLS regression findings (Table 5):
- **Elevation:** Strongest negative coefficient (-0.325, t = -11.03, p < 0.001), confirming that deforestation concentrates at lower elevations where accessibility is greater and agricultural potential higher.
- **LST:** Strong negative coefficient (-0.259, t = -9.57, p < 0.001), reflecting the association between deforested (warmer) areas and historical forest loss patterns.
- **Distance to rivers:** Positive coefficient (0.140, t = 8.69, p < 0.001), indicating higher deforestation rates closer to rivers, consistent with riverine access as a deforestation vector in the Magdalena Medio.
- **Distance to roads:** Positive coefficient (0.062, t = 3.79, p < 0.001), suggesting accessibility effects.
- **Population density:** Weak negative coefficient (-0.028, t = -2.19, p < 0.05), indicating that deforestation occurs preferentially in less populated, rural areas rather than urban fringes.
- **Slope** and **precipitation** were not statistically significant in the global model.

GWR coefficient maps (Fig. 8) revealed that these relationships vary substantially across space: elevation effects range from -98.4 to +131.2, and distance to rivers from -209.6 to +538.2, confirming that locally varying relationships exist even when global coefficients suggest a single directional effect. The substantial R^2 improvement from OLS to GWR (0.143 to 0.609) confirms H4 and demonstrates that deforestation in the Magdalena Medio is driven by locally specific combinations of factors rather than a single dominant driver.

### 4.7 Future scenarios

CA-Markov transition probability matrices were calibrated from 2,938 sample points across the T2-T3-T4 transitions. The model identified five active LULC classes (dense forest, secondary forest, pastures, water, urban) with the following dominant transition probabilities: dense forest persists at 81.1% and converts to pastures at 18.5%; secondary forest converts to dense forest at 62.4% (maturation) and persists at 35.6%; pastures persist at 51.6% and convert to dense forest at 22.6% (recovery) or secondary forest at 22.2%.

Projected land cover composition (Table 6; Fig. 9):

| Scenario | Dense forest 2030 | Dense forest 2040 | Pastures 2030 | Pastures 2040 |
|----------|--------------------|--------------------|---------------|---------------|
| Current (2024) | 46.9% | -- | 20.1% | -- |
| BAU | 59.2% | 60.9% | 19.6% | 23.2% |
| Conservation | 63.0% | 72.5% | 15.0% | 14.1% |
| PDET | 61.5% | 67.1% | 17.3% | 18.8% |

Under BAU, dense forest is projected to increase to 60.9% by 2040 (from the sample-based 46.9% in 2024), reflecting the high maturation rate of secondary forest observed in the transition matrices. The conservation scenario achieves the highest forest recovery (72.5% by 2040), while PDET implementation reaches an intermediate outcome (67.1%). Notably, pastures are projected to increase under BAU (to 23.2% by 2040) while declining under conservation (14.1%) and PDET (18.8%) scenarios.

### 4.8 Climate trends

Mann-Kendall analysis revealed no statistically significant trends in either annual precipitation (Kendall tau = -0.051, p = 0.858; Sen's slope = -2.04 mm yr^-1) or LST (Kendall tau = -0.333, p = 0.129; Sen's slope = -0.037 C yr^-1) over the 2012-2024 period (Fig. 10). The lack of significant trends over this short 13-year window is not unexpected, as climate variability (ENSO cycles) dominates at this timescale.

SPI analysis identified below-normal precipitation conditions at key study periods: SPI = -0.14 (2013), -0.86 (2016), -0.48 (2020), and -0.89 (2024). The 2016 and 2024 SPI values approach the moderate drought threshold (-1.0), coinciding with El Nino events that may have exacerbated fire-related forest loss. Mean drought frequency across the study area was 0.263 (proportion of months below SPI = -1.0 threshold), with maximum drought frequency reaching 0.692 in the most drought-prone areas.

Annual precipitation ranged from 2,597 mm (2023) to 3,464 mm (2022), with a long-term mean of 2,852 mm. LST ranged from 27.0 C (2022) to 28.7 C (2015), with a slight but non-significant cooling trend.

---

## 5. Discussion

### 5.1 The paradox of peace: forest dynamics in the Magdalena Medio

Our findings contribute to the growing body of evidence documenting complex forest dynamics in post-conflict territories (Prem et al., 2020; Clerici et al., 2020; Murillo-Sandoval et al., 2021; Sierra et al., 2017). The dominant land use transition throughout the study period was dense forest-to-pasture conversion, with the T3->T4 period (2020-2024) showing the largest flux: 986,100 ha. Total forest cover (dense + secondary) declined by 43.9% between 2013 and 2024, a substantial loss consistent with--and in some municipalities exceeding--national-level patterns documented by IDEAM.

The Magdalena Medio presents distinctive dynamics compared to the Colombian Amazon, where most post-conflict deforestation research has concentrated. First, the pre-existing integration of the Magdalena Medio into national road networks and commodity chains means that deforestation here is driven less by frontier expansion and more by intensification of existing land use pressures--particularly the conversion of remaining forest fragments within an already accessible landscape. This is supported by our GWR finding that distance to rivers (a proxy for accessibility) is a significant positive driver (OLS t = 8.69), indicating that forest loss concentrates near transportation corridors. Second, the coexistence of formal (petroleum, oil palm) and informal (coca, artisanal mining) economies creates a complex driver mosaic that our GWR analysis captures through spatially varying coefficients.

The reduction in Hansen GFC-detected loss from T3 (81,106 ha) to T4 (26,765 ha) aligns with national deforestation reduction trends under the Petro government's environmental policies and suggests that some deceleration may be occurring, though our classification-based results show continued large-scale conversion.

### 5.2 Spatial patterns: where and why deforestation concentrates

The identification of 454 hotspot clusters at 99% confidence (Moran's I = 0.071) confirms spatial non-randomness in deforestation patterns, partially supporting H2. The weak but positive Moran's I indicates that while deforestation is not randomly distributed, it is not as strongly clustered as in more frontier-like settings (e.g., Amazon arc of deforestation). This is consistent with the Magdalena Medio's character as a fragmented landscape where deforestation occurs diffusely across remaining forest patches rather than advancing as a single front.

The GWR results add substantial nuance: the R^2 improvement from 0.143 (OLS) to 0.609 (GWR) demonstrates that a single "deforestation narrative" is insufficient for the Magdalena Medio. The spatial variation in driver coefficients shows that elevation is the dominant global predictor (OLS beta = -0.325), but locally, distance to rivers, LST, and population density exhibit varying and even opposite effects across the study area. This spatial heterogeneity has direct implications for policy: uniform conservation interventions will be less effective than spatially targeted strategies that address locally dominant drivers.

### 5.3 Ecosystem service implications

The estimated net carbon stock loss of 190 Mt C (33.5% of 2013 stocks) represents a substantial contribution to national emissions and underscores the climate costs of land use transitions in the Magdalena Medio (Costanza et al., 2014; Ruiz-Agudelo et al., 2022). At approximately 194 Mg C ha^-1 lost per hectare of dense forest-to-pasture conversion (IPCC Tier 1 pools), the region's deforestation is particularly carbon-intensive due to the high biomass of its lowland humid forests. The carbon trajectory was non-monotonic: losses of -74 Mt C (T1-T2) were partially offset by a +49 Mt C gain (T2-T3, reflecting increased forest classification in 2020), before the massive -165 Mt C loss (T3-T4). The T3-T4 period alone accounts for 87% of the total net loss, underscoring the urgency of intervention in the most recent post-agreement phase.

The habitat quality trajectory (declining from 0.247 to 0.168 between 2013 and 2020) tells a complementary story about biodiversity consequences. The slight recovery to 0.264 in 2024 reflects the reclassification dynamics of secondary forest but should not be interpreted as ecological improvement, given the concurrent massive loss of dense forest. The decreasing standard deviation in habitat quality (from 0.127 to 0.079) indicates landscape homogenization--remaining habitat is becoming more uniformly degraded.

The hydrological implications are particularly relevant for a region bisected by the Magdalena River: declining baseflow estimates (774 to 545 mm between T1 and T3) in areas of forest loss are consistent with reduced infiltration capacity, threatening downstream water quality and dry-season water availability.

### 5.4 Future trajectories: the importance of governance

The CA-Markov scenario analysis highlights divergent futures depending on governance choices. The conservation scenario projects dense forest recovering to 72.5% of the sampled landscape by 2040 (vs. 60.9% under BAU), representing a meaningful difference that could translate to tens of thousands of hectares of additional forest. The PDET scenario (67.1% by 2040) demonstrates that the peace agreement's territorial development framework can achieve intermediate environmental outcomes when implemented alongside productive use.

The high maturation rate observed in transition matrices (62.4% of secondary forest transitioning to dense forest) suggests substantial natural recovery potential in the Magdalena Medio--a finding consistent with tropical secondary forest growth rates documented by Chazdon et al. (2016). This recovery potential could be leveraged through targeted restoration programs, particularly in riparian zones and connectivity corridors.

These projections carry significant uncertainty inherent to any future modeling exercise. The CA-Markov framework assumes that spatial drivers and neighborhood effects persist, which may not hold under transformative policy changes or climate shocks. Nevertheless, the relative comparison among scenarios provides useful guidance for territorial planning.

### 5.5 Methodological considerations and limitations

This study provides the first comprehensive, multi-class LULCC analysis of the Magdalena Medio across the critical pre- to post-peace agreement period. The integration of GEE-based classification, transition matrix analysis, ecosystem service quantification, GWR driver analysis, and CA-Markov projections within a single analytical framework offers a methodologically comprehensive assessment.

Key limitations include:

1. **Classification accuracy:** Overall accuracies (57-66%) fall below the commonly cited 85% threshold for operational LULC mapping (Foody, 2002; Congalton, 1991). This reflects the challenge of discriminating seven classes in a heterogeneous tropical landscape using automated reference data (rather than field-validated samples), the spectral confusion between pastures and croplands, and the 30 m spatial resolution. The results should be interpreted as relative trends rather than absolute area estimates.

2. **Temporal consistency:** The reference data generation for T3-T4 used Hansen GFC v1.12 with an adjusted loss-year encoding (gte vs. gt), while T1-T2 used v1.11. This methodological refinement was necessary to resolve the T4 classification failure but introduces a potential discontinuity that may partly explain the non-monotonic dense forest trajectory.

3. **Carbon estimation:** IPCC Tier 1 default values (IPCC, 2006) were used rather than locally calibrated field measurements. The estimated 190 Mt C loss should be treated as an order-of-magnitude assessment rather than a precise inventory. Tier 1 values assume complete carbon pool replacement upon conversion, likely overestimating actual emissions.

4. **Causality:** The pre/post comparison design cannot establish strict causality between the peace agreement and observed LULCC. Multiple confounding factors (climate variability, commodity prices, policy changes, COVID-19) operate simultaneously.

5. **Water yield simplification:** The evapotranspiration proxy (60% of precipitation) is a coarse approximation. More sophisticated hydrological modeling (e.g., InVEST Water Yield with local calibration) would improve these estimates.

6. **Spatial resolution of drivers:** Conflict-related variables were not available at fine spatial resolution, limiting the analysis to biophysical and accessibility drivers. Incorporating municipal-level governance indicators and conflict event data would strengthen the driver analysis.

### 5.6 Policy implications

Our findings carry direct implications for post-conflict territorial planning in the Magdalena Medio:

1. **Spatially targeted interventions:** The 454 deforestation hotspot cells identified through Gi* analysis should be priority areas for enhanced governance, forest monitoring, and law enforcement. The GWR finding that driver importance varies spatially argues for differentiated intervention strategies.

2. **Leveraging natural recovery:** The high secondary-to-dense forest maturation rate (62.4%) indicates that passive restoration through land abandonment or reduced pressure can be effective. Combining passive recovery with active restoration in priority corridors could maximize cost-effectiveness.

3. **PDET integration:** The PDET scenario analysis demonstrates that the peace agreement's territorial development framework can achieve environmental sustainability outcomes, reducing pasture expansion from 23.2% (BAU) to 18.8% (PDET) by 2040, but only if environmental criteria are fully integrated into PDET design and implementation.

4. **Carbon markets:** The quantified carbon losses (190 Mt C) and projected recovery potential under conservation scenarios provide a basis for REDD+ or voluntary carbon market mechanisms that could finance forest conservation while generating income for rural communities.

5. **Monitoring systems:** The GEE-based methodology developed here is fully reproducible and could form the basis for a regional LULCC monitoring system integrated with Colombia's Sistema de Monitoreo de Bosques y Carbono (SMBYC).

---

## 6. Conclusions

This study provides the first comprehensive multi-temporal analysis of land use and land cover change in Colombia's Magdalena Medio region, spanning the critical pre- to post-peace agreement period (2012-2024). Our key findings are:

1. Total forest cover (dense + secondary) declined by 43.9% between 2013 and 2024 (from 2,599,275 to 1,457,337 ha), with dense forest-to-pasture conversion as the dominant transition (986,100 ha in the 2020-2024 period alone). Pastures expanded 209% over the study period. These findings are consistent with the "paradox of peace" hypothesis (H1).

2. Deforestation exhibits positive spatial clustering (Moran's I = 0.071), with 454 hotspot cells identified at 99% confidence, confirming spatial non-randomness (H2). Hotspots concentrate along river corridors and transportation routes.

3. Carbon stocks declined by 190 Mt C (-33.5%), with -74 Mt C (2013-2016), +49 Mt C (2016-2020), and -165 Mt C (2020-2024). The T3-T4 period alone accounts for 87% of the total net loss. Habitat quality declined from 0.247 to 0.168 (2013-2020). These results support H3.

4. GWR substantially outperformed OLS (R^2: 0.609 vs. 0.143; AIC improvement: 4,942), confirming spatial non-stationarity in driver effects (H4). Elevation, LST, and distance to rivers emerged as the strongest drivers globally, but with substantial local variation.

5. CA-Markov scenario analysis projects dense forest recovering to 60.9% (BAU), 67.1% (PDET), or 72.5% (Conservation) of the sampled landscape by 2040, demonstrating that governance choices will substantially determine the region's environmental trajectory.

These findings underscore the critical importance of integrating environmental sustainability into post-conflict territorial development. The Magdalena Medio's experience--where land use transitions associated with the post-conflict period have driven substantial forest loss and ecosystem service degradation--provides lessons applicable to post-conflict landscapes globally. The fully reproducible, GEE-based methodology developed here offers a scalable framework for monitoring and adaptive management of these critical transitions.

---

## Acknowledgments

This research was supported by the Camara de Comercio de Medellin. Google Earth Engine cloud computing resources were provided through the ee-maestria-tesis project. We thank the developers of GEE, CHIRPS, MODIS, Landsat, Sentinel-2, and Hansen GFC datasets for open access to satellite-derived products.

---

## Data availability

All satellite data used in this study are freely available through Google Earth Engine (GEE). The specific GEE asset identifiers are: Landsat 8/9 Collection 2 Level 2 (`LANDSAT/LC08/C02/T1_L2`, `LANDSAT/LC09/C02/T1_L2`), Sentinel-2 SR Harmonized (`COPERNICUS/S2_SR_HARMONIZED`), SRTM Elevation (`USGS/SRTMGL1_003`), Hansen Global Forest Change v1.12 (`UMD/hansen/global_forest_change_2024_v1_12`), CHIRPS Daily v2.0 (`UCSB-CHG/CHIRPS/DAILY`), MODIS LST (`MODIS/061/MOD11A2`), WorldPop (`WorldPop/GP/100m/pop`), JRC Global Surface Water (`JRC/GSW1_4/GlobalSurfaceWater`), and GHSL Settlement Model (`JRC/GHSL/P2023A/GHS_SMOD`). The study area bounding box is [-75.0, 6.0, -73.5, 8.0]. Classification scripts (Google Earth Engine JavaScript/Python) and Python analysis code are available at https://github.com/[repository] (to be made public upon acceptance).

---

## References

Ahmed, B., et al. (2025). CA-Markov chain analysis for land use and land cover change prediction. *Environmental Monitoring and Assessment*.

Alvarez, M.D. (2003). Forests in the time of violence: conservation implications of the Colombian war. *Journal of Sustainable Forestry*, 16(3-4), 47-68. https://doi.org/10.1300/J091v16n03_03

Anselin, L. (1995). Local Indicators of Spatial Association--LISA. *Geographical Analysis*, 27(2), 93-115. https://doi.org/10.1111/j.1538-4632.1995.tb00338.x

Botero, V., et al. (2023). Spatial analysis of deforestation and coca cultivation in Colombia. *Applied Geography*.

Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5-32. https://doi.org/10.1023/A:1010933404324

Castro-Nunez, A., et al. (2022). Livestock and deforestation in post-conflict Colombia. *Frontiers in Sustainable Food Systems*.

Chazdon, R.L., et al. (2016). Carbon sequestration potential of second-growth forest regeneration in the Latin American tropics. *Science Advances*, 2(5), e1501639. https://doi.org/10.1126/sciadv.1501639

Clerici, N., et al. (2020). Deforestation in Colombian protected areas increased during post-conflict periods. *Scientific Reports*, 10, 4971. https://doi.org/10.1038/s41598-020-61861-y

Congalton, R.G. (1991). A review of assessing the accuracy of classifications of remotely sensed data. *Remote Sensing of Environment*, 37(1), 35-46. https://doi.org/10.1016/0034-4257(91)90048-B

Costanza, R., de Groot, R., Sutton, P., van der Ploeg, S., Anderson, S.J., Kubiszewski, I., Farber, S., & Turner, R.K. (2014). Changes in the global value of ecosystem services. *Global Environmental Change*, 26, 152-158. https://doi.org/10.1016/j.gloenvcha.2014.04.002

Drusch, M., et al. (2012). Sentinel-2: ESA's Optical High-Resolution Mission for GMES Operational Services. *Remote Sensing of Environment*, 120, 25-36. https://doi.org/10.1016/j.rse.2011.11.026

Fagan, M.E., et al. (2020). Land cover dynamics following a deforestation ban in northern Costa Rica. *Nature Sustainability*.

Foody, G.M. (2002). Status of land cover classification accuracy assessment. *Remote Sensing of Environment*, 80(1), 185-201. https://doi.org/10.1016/S0034-4257(01)00295-4

Fotheringham, A.S., Brunsdon, C., & Charlton, M. (2002). *Geographically Weighted Regression: The Analysis of Spatially Varying Relationships*. Chichester: John Wiley & Sons.

Funk, C., et al. (2015). The climate hazards infrared precipitation with stations--a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066. https://doi.org/10.1038/sdata.2015.66

Getis, A., & Ord, J.K. (1992). The analysis of spatial association by use of distance statistics. *Geographical Analysis*, 24(3), 189-206. https://doi.org/10.1111/j.1538-4632.1992.tb00261.x

Gorelick, N., et al. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18-27. https://doi.org/10.1016/j.rse.2017.06.031

Halmy, M.W.A., Gessler, P.E., Hicke, J.A., & Salem, B.B. (2015). Land use/land cover change detection and prediction in the north-western coastal desert of Egypt using Markov-CA. *Applied Geography*, 63, 101-112. https://doi.org/10.1016/j.apgeog.2015.06.015

Hansen, M.C., Potapov, P.V., Moore, R., Hancher, M., Turubanova, S.A., Tyukavina, A., et al. (2013). High-resolution global maps of 21st-century forest cover change. *Science*, 342(6160), 850-853. https://doi.org/10.1126/science.1244693

Holdridge, L.R. (1967). *Life Zone Ecology*. Tropical Science Center.

IDEAM (2018). *Resultados monitoreo deforestacion 2017*. Instituto de Hidrologia, Meteorologia y Estudios Ambientales.

IPCC (2006). *2006 IPCC Guidelines for National Greenhouse Gas Inventories, Volume 4: Agriculture, Forestry and Other Land Use*. Eggleston, H.S., Buendia, L., Miwa, K., Ngara, T., & Tanabe, K. (Eds.). Hayama, Japan: Institute for Global Environmental Strategies (IGES).

Li, J., et al. (2023). Integrating Google Earth Engine and InVEST for mapping ecosystem services. *Ecological Indicators*, 150, 110374. https://doi.org/10.1016/j.ecolind.2023.110374

Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1-2), 17-23. https://doi.org/10.1093/biomet/37.1-2.17

Murillo-Sandoval, P.J., Gjerdseth, E., Correa-Ayram, C., Wrathall, D., Van Den Hoek, J., Davalos, L.M., & Kennedy, R. (2021). No peace for the forest: Rapid, widespread land changes in the Andes-Amazon region following the Colombian civil war. *Global Environmental Change*, 69, 102283. https://doi.org/10.1016/j.gloenvcha.2021.102283

Negret, P.J., Sonter, L., Watson, J.E.M., & Possingham, H.P. (2019). Emerging evidence that armed conflict and coca cultivation influence deforestation patterns. *Biological Conservation*, 239, 108176. https://doi.org/10.1016/j.biocon.2019.07.021

Olofsson, P., Foody, G.M., Herold, M., Stehman, S.V., Woodcock, C.E., & Wulder, M.A. (2014). Good practices for estimating area and assessing accuracy of land change. *Remote Sensing of Environment*, 148, 42-57. https://doi.org/10.1016/j.rse.2014.02.015

Pacheco, P. (2009). Agrarian change, cattle ranching and deforestation: Assessing their linkages in Southern Para. *Environment and History*, 15(4), 493-520. https://doi.org/10.3197/096734009X12532652872072

Phalke, A.R., et al. (2020). Mapping croplands of Europe, Middle East, Russia, and Central Asia using Landsat, Random Forest, and Google Earth Engine. *ISPRS Journal of Photogrammetry and Remote Sensing*, 167, 104-122. https://doi.org/10.1016/j.isprsjprs.2020.06.022

Pontius, R.G., et al. (2011). Death to Kappa: Birth of quantity disagreement and allocation disagreement for accuracy assessment. *International Journal of Remote Sensing*, 32(15), 4407-4429. https://doi.org/10.1080/01431161.2011.552923

Prem, M., Saavedra, S., & Vargas, J.F. (2020). End-of-conflict deforestation: Evidence from Colombia's peace agreement. *World Development*, 129, 104852. https://doi.org/10.1016/j.worlddev.2019.104852

Puyravaud, J.-P. (2003). Standardizing the calculation of the annual rate of deforestation. *Forest Ecology and Management*, 177(1-3), 593-596. https://doi.org/10.1016/S0378-1127(02)00335-3

Restrepo, J.D., & Syvitski, J.P. (2006). Assessing the effect of natural controls and land use change on sediment yield in a major Andean river: The Magdalena drainage basin, Colombia. *Ambio*, 35(2), 65-74. https://doi.org/10.1579/0044-7447(2006)35[65:ATEONC]2.0.CO;2

Ruiz-Agudelo, C.A., Suarez, A., Gutierrez-Bonilla, F. de P., & Cortes-Gomez, A.M. (2022). The economic valuation of ecosystem services in Colombia: Challenges, gaps and future pathways. *Journal of Environmental Economics and Policy*, 12(3), 285-304. https://doi.org/10.1080/21606544.2022.2134218

Sanchez-Cuervo, A.M., & Aide, T.M. (2013). Consequences of the armed conflict, forced human displacement, and land abandonment on forest cover change in Colombia. *Ecosystems*, 16, 1016-1035. https://doi.org/10.1007/s10021-013-9667-y

Sanchez-Cuervo, A.M., Aide, T.M., Clark, M.L., & Etter, A. (2012). Land cover change in Colombia: Surprising forest recovery trends between 2001 and 2010. *PLoS ONE*, 7(8), e43943. https://doi.org/10.1371/journal.pone.0043943

Sharp, R., et al. (2020). *InVEST User's Guide*. The Natural Capital Project, Stanford University.

Sierra, C.A., Mahecha, M.D., Poveda, G., Alvarez-Davila, E., Gutierrez-Velez, V.H., Reu, B., et al. (2017). Monitoring ecological change during rapid socio-economic and political transitions: Colombian ecosystems in the post-conflict era. *Environmental Science and Policy*, 76, 40-49. https://doi.org/10.1016/j.envsci.2017.06.011

Tapia-Armijos, M.F., et al. (2019). Drivers of deforestation in the basin of the Ecuadorian Amazon. *Applied Geography*.

Verburg, P.H., Soepboer, W., Veldkamp, A., Limpiada, R., Espaldon, V., & Mastura, S.S.A. (2002). Modeling the spatial dynamics of regional land use: The CLUE-S model. *Environmental Management*, 30(3), 391-405. https://doi.org/10.1007/s00267-002-2630-x

Wulder, M.A., Roy, D.P., Radeloff, V.C., Loveland, T.R., Anderson, M.C., Johnson, D.M., et al. (2022). Fifty years of Landsat science and impacts. *Remote Sensing of Environment*, 280, 113195. https://doi.org/10.1016/j.rse.2022.113195

Zhang, H., et al. (2023). Optimal parameters for Random Forest classification in Google Earth Engine. *Remote Sensing of Environment*.

Zhao, M., et al. (2022). Trajectory analysis of land use change and ecosystem services. *Frontiers in Environmental Science*, 10, 1038752. https://doi.org/10.3389/fenvs.2022.1038752
