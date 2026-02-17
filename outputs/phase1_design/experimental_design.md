# FASE 1.3: Diseno Experimental
## Paper: LULCC y Servicios Ecosistemicos en el Magdalena Medio Post-Acuerdo de Paz

**Autor:** Cristian Espinal - Camara de Comercio de Medellin
**Fecha:** 2026-02-13

---

## 1. AREA DE ESTUDIO

### 1.1 Definicion Geografica

**Region:** Magdalena Medio, Colombia
**Extension aproximada:** ~30,000 km2
**Coordenadas generales:** 6.0N - 8.0N latitud, 73.5W - 75.0W longitud
**Eje geografico:** Valle medio del Rio Magdalena, entre las cordilleras Central y Oriental

### 1.2 Municipios Incluidos (30 municipios nucleares)

| Departamento | Municipios | # |
|-------------|-----------|---|
| **Santander** | Barrancabermeja, Puerto Wilches, Sabana de Torres, Cimitarra, Puerto Parra, Rionegro, San Vicente de Chucuri, El Carmen de Chucuri, Betulia, Bolivar, Landazuri, Simacota, Santa Helena del Opon, El Penon | 14 |
| **Antioquia** | Puerto Berrio, Yondo, Puerto Nare, Puerto Triunfo, Caracoli, Maceo | 6 |
| **Bolivar** | San Pablo, Simiti, Santa Rosa del Sur, Cantagallo, Morales, Arenal | 6 |
| **Cesar** | Aguachica, Gamarra, San Martin, San Alberto | 4 |

**Capital regional de facto:** Barrancabermeja (Santander)

### 1.3 Caracteristicas Biofisicas

| Parametro | Valor |
|-----------|-------|
| Altitud | 50 - 500 msnm (predominante <200 m) |
| Temperatura media | 27 - 30 C |
| Precipitacion anual | 1,800 - 3,500 mm |
| Regimen de lluvias | Bimodal (abr-may, sep-nov) |
| Periodo seco | Dic-Mar, Jun-Ago (relativo) |
| Ecoregion | Bosque humedo tropical del Magdalena-Uraba |
| Cuenca principal | Rio Magdalena (tramo medio) |
| Suelos predominantes | Aluviales, con sectores de alta fertilidad |

### 1.4 Caracteristicas Socioeconomicas

| Parametro | Valor |
|-----------|-------|
| Poblacion aproximada | ~800,000 habitantes (30 municipios) |
| Urbanizacion | ~56% urbana, 44% rural |
| Economia | Petroleo (Barrancabermeja), ganaderia extensiva, mineria aurifera, palma de aceite, coca |
| Conflicto armado | Region historica de presencia FARC, ELN, paramilitares |
| PDET | Si, municipios priorizados en acuerdo de paz |
| Infraestructura | Ruta del Sol, rio Magdalena navegable, refineria Ecopetrol |

### 1.5 Justificacion de la Seleccion del Area

1. **Hotspot de deforestacion confirmado** por Sanchez-Cuervo & Aide (2013) usando Getis-Ord Gi*
2. **Ausencia de areas protegidas** significativas en la region (gap de conservacion)
3. **Region PDET** priorizada en el acuerdo de paz
4. **Transicion economica activa** (petroleo + ganaderia + palma + mineria)
5. **Gap de conocimiento critico**: no existe estudio LULCC multi-clase post-2016
6. **Diversidad de drivers**: multiples presiones simultaneas (ganaderia, coca, mineria, palma)

---

## 2. DISENO TEMPORAL

### 2.1 Periodos de Analisis

```
    2012        2016        2020        2024
      |-----------|-----------|-----------|
      PRE-ACUERDO  POST-1      POST-2
      (conflicto)  (transicion) (implementacion)
```

| Periodo | Anios | Justificacion | Contexto politico |
|---------|-------|---------------|-------------------|
| **T1: Pre-acuerdo** | 2012-2015 | Linea base durante conflicto activo | Negociaciones en La Habana (2012-2016); FARC controlan territorios |
| **T2: Transicion** | 2016-2019 | Periodo inmediato post-cese al fuego | Firma del acuerdo (nov 2016); desmovilizacion FARC; vacios de gobernanza |
| **T3: Implementacion** | 2020-2024 | Consolidacion e implementacion PDET | COVID-19 (2020); gobierno Petro (2022); Paz Total |

### 2.2 Mapas de Clasificacion LULC (4 cortes temporales)

| Mapa | Anio central | Ventana de imagenes | Sensor principal | Sensor complementario |
|------|-------------|---------------------|------------------|----------------------|
| **M1** | 2013 | 2012-01-01 a 2014-12-31 | Landsat 8 OLI | - |
| **M2** | 2016 | 2015-01-01 a 2017-06-30 | Landsat 8 OLI | Sentinel-2A (desde jun 2015) |
| **M3** | 2020 | 2019-01-01 a 2021-06-30 | Landsat 8 OLI | Sentinel-2A/B |
| **M4** | 2024 | 2023-01-01 a 2024-12-31 | Landsat 8/9 OLI-2 | Sentinel-2A/B |

**Nota:** Se usa ventana de ~2 anios para maximizar pixeles libres de nubes en zona tropical.

### 2.3 Series Temporales Continuas

| Analisis | Periodo | Frecuencia | Fuente |
|----------|---------|------------|--------|
| NDVI/EVI trend | 2012-2024 | Mensual | MODIS MOD13Q1 (250m) |
| Precipitacion | 2012-2024 | Mensual | CHIRPS (5km) |
| Temperatura superficial | 2012-2024 | 8-dias | MODIS MOD11A2 (1km) |
| Forest loss | 2012-2023 | Anual | Hansen GFC v1.11 (30m) |
| Alertas deforestacion | 2016-2024 | Trimestral | IDEAM AT-D |

---

## 3. VARIABLES DEL ESTUDIO

### 3.1 Variable Dependiente Principal

**Cobertura y uso del suelo (LULC)** - 7 clases:

| Clase | Codigo | Descripcion | Fuente validacion |
|-------|--------|-------------|-------------------|
| Bosque denso | 1 | Cobertura arborea >60%, dosel cerrado | MapBiomas, Hansen GFC |
| Bosque fragmentado/secundario | 2 | Cobertura arborea 30-60%, sucesion | MapBiomas |
| Pasturas/Herbazales | 3 | Ganaderia extensiva, pastos naturales e introducidos | MapBiomas, UPRA |
| Cultivos | 4 | Agricultura (palma, arroz, coca, otros) | MapBiomas, UPRA |
| Cuerpos de agua | 5 | Rios, cienagas, humedales, embalses | JRC Water, HydroSHEDS |
| Areas urbanas/Infraestructura | 6 | Centros poblados, vias, zonas industriales | WorldPop, GHSL |
| Suelo desnudo/Mineria | 7 | Areas expuestas, mineria, degradacion | Visual + indices espectrales |

### 3.2 Variables Dependientes Derivadas

| Variable | Unidad | Derivada de | Metodo |
|----------|--------|------------|--------|
| Tasa de deforestacion | ha/anio, %/anio | Matrices de transicion | Post-clasificacion |
| Tasa de cambio neto por clase | ha/anio | Matrices de transicion | Post-clasificacion |
| Stock de carbono | Mg C / ha | LULC + pools de carbono | InVEST Carbon |
| Rendimiento hidrico | mm/anio | LULC + clima + suelo | InVEST SWY |
| Calidad de habitat | Indice 0-1 | LULC + amenazas | InVEST HQ |
| Numero de parches (NP) | count | LULC bosque | FRAGSTATS |
| Densidad de bordes (ED) | m/ha | LULC bosque | FRAGSTATS |
| Indice de agregacion (AI) | % | LULC bosque | FRAGSTATS |
| Indice de cohesion (COHESION) | % | LULC bosque | FRAGSTATS |

### 3.3 Variables Independientes (Drivers)

| Variable | Tipo | Resolucion | Fuente GEE | Asset ID o Coleccion |
|----------|------|-----------|------------|---------------------|
| Elevacion (DEM) | Continua | 30 m | SRTM | USGS/SRTMGL1_003 |
| Pendiente | Continua | 30 m | Derivada SRTM | ee.Terrain.slope() |
| Distancia a carreteras | Continua | - | OpenStreetMap / GRIP | projects/sat-io/open-datasets/grip-roads |
| Distancia a rios | Continua | - | HydroSHEDS | WWF/HydroSHEDS |
| Distancia a centros urbanos | Continua | - | GHSL | JRC/GHSL/P2023A/GHS_SMOD |
| Precipitacion media anual | Continua | 5 km | CHIRPS | UCSB-CHG/CHIRPS/DAILY |
| Temperatura superficial media | Continua | 1 km | MODIS LST | MODIS/061/MOD11A2 |
| Densidad poblacional | Continua | 100 m | WorldPop | WorldPop/GP/100m/pop |
| Indice de aridez | Continua | 1 km | CGIAR-CSI | - |
| Aptitud ganadera | Categorica | 1 km | FAO GAEZ / UPRA | - |
| Presencia historica FARC | Categorica | Municipal | CERAC / datos externos | Shapefile externo |
| Distancia a areas protegidas | Continua | - | WDPA | WCMC/WDPA |

### 3.4 Variables de Control

| Control | Justificacion | Implementacion |
|---------|---------------|----------------|
| **Topografia** | Controlar efecto de altitud/pendiente sobre cobertura | Incluir DEM y pendiente en modelos |
| **Accesibilidad** | Controlar sesgo por distancia a infraestructura | Distancia a vias y centros urbanos |
| **Clima** | Controlar variabilidad climatica interanual | Series CHIRPS + analisis de anomalias |
| **COVID-19** | Controlar efecto pandemia 2020-2021 sobre dinamicas | Analisis sub-periodo T2a (2016-2019) vs T2b (2020-2021) |
| **Tipo de suelo** | Controlar aptitud agricola intrinseca | Datos de suelos IGAC / SoilGrids |
| **MapBiomas** | Validacion cruzada independiente | Comparar clasificacion propia vs MapBiomas Colombia Collection 2.0 |

---

## 4. FUENTES DE DATOS SATELITALES

### 4.1 Colecciones GEE Principales

| # | Coleccion | Asset ID GEE | Bandas/Variables | Resolucion | Uso |
|---|-----------|-------------|------------------|-----------|-----|
| 1 | Landsat 8 SR | LANDSAT/LC08/C02/T1_L2 | B2-B7, QA_PIXEL | 30 m | Clasificacion LULC |
| 2 | Landsat 9 SR | LANDSAT/LC09/C02/T1_L2 | B2-B7, QA_PIXEL | 30 m | Clasificacion LULC (2022+) |
| 3 | Sentinel-2 SR | COPERNICUS/S2_SR_HARMONIZED | B2-B12, SCL | 10-20 m | Clasificacion LULC (2016+) |
| 4 | MODIS NDVI | MODIS/061/MOD13Q1 | NDVI, EVI | 250 m | Series temporales vegetacion |
| 5 | MODIS LST | MODIS/061/MOD11A2 | LST_Day, LST_Night | 1 km | Temperatura superficial |
| 6 | CHIRPS Daily | UCSB-CHG/CHIRPS/DAILY | precipitation | 5 km | Precipitacion |
| 7 | ERA5 Monthly | ECMWF/ERA5_LAND/MONTHLY_AGGR | temperature_2m, evaporation | 11 km | Variables climaticas |
| 8 | SRTM DEM | USGS/SRTMGL1_003 | elevation | 30 m | Topografia |
| 9 | Hansen GFC | UMD/hansen/global_forest_change_2023_v1_11 | treecover, loss, gain | 30 m | Validacion forestal |
| 10 | JRC Water | JRC/GSW1_4/GlobalSurfaceWater | occurrence, recurrence | 30 m | Cuerpos de agua |
| 11 | WorldPop | WorldPop/GP/100m/pop | population | 100 m | Densidad poblacional |
| 12 | GHSL | JRC/GHSL/P2023A/GHS_SMOD | smod | 1 km | Asentamientos humanos |
| 13 | MapBiomas | projects/mapbiomas_af_trinacional/public/collection1/mapbiomas_colombia_collection1_integration_v1 | classification | 30 m | Validacion cruzada |
| 14 | SoilGrids | projects/soilgrids-isric | clay, sand, soc, ph | 250 m | Propiedades del suelo |

### 4.2 Indices Espectrales a Calcular

| Indice | Formula | Uso |
|--------|---------|-----|
| NDVI | (NIR - Red) / (NIR + Red) | Vigor vegetacion |
| EVI | 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1) | Vegetacion densa |
| NDWI | (Green - NIR) / (Green + NIR) | Cuerpos de agua |
| NDBI | (SWIR1 - NIR) / (SWIR1 + NIR) | Areas urbanas |
| BSI | ((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue)) | Suelo desnudo |
| NBR | (NIR - SWIR2) / (NIR + SWIR2) | Areas quemadas/degradadas |
| SAVI | ((NIR - Red) / (NIR + Red + 0.5)) * 1.5 | Vegetacion con suelo expuesto |
| MNDWI | (Green - SWIR1) / (Green + SWIR1) | Agua mejorado |

### 4.3 Estrategia de Compositing

```
Para cada periodo temporal:

1. Filtrar coleccion por ROI y rango de fechas
2. Cloud masking:
   - Landsat: QA_PIXEL bit flags (cloud, cloud shadow, cirrus)
   - Sentinel-2: SCL band (clases 3,8,9,10,11 = nubes/sombras)
3. Calcular indices espectrales para cada imagen
4. Generar composite de MEDIANA por pixel
5. Si gaps persisten: ampliar ventana temporal +6 meses
6. Bandas del composite final:
   - 6 bandas reflectancia (B, G, R, NIR, SWIR1, SWIR2)
   - 8 indices espectrales
   - 3 variables topograficas (elev, slope, aspect)
   - Total: 17 features para clasificacion RF
```

---

## 5. PROTOCOLO DE CLASIFICACION

### 5.1 Muestreo de Entrenamiento

| Parametro | Valor |
|-----------|-------|
| Metodo | Muestreo estratificado aleatorio |
| Minimo por clase | 500 puntos (total minimo: 3,500) |
| Fuentes de referencia | Google Earth Pro (alta resolucion), MapBiomas, datos IDEAM |
| Distribucion espacial | Proporcional al area de cada clase + sobremuestreo de clases raras |
| Consistencia temporal | Mismos sitios re-evaluados para cada periodo cuando sea posible |

### 5.2 Algoritmo de Clasificacion

| Parametro RF | Valor |
|-------------|-------|
| Algoritmo | ee.Classifier.smileRandomForest() |
| Numero de arboles (ntree) | 500 |
| Variables por split (mtry) | sqrt(n_features) = ~4 |
| Min samples leaf | 5 |
| Features de entrada | 17 (6 bandas + 8 indices + 3 topo) |
| Modo | Clasificacion supervisada |

### 5.3 Validacion y Accuracy Assessment

| Metodo | Descripcion |
|--------|-------------|
| **Hold-out** | 70% entrenamiento / 30% validacion (estratificado) |
| **Metricas** | Overall Accuracy, Kappa, F1-score por clase, User's/Producer's Accuracy |
| **Cross-validation temporal** | Entrenar con T1+T2, validar con T3 (y viceversa) |
| **Cross-validation espacial** | K-fold espacial (5 bloques geograficos) |
| **Validacion independiente** | Comparar con MapBiomas Colombia Collection 2.0 (mismos anios) |
| **Validacion con Hansen** | Comparar clase bosque vs. Hansen treecover2000 + loss |
| **Incertidumbre** | Mapa de probabilidad de clase (RF class probabilities) |

**Criterios de aceptacion:**
- OA >= 85%
- Kappa >= 0.80
- F1 por clase >= 0.75 (todas las clases)
- Concordancia con MapBiomas >= 80% (nivel 1)

---

## 6. PROTOCOLO DE ANALISIS DE CAMBIO

### 6.1 Matrices de Transicion

```
3 matrices de transicion:
  T1->T2: 2013 a 2016 (pre -> transicion)
  T2->T3: 2016 a 2020 (transicion -> post-1)
  T3->T4: 2020 a 2024 (post-1 -> post-2)

Cada matriz: 7x7 clases
Unidades: hectareas y porcentaje
Estadisticas: persistencia, ganancias, perdidas, cambio neto, swap
```

### 6.2 Tasas de Cambio

| Metrica | Formula | Referencia |
|---------|---------|-----------|
| Tasa anual de deforestacion | r = (1/t) * ln(A2/A1) | FAO (Puyravaud 2003) |
| Tasa neta de cambio por clase | (A_t2 - A_t1) / A_t1 * 100 | Standard |
| Intensidad de transicion | Pontius et al. 2004 framework | Intensity Analysis |

### 6.3 Deteccion Continua (LandTrendr)

| Parametro LandTrendr | Valor |
|---------------------|-------|
| maxSegments | 6 |
| spikeThreshold | 0.9 |
| vertexCountOvershoot | 3 |
| preventOneYearRecovery | true |
| recoveryThreshold | 0.25 |
| pvalThreshold | 0.05 |
| bestModelProportion | 0.75 |
| minObservationsNeeded | 6 |
| Indice | NBR (Normalized Burn Ratio) |
| Periodo | 2012-2024 |

---

## 7. PROTOCOLO DE SERVICIOS ECOSISTEMICOS

### 7.1 InVEST Carbon Storage

| Parametro | Descripcion |
|-----------|-------------|
| Modelo | Carbon Storage and Sequestration |
| Inputs | Mapa LULC + tabla de pools de carbono por clase |
| Pools de carbono | C_above (biomasa aerea), C_below (raices), C_soil (suelo), C_dead (necromasa) |
| Fuente de valores | IPCC Tier 1 defaults para bosque tropical humedo + literatura regional |

**Valores de referencia (Mg C/ha):**

| Clase LULC | C_above | C_below | C_soil | C_dead | Total |
|-----------|---------|---------|--------|--------|-------|
| Bosque denso | 120 | 30 | 80 | 12 | 242 |
| Bosque secundario | 60 | 15 | 65 | 6 | 146 |
| Pasturas | 5 | 3 | 40 | 0.5 | 48.5 |
| Cultivos | 8 | 2 | 35 | 0.5 | 45.5 |
| Agua | 0 | 0 | 0 | 0 | 0 |
| Urbano | 2 | 0 | 20 | 0 | 22 |
| Suelo desnudo | 0 | 0 | 15 | 0 | 15 |

*Nota: Valores IPCC Tier 1 para bosque tropical humedo; se ajustaran con literatura regional.*

### 7.2 InVEST Seasonal Water Yield

| Parametro | Fuente |
|-----------|--------|
| Precipitacion | CHIRPS mensual |
| Evapotranspiracion referencia | MODIS MOD16A2 |
| Tipo de suelo / Ksat | SoilGrids ISRIC |
| Curva numero | Valores por clase LULC + grupo hidrologico |
| DEM | SRTM 30m |

### 7.3 InVEST Habitat Quality

| Parametro | Valor |
|-----------|-------|
| Amenazas | Pasturas (dist_max=8km), Cultivos (dist_max=6km), Urbano (dist_max=10km), Carreteras (dist_max=3km), Mineria (dist_max=5km) |
| Sensibilidad por habitat | Bosque denso: max sensibilidad; Bosque secundario: alta; otros: baja |
| Half-saturation constant | 0.5 |
| Decay function | Exponential |

---

## 8. PROTOCOLO DE ANALISIS ESTADISTICO

### 8.1 Comparacion Temporal

| Test | Aplicacion | Datos |
|------|-----------|-------|
| Mann-Whitney U | Comparar tasas de cambio pre vs. post | Tasas municipales |
| Chi-cuadrado | Comparar matrices de transicion entre periodos | Matrices de transicion |
| t-test pareado | Comparar ES pre vs. post por unidad espacial | Valores ES por subcuenca/municipio |
| ANOVA + Tukey HSD | Comparar entre los 3 periodos | Variables ES y cambio |

### 8.2 Analisis Espacial

| Metodo | Parametros | Software |
|--------|-----------|---------|
| Moran's I global | Queen contiguity weight, 999 permutaciones | PySAL / GeoDa |
| Getis-Ord Gi* | Fixed distance band, FDR correction | PySAL / ArcPy |
| Kernel density | Bandwidth = 5km, Gaussian kernel | scipy.stats |
| FRAGSTATS | 8-cell neighborhood, class + landscape level | pylandstats / FRAGSTATS |

### 8.3 Modelamiento de Drivers (GWR)

| Parametro | Valor |
|-----------|-------|
| Modelo base | OLS global (referencia) |
| Modelo espacial | GWR con kernel adaptativo |
| Seleccion bandwidth | AICc minimization |
| Variable dependiente | Tasa de deforestacion 2016-2024 (%) |
| Variables independientes | 10 drivers (ver seccion 3.3) |
| Seleccion de variables | Stepwise + VIF < 7.5 |
| Diagnostico | AICc, R2 local, residuales Moran |
| Comparacion | AICc(GWR) vs AICc(OLS) |

### 8.4 Modelamiento Predictivo (CA-Markov)

| Parametro | Valor |
|-----------|-------|
| Calibracion | LULC 2013 + 2016 -> prediccion 2020 |
| Validacion | Prediccion 2020 vs. LULC 2020 observado |
| Proyeccion | 2030, 2040 |
| Escenarios | 3: BAU, Conservacion, Desarrollo PDET |
| Iteraciones CA | 10 por anio |
| Factores de restriccion | Areas protegidas, cuerpos de agua, pendiente >45 grados |
| Software | MOLUSCE plugin o TerrSet / Python CA-Markov |

**Escenarios definidos:**

| Escenario | Reglas de transicion |
|-----------|---------------------|
| **BAU (Business as Usual)** | Probabilidades de transicion 2016-2024 se mantienen constantes |
| **Conservacion** | Deforestacion reducida 50%; restauracion bosque secundario +30% |
| **PDET (Desarrollo planificado)** | Expansion agricola controlada; restauracion de bosque ripario; ganaderia intensificada (menos area) |

---

## 9. UNIDADES ESPACIALES DE ANALISIS

| Nivel | Unidad | Cantidad aprox. | Uso |
|-------|--------|----------------|-----|
| 1 - Region | Magdalena Medio completo | 1 | Estadisticas generales |
| 2 - Departamento | Santander, Antioquia, Bolivar, Cesar | 4 | Comparacion departamental |
| 3 - Municipio | 30 municipios | 30 | Analisis de drivers, GWR |
| 4 - Subcuenca | Cuencas tributarias Magdalena | ~15-20 | Analisis hidrologico |
| 5 - Grilla | Celdas regulares 1x1 km | ~30,000 | Hotspot analysis, densidad |
| 6 - Pixel | 30x30 m | ~33 millones | Clasificacion, cambio |

---

## 10. CRONOGRAMA DE EJECUCION

| Fase | Actividad | Duracion estimada | Dependencias |
|------|-----------|-------------------|-------------|
| 2.1 | Adquisicion y preprocesamiento GEE | 3-4 dias | Config GEE lista |
| 2.2 | Muestreo de entrenamiento | 2-3 dias | Composites listos |
| 2.3 | Clasificacion RF + validacion | 2-3 dias | Training data + composites |
| 2.4 | Matrices de transicion + LandTrendr | 2-3 dias | Mapas LULC clasificados |
| 2.5 | FRAGSTATS + hotspots | 1-2 dias | Mapas LULC |
| 2.6 | InVEST (Carbon, Water, Habitat) | 2-3 dias | Mapas LULC + datos auxiliares |
| 2.7 | Analisis climatico | 1-2 dias | CHIRPS + MODIS procesados |
| 3.1 | Analisis estadistico espacial | 2-3 dias | Todas las variables |
| 3.2 | GWR + drivers | 1-2 dias | Variables dependientes + independientes |
| 3.3 | CA-Markov escenarios | 2-3 dias | Matrices de transicion + drivers |
| 4.1 | Figuras y mapas | 2-3 dias | Todos los resultados |
| 5.1 | Redaccion paper | 5-7 dias | Todos los analisis + figuras |
| 6-8 | Validacion + submission prep | 3-4 dias | Manuscrito completo |

**Total estimado: 26-38 dias de trabajo**

---

## 11. ESTRUCTURA DE ARCHIVOS DEL PROYECTO

```
magdalena_medio_research/
|-- .env                          # GEE Project ID
|-- gee_config.py                 # Configuracion GEE (ACTUALIZAR bbox)
|-- requirements.txt              # Dependencias Python
|
|-- scripts/
|   |-- 01_preprocessing.py       # Cloud masking, composites
|   |-- 02_training_samples.py    # Generacion muestras entrenamiento
|   |-- 03_classification.py      # RF classification en GEE
|   |-- 04_accuracy_assessment.py # Validacion y metricas
|   |-- 05_change_detection.py    # Matrices transicion + LandTrendr
|   |-- 06_fragmentation.py       # FRAGSTATS metricas
|   |-- 07_hotspot_analysis.py    # Moran I + Getis-Ord
|   |-- 08_ecosystem_services.py  # InVEST models
|   |-- 09_climate_analysis.py    # CHIRPS + MODIS trends
|   |-- 10_gwr_drivers.py         # Regresion geograficamente ponderada
|   |-- 11_ca_markov.py           # Escenarios futuros
|   |-- 12_visualization.py       # Figuras publication-ready
|   |-- utils.py                  # Funciones auxiliares
|
|-- data/
|   |-- literature/               # Base de datos bibliografica
|   |-- shapefiles/               # Limites municipales, subcuencas
|   |-- training_samples/         # Puntos de entrenamiento/validacion
|   |-- rasters/                  # Datos raster procesados
|   |-- tables/                   # Datos tabulares
|
|-- outputs/
|   |-- phase1_design/            # Documentos de diseno
|   |-- phase2_gee/               # Resultados GEE
|   |-- phase3_stats/             # Resultados estadisticos
|   |-- phase4_figures/           # Figuras finales
|
|-- figures/                      # Figuras publication-ready (300 DPI)
|-- tables/                       # Tablas del paper
|-- supplementary/                # Material suplementario
|-- logs/                         # Logs de procesamiento
```

---

## 12. ACTUALIZACION REQUERIDA: gee_config.py

El bounding box actual es demasiado amplio. Se debe actualizar:

**Actual:** [-75.5, 5.5] a [-73.5, 8.5] (~120,000 km2)
**Propuesto:** Usar limites municipales reales de los 30 municipios

La estrategia optima es:
1. Cargar shapefile de municipios DANE de Colombia
2. Filtrar los 30 municipios del Magdalena Medio
3. Usar la union (dissolve) como area de estudio
4. Alternativa: bbox ajustado [-75.0, 6.0] a [-73.5, 8.0] + clip con municipios

---

## 13. CONSIDERACIONES ETICAS Y LIMITACIONES

### Consideraciones Eticas
- Datos satelitales de acceso abierto (no requieren permisos especiales)
- No se recolectan datos de individuos
- Datos de conflicto armado usados de forma agregada (municipal)
- Resultados deben considerar sensibilidad politica de la region

### Limitaciones Anticipadas
1. **Nubosidad:** Zona tropical con alta cobertura nubosa; mitigado con composites multi-anio y multi-sensor
2. **Resolucion temporal:** 4 cortes temporales pueden perder cambios intra-periodo; mitigado con LandTrendr continuo
3. **Resolucion espacial:** 30m puede no detectar cambios pequenos (<0.1 ha); mitigado con Sentinel-2 a 10m
4. **Ground truth:** Acceso limitado a datos de campo; mitigado con validacion cruzada MapBiomas + Hansen + Google Earth
5. **Causalidad:** Analisis correlacional, no causal; mitigado con diseno cuasi-experimental pre/post
6. **Escala:** Analisis regional, no local; resultados no extrapolables directamente a parcela
7. **Datos de conflicto:** Variables proxy (presencia municipal FARC), no mediciones directas de control territorial
