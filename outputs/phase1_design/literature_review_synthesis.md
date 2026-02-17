# FASE 1.1: Revision de Literatura Sistematica
## Cambio de Uso del Suelo en el Magdalena Medio Post-Acuerdo de Paz

**Autor:** Cristian Espinal - Camara de Comercio de Medellin
**Fecha:** 2026-02-13
**Total referencias identificadas:** 25 papers principales

---

## 1. RESUMEN DE BUSQUEDA

Se realizaron 8 busquedas sistematicas cubriendo 4 ejes tematicos:
- Eje 1: LULCC post-conflicto en Colombia (6 papers clave)
- Eje 2: Servicios ecosistemicos y teledeteccion tropical (5 papers)
- Eje 3: Metodos de clasificacion y deteccion de cambios (6 papers)
- Eje 4: Estadistica espacial y modelamiento predictivo (5 papers)
- Transversales: Fragmentacion, hidrologia, clima (3 papers)

---

## 2. HALLAZGOS CLAVE POR EJE TEMATICO

### 2.1 Deforestacion Post-Conflicto en Colombia

**Hallazgo principal:** El acuerdo de paz de 2016 genero un incremento paradojico en la deforestacion.

- **Prem et al. (2020, World Development):** Las areas controladas por las FARC experimentaron un aumento diferencial en deforestacion despues del cese al fuego. El efecto se atenua con mayor presencia estatal y capacidad judicial.

- **Clerici et al. (2020, Scientific Reports):** 79% de las areas protegidas analizadas (31/39) mostraron aumento en deforestacion post-acuerdo: incremento del 177% en la tasa y 330 km2 de perdida forestal adicional.

- **Castro-Nunez et al. (2022, Frontiers):** La deforestacion se vincula a ganaderia extensiva y especulacion de tierras (no a agricultura de pequena escala). Incremento del 40% en conversion bosque-agricultura.

- **Fagan et al. (2020, Oregon State/Nature Sustainability):** La lenta implementacion de la gobernanza ambiental, mercados ilegales de tierra, y ganaderia ilicita son los mecanismos explicativos clave.

**DATO CRITICO PARA MAGDALENA MEDIO:** Sanchez-Cuervo & Aide (2013, Ecosphere) identificaron el Magdalena Medio como **hotspot de deforestacion** usando Getis-Ord Gi*, y notaron la **ausencia de areas protegidas** en la region.

### 2.2 Servicios Ecosistemicos y Teledeteccion

- **Li et al. (2023, Ecological Indicators):** Demostraron la integracion GEE + InVEST para cuantificar cambios en almacenamiento de carbono vinculados a transiciones LULC. Workflow directamente replicable.

- **Negrete-Cardoso et al. (2025, Regional Environmental Change):** Trade-offs entre servicios ecosistemicos en los Andes colombianos; patron espacial de sinergias/trade-offs varia por tipo de cobertura.

- **Zhao et al. (2022, Frontiers):** Analisis de trayectoria + InVEST captura mejor los cambios no lineales en carbono que comparaciones simples antes-despues.

- **Restrepo & Syvitski (2006, AMBIO):** La erosion en la cuenca del Magdalena aumento 33% entre 1972-2010 por deforestacion. Impactos directos en servicios hidrologicos.

### 2.3 Metodos de Clasificacion y Deteccion de Cambios

**Random Forest en GEE:**
- Phalke et al. (2020): OA > 84.3%, composites temporales de verano alcanzan 89.8%
- Zhang et al. (2023): ntree=500 optimo; composites de mediana superiores; multi-temporal mejora OA >90%
- Reyes-Palomeque et al. (2023): Fusion Sentinel-1+2 en Colombia alcanza 68-82% OA

**Deteccion de cambios:**
- Pinto et al. (2022): LandTrendr validado para bosque tropical colombiano (Amazonia)
- Schultz et al. (2021): BFAST con componente estacional critico para bosques tropicales
- Carrasco et al. (2019): Integracion multi-sensor (optico+SAR) mejora OA a 82-91%

### 2.4 Estadistica Espacial y Modelamiento

**Hotspots y autocorrelacion:**
- Sanchez-Cuervo (2013): Moran I + Getis-Ord identificaron Magdalena Medio como cluster de deforestacion (Moran I significativo, P < 0.00001)
- Botero et al. (2023): LISA confirmo autocorrelacion espacial positiva entre deforestacion y coca en 57 municipios High-High

**Drivers y prediccion:**
- Tapia-Armijos et al. (2019): GWR + Random Forest logro predicciones >74% para drivers de deforestacion en Amazonia ecuatoriana
- Ahmed et al. (2025): CA-Markov alcanzo OA 93.6%, Kappa 0.92 para escenarios futuros LULC

**Fragmentacion:**
- Taubert et al. (2025, Science): 58-80% de bosques tropicales se fragmentaron mas entre 2000-2020; agricultura itinerante causa 61% en tropicos

---

## 3. GAPS IDENTIFICADOS

### 3.1 Gaps Geograficos
| Gap | Detalle |
|-----|---------|
| **G1** | No existe estudio de LULCC completo (multi-clase) para el Magdalena Medio post-2016 |
| **G2** | Los estudios post-conflicto se enfocan en Amazonia y Pacifico, no en Magdalena Medio |
| **G3** | Magdalena Medio identificado como hotspot pero sin analisis detallado reciente |

### 3.2 Gaps Metodologicos
| Gap | Detalle |
|-----|---------|
| **G4** | Estudios post-conflicto usan Hansen GFC (solo deforestacion), no clasificacion LULC completa |
| **G5** | No se ha aplicado LandTrendr/BFAST al Magdalena Medio |
| **G6** | No existe integracion GEE + InVEST para servicios ecosistemicos en Magdalena Medio |
| **G7** | Analisis de fragmentacion con metricas de conectividad no aplicado en la region |

### 3.3 Gaps Tematicos
| Gap | Detalle |
|-----|---------|
| **G8** | No hay cuantificacion de perdida de carbono por deforestacion post-acuerdo en Magdalena Medio |
| **G9** | No se han modelado escenarios futuros LULC para la region |
| **G10** | Falta analisis integrado de drivers socioeconomicos + ambientales con GWR en Magdalena Medio |
| **G11** | No hay evaluacion de impactos hidrologicos de cambio LULC post-conflicto |

---

## 4. TABLA COMPARATIVA DE ENFOQUES METODOLOGICOS

| Estudio | Sensor | Clasificacion | Deteccion Cambio | Servicios Ecosist. | Estad. Espacial | Prediccion | Area |
|---------|--------|--------------|-------------------|--------------------|--------------------|------------|------|
| Prem 2020 | Hansen GFC | No | Deforest. binaria | No | DiD | No | Colombia FARC |
| Clerici 2020 | Hansen GFC | No | Deforest. binaria | No | Descriptiva | No | Areas protegidas |
| Castro-Nunez 2022 | Hansen GFC | No | Deforest. binaria | No | Panel data | No | Colombia |
| Sanchez-Cuervo 2013 | MODIS VCF | No | Woody cover | No | Moran I, Gi* | No | Colombia nacional |
| Phalke 2020 | Landsat 8 | RF en GEE | Post-clasificacion | No | No | No | Global |
| Li 2023 | GEE multi | RF en GEE | Trayectoria | InVEST Carbon | No | No | China |
| Tapia-Armijos 2019 | Multi | - | Deforest. | No | GWR+RF | No | Ecuador Amazon |
| Ahmed 2025 | Multi | Supervisada | Post-clasificacion | No | No | CA-Markov | Asia |
| **NUESTRO ESTUDIO** | **Landsat+Sentinel** | **RF en GEE** | **LandTrendr+matriz** | **InVEST (C,H2O,Hab)** | **Moran,GWR,Gi*** | **CA-Markov** | **Magdalena Medio** |

---

## 5. POSICIONAMIENTO DEL PAPER

### Contribucion unica:
1. **Primer estudio LULCC completo (multi-clase) del Magdalena Medio post-acuerdo de paz** - llenando G1, G2, G3
2. **Integracion metodologica sin precedentes** para la region: RF en GEE + LandTrendr + InVEST + GWR + CA-Markov (G4-G7)
3. **Cuantificacion de perdida de servicios ecosistemicos** (carbono, agua, habitat) directamente vinculada al post-conflicto (G8, G11)
4. **Escenarios predictivos 2030-2040** con implicaciones para politica publica (G9)
5. **Analisis espacialmente explicito de drivers** usando GWR, superando modelos globales (G10)

### Valor agregado:
- Reproducibilidad total (codigo GEE + Python abierto)
- Region critica poco estudiada con datos nuevos
- Alineacion con SDGs 13, 15, 16
- Relevancia para politica de paz y desarrollo territorial

---

## 6. REFERENCIAS CLAVE PARA CADA SECCION DEL PAPER

### Introduction (minimo 40 refs):
- Post-conflict landscapes: Prem 2020, Clerici 2020, Castro-Nunez 2022, Fagan 2020
- LULCC tropicales: Sanchez-Cuervo 2013, Sanchez-Cuervo 2012, Etter 2006
- Servicios ecosistemicos: Li 2023, Negrete-Cardoso 2025, Zhao 2022
- Metodologia RS: Phalke 2020, Zhang 2023, Pinto 2022
- Estadistica espacial: Tapia-Armijos 2019, Botero 2023

### Methods:
- Clasificacion RF: Phalke 2020, Zhang 2023, Reyes-Palomeque 2023
- Change detection: Pinto 2022, Schultz 2021, Carrasco 2019
- Servicios ecosistemicos: Li 2023, Zhao 2022
- Estadistica espacial: Tapia-Armijos 2019, Perez-Vega 2012
- Prediccion: Ahmed 2025

### Discussion:
- Contexto post-conflicto: Prem 2020, Clerici 2020, Castro-Nunez 2022
- Fragmentacion: Taubert 2025, Brinck 2017
- Impactos hidrologicos: Restrepo 2006, Restrepo 2012
- Drivers: Botero 2023, Clerici 2021

---

## 7. BUSQUEDAS COMPLEMENTARIAS RECOMENDADAS

Para completar las ~40 referencias necesarias, se recomienda buscar:
1. MapBiomas Colombia (datos LULC oficiales)
2. IDEAM reportes de deforestacion Colombia
3. Programa Bosques de Paz / Colombia Sostenible
4. Papers sobre indice de calidad de habitat InVEST en tropicos
5. Estudios de LST (Land Surface Temperature) en Colombia
6. Papers Mann-Kendall + CHIRPS precipitacion Colombia
7. Estudios de coca/conflicto en Magdalena Medio (sociologia)
8. Politicas publicas: PDET (Planes de Desarrollo con Enfoque Territorial)
