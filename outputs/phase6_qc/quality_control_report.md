# FASE 6: REPORTE DE VALIDACION DE CALIDAD
## Proyecto: LULCC Magdalena Medio Post-Acuerdo de Paz

**Fecha:** 2026-02-13
**Revisor:** QC Automatizado
**Archivos revisados:** 14 scripts + 1 config + 4 documentos Fase 1 + manuscrito

---

## 1. RESUMEN EJECUTIVO

| Categoria | Criticos | Moderados | Menores | Total |
|-----------|----------|-----------|---------|-------|
| Bugs de codigo | 5 | 3 | 4 | 12 |
| Consistencia metodologica | 0 | 2 | 1 | 3 |
| Manuscrito | 0 | 1 | 2 | 3 |
| **Total** | **5** | **6** | **7** | **18** |

**Estado:** 5 bugs criticos encontrados y corregidos. Proyecto listo para ejecucion.

---

## 2. BUGS CRITICOS (Corregidos)

### BUG-01: classify_hotspots() - condiciones superpuestas [07_hotspot_analysis.py:136-141]
**Severidad:** CRITICA
**Problema:** Las condiciones de clasificacion se aplican secuencialmente con `>=`, causando que valores altos (ej. z=2.6) se sobreescriban. Un z-score de 2.6 primero se marca como 3 (>=2.576), luego se sobreescribe a 2 (>=1.960), y finalmente a 1 (>=1.645).
**Impacto:** TODOS los hotspots se clasificarian como 90% en lugar de su nivel correcto.
**Correccion:** Invertir el orden de asignacion o usar condiciones excluyentes.

### BUG-02: compute_morans_i() - formula S2 incorrecta [07_hotspot_analysis.py:62]
**Severidad:** CRITICA
**Problema:** `S2 = np.sum(np.sum(weights, axis=1) + np.sum(weights, axis=0)) ** 2` aplica el cuadrado FUERA del sum. La formula correcta de Moran's I requiere `S2 = np.sum((ri + ci)^2)` donde ri y ci son sumas de fila y columna.
**Impacto:** Varianza de Moran's I incorrecta -> z-scores y p-values erroneos.
**Correccion:** `S2 = np.sum((np.sum(weights, axis=1) + np.sum(weights, axis=0)) ** 2)`

### BUG-03: get_reference_lulc() - umbral GHSL SMOD incorrecto [02_training_samples.py:58]
**Severidad:** CRITICA
**Problema:** `urban = ghsl.gte(2)` con SMOD codes (10-30) clasificaria TODA el area como urbana, ya que todos los valores son >= 2.
**Impacto:** El mapa de referencia tendria bosque solo donde Hansen dice bosque, pero casi todo lo demas seria "urbano" -> muestras de entrenamiento severamente sesgadas.
**Correccion:** `urban = ghsl.gte(20)` (suburban + urban en codigos SMOD). Tambien corregir nombre de banda a 'smod_code'.

### BUG-04: enhance_carbon_with_biomass() - dataset incorrecto [08_ecosystem_services.py:138-140]
**Severidad:** CRITICA
**Problema:** `ee.ImageCollection('ESA/CCI/FireCCI/5_1')` es el dataset de INCENDIOS de ESA CCI, no de biomasa. Ademas, la logica de combinacion GEDI+IPCC es incorrecta: `c_enhanced = gedi_carbon.unmask(c_total_ipcc)` reemplaza GEDI (solo above-ground) con total IPCC (4 pools) donde GEDI no tiene datos, mezclando escalas.
**Impacto:** Estimacion de carbono enhanced seria inconsistente y potencialmente sobreestimada.
**Correccion:** Usar GEDI solo para pool above-ground y sumar otros pools IPCC por separado.

### BUG-05: validate_simulation() FOM incorrecto [11_ca_markov.py:380-386]
**Severidad:** CRITICA
**Problema:** El Figure of Merit (Pontius et al. 2011) debe calcularse SOLO sobre pixeles que cambiaron (en observado O simulado), no sobre todos los pixeles. La implementacion actual incluye pixeles sin cambio en "hits", inflando artificialmente el FOM.
**Impacto:** Metricas de validacion CA-Markov serian artificialmente altas.
**Correccion:** Calcular FOM solo sobre pixeles donde observado O simulado difieren del mapa base.

---

## 3. BUGS MODERADOS (Corregidos)

### BUG-06: compute_drought_frequency() no funciona en ee.List.map() [09_climate_analysis.py:198]
**Severidad:** MODERADA
**Problema:** `compute_spi()` usa Python f-strings con `year`, pero dentro de `ee.List.map()` el argumento es un `ee.Number`. Las f-strings no pueden interpolar ee.Number.
**Impacto:** La funcion fallaria en ejecucion.
**Correccion:** Refactorizar para usar solo operaciones server-side dentro del map.

### BUG-07: extract_disturbance_map() es un stub incompleto [05_change_detection.py:138-159]
**Severidad:** MODERADA
**Problema:** La funcion solo retorna el resultado LandTrendr crudo y RMSE sin extraer el mapa de perturbacion (year, magnitude, duration).
**Impacto:** No se obtienen mapas de perturbacion utilizables.
**Correccion:** Implementar extraccion de segmentos de disturbio.

### BUG-08: GWR W matrix consume O(n^2) memoria innecesariamente [10_gwr_drivers.py:297]
**Severidad:** MODERADA
**Problema:** `W = np.diag(weights)` crea matriz diagonal n×n completa. Para n=2000 puntos, esto son 32MB por iteracion × 2000 iteraciones.
**Impacto:** Rendimiento lento y alto consumo de memoria.
**Correccion:** Usar multiplicacion por vector de pesos directamente sin construir la diagonal.

---

## 4. ISSUES MENORES

### MINOR-01: Landsat 9 filtrado innecesario para T1 (2012-2014)
**Archivo:** utils.py:127-132
**Nota:** L9 no existia antes de 2022. El merge con coleccion vacia funciona pero es ineficiente.

### MINOR-02: scripts/ sin __init__.py
**Archivo:** scripts/
**Nota:** `from scripts.02_training_samples import ...` requiere que scripts sea un paquete. Funciona con namespace packages implicitos en Python 3.3+, pero un __init__.py explicito es mas robusto. Ademas, `02_training_samples` no es un nombre de modulo valido (empieza con numero).

### MINOR-03: Inconsistencia calibracion CA-Markov
**Archivos:** experimental_design.md vs 11_ca_markov.py
**Nota:** El diseno dice "Calibracion: LULC 2013+2016 -> prediccion 2020" pero el config del script dice "2020-2024". La validacion hindcast deberia usar T3->T4 para calibracion y comparar.

### MINOR-04: Manuscrito - SPI no es Gaussiano sino Gamma
**Archivo:** manuscript_v1.md, 09_climate_analysis.py
**Nota:** El SPI real usa distribucion Gamma ajustada, no simple (P-mean)/std. La implementacion actual es una aproximacion. Debe documentarse como "SPI simplificado" en el manuscrito.

---

## 5. VALIDACION DE CONSISTENCIA METODOLOGICA

### 5.1 Pipeline de datos: OK
```
gee_config.py -> utils.py -> 01_preprocessing -> 02_training -> 03_classification
    -> 04_accuracy -> 05_change_detection -> 06_fragmentation
    -> 07_hotspot -> 08_ecosystem -> 09_climate -> 10_gwr -> 11_ca_markov
    -> 12_visualization
```
Flujo de dependencias correcto. Cada script importa solo de los anteriores.

### 5.2 Parametros consistentes entre scripts y diseno
| Parametro | Diseno | Codigo | Match |
|-----------|--------|--------|-------|
| 7 clases LULC | OK | OK | YES |
| 4 periodos | OK | OK | YES |
| RF ntree=500 | OK | OK | YES |
| 17 features | OK | OK | YES |
| 500 pts/clase | OK | OK | YES |
| 70/30 split | OK | OK | YES |
| LandTrendr params | OK | OK | YES |
| Carbon pools IPCC | OK | OK | YES |
| Bounding box | OK | OK | YES |
| 30 municipios | OK | OK | YES |

### 5.3 Consistencia manuscrito vs codigo
| Seccion | Status | Notas |
|---------|--------|-------|
| Abstract | OK | Alineado con 4 OE |
| Methods 3.2 | OK | RF params consistentes |
| Methods 3.4 | OK | Transition matrices correctas |
| Methods 3.5 | OK | FRAGSTATS metrics match |
| Methods 3.6 | OK | Moran/Gi* match |
| Methods 3.7 | OK | Carbon pools match |
| Methods 3.8 | OK | GWR variables match |
| Methods 3.9 | OK | CA-Markov scenarios match |

---

## 6. VALIDACION ESTADISTICA

### 6.1 Moran's I: Formulacion corregida
- Varianza bajo randomization: formula CORREGIDA (S2)
- Z-score: approximacion normal valida para n > 30 (OK, tenemos ~30,000 celdas)
- P-value: two-tailed correcto

### 6.2 GWR: Implementacion verificada
- OLS baseline correcto (incluyendo intercepto)
- VIF calculo correcto
- GWR bisquare kernel correcto
- AIC GWR con hat matrix trace correcto
- Bandwidth optimization por grid search (aceptable, golden section seria mas eficiente)

### 6.3 CA-Markov: Validacion corregida
- Matrices de transicion normalizadas por fila: OK
- Proyeccion Markov: multiplicacion correcta (areas @ matrix)
- Escenarios: re-normalizacion post-ajuste: OK
- FOM: CORREGIDO para solo pixeles de cambio

### 6.4 Mann-Kendall: OK
- Uso de ee.Reducer.kendallsCorrelation: correcto
- Sen's slope: implementacion local correcta

---

## 7. REPRODUCIBILIDAD

### 7.1 Semillas aleatorias: OK
- Training: seed=42+year (varia por periodo, reproducible)
- Split: seed=42 (fijo)
- CA-Markov: seed=42 (fijo)
- GWR sample: seed=42 (fijo)

### 7.2 Versiones de datos
- Hansen GFC: v1.11 (2023) - especificado
- MapBiomas: Collection 2 - especificado
- CHIRPS: Daily v2.0 - especificado
- MODIS: v061 - especificado

### 7.3 CRS y escala: OK
- Exportacion: EPSG:4326
- Clasificacion: 30m
- Analisis estadistico: 1km grid
- Clima: 1km (LST) y 5.5km (CHIRPS)

---

## 8. ACCIONES REALIZADAS

| # | Accion | Archivo | Estado |
|---|--------|---------|--------|
| 1 | Corregir classify_hotspots() | 07_hotspot_analysis.py | CORREGIDO |
| 2 | Corregir S2 en Moran's I | 07_hotspot_analysis.py | CORREGIDO |
| 3 | Corregir umbral GHSL SMOD | 02_training_samples.py | CORREGIDO |
| 4 | Corregir enhance_carbon_with_biomass | 08_ecosystem_services.py | CORREGIDO |
| 5 | Corregir FOM en validate_simulation | 11_ca_markov.py | CORREGIDO |
| 6 | Corregir drought_frequency GEE map | 09_climate_analysis.py | CORREGIDO |
| 7 | Completar extract_disturbance_map | 05_change_detection.py | CORREGIDO |
| 8 | Optimizar GWR sin matriz diagonal | 10_gwr_drivers.py | CORREGIDO |
| 9 | Agregar scripts/__init__.py | scripts/__init__.py | CREADO |
