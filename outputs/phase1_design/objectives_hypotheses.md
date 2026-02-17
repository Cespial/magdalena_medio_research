# FASE 1.2: Objetivos e Hipotesis
## Paper: Cambio de Uso del Suelo y Servicios Ecosistemicos en el Magdalena Medio Post-Acuerdo de Paz

**Autor:** Cristian Espinal - Camara de Comercio de Medellin
**Fecha:** 2026-02-13

---

## 1. TITULO PROPUESTO DEL PAPER

**Ingles (titulo principal):**
> "Post-conflict land use transitions and ecosystem service loss in Colombia's Magdalena Medio: A multi-temporal remote sensing analysis (2012-2024)"

**Titulo alternativo (mas corto):**
> "Land use change and ecosystem service degradation in a post-conflict tropical landscape: Evidence from Colombia's Magdalena Medio"

---

## 2. OBJETIVO GENERAL

Analizar las dinamicas espacio-temporales de cambio de uso y cobertura del suelo (LULCC) en la region del Magdalena Medio colombiano durante el periodo 2012-2024, cuantificando su impacto sobre los servicios ecosistemicos y evaluando los patrones espaciales de cambio asociados a la transicion post-acuerdo de paz, con el fin de generar insumos para la planificacion territorial sostenible y la implementacion efectiva de los PDET.

**Alineacion con SDGs:**
- **SDG 15 (Life on Land):** Indicador 15.1.1 (proporcion de superficie forestal), 15.3.1 (proporcion de tierras degradadas)
- **SDG 13 (Climate Action):** Cuantificacion de perdida de almacenes de carbono por deforestacion
- **SDG 16 (Peace, Justice and Strong Institutions):** Analisis del nexo paz-medio ambiente en contexto post-conflicto
- **SDG 6 (Clean Water):** Evaluacion de impactos en regulacion hidrica

---

## 3. OBJETIVOS ESPECIFICOS

### OE1: Clasificar y mapear los cambios de cobertura y uso del suelo
**Enunciado:** Clasificar la cobertura y uso del suelo en el Magdalena Medio para cuatro periodos (2012, 2016, 2020, 2024) mediante algoritmos de Random Forest en Google Earth Engine, utilizando imagenes Landsat 8-9 y Sentinel-2, y generar matrices de transicion que caractericen las trayectorias de cambio pre y post-acuerdo de paz.

**Variables:**
- Dependiente: Tipo de cobertura (7 clases: bosque primario, bosque secundario, pasturas, cultivos, cuerpos de agua, areas urbanas, suelo desnudo)
- Independiente: Periodo temporal (pre-acuerdo vs. post-acuerdo)

**Metricas de exito:**
- Overall Accuracy >= 85%
- Kappa coefficient >= 0.80
- F1-score por clase >= 0.75
- Matrices de transicion completas para 3 intervalos

**Indicadores SDG:** 15.1.1, 15.3.1

---

### OE2: Identificar patrones espaciales y hotspots de cambio
**Enunciado:** Identificar y caracterizar los patrones espaciales de deforestacion y cambio de uso del suelo mediante analisis de autocorrelacion espacial (Moran's I), deteccion de hotspots (Getis-Ord Gi*), y metricas de fragmentacion del paisaje, comparando los periodos pre-acuerdo (2012-2016) y post-acuerdo (2017-2024).

**Variables:**
- Dependiente: Tasa de cambio de cobertura (ha/anio, %/anio)
- Independientes: Distancia a vias, distancia a centros poblados, elevacion, pendiente, presencia historica de grupos armados, densidad poblacional

**Metricas de exito:**
- Moran's I significativo (p < 0.05)
- Hotspots con confianza >= 95%
- Metricas FRAGSTATS: NP, PD, LPI, ED, COHESION, AI

**Indicadores SDG:** 15.1.1, 15.2.1

---

### OE3: Cuantificar la perdida de servicios ecosistemicos
**Enunciado:** Estimar y comparar la provision de servicios ecosistemicos clave (almacenamiento de carbono, regulacion hidrica, calidad de habitat) entre los periodos pre y post-acuerdo de paz, utilizando modelos InVEST integrados con datos de cobertura derivados de GEE.

**Variables:**
- Dependientes: Stock de carbono (Mg C/ha), rendimiento hidrico (mm/anio), indice de calidad de habitat (0-1)
- Independiente: Tipo de cobertura y periodo temporal

**Modelos a utilizar:**
| Servicio Ecosistemico | Modelo InVEST | Datos de entrada |
|----------------------|---------------|------------------|
| Almacenamiento de carbono | Carbon Storage and Sequestration | LULC maps + carbon pools por clase |
| Regulacion hidrica | Seasonal Water Yield | LULC + precipitacion CHIRPS + ET + suelos |
| Calidad de habitat | Habitat Quality | LULC + amenazas + sensibilidad por habitat |

**Metricas de exito:**
- Cambio en carbono total (Mg C) con intervalos de confianza
- Cambio en rendimiento hidrico (mm) por subcuenca
- Cambio en indice de habitat promedio por municipio

**Indicadores SDG:** 13.2.2, 15.1.1, 6.6.1

---

### OE4: Modelar escenarios futuros y analizar drivers de cambio
**Enunciado:** Analizar los factores socioambientales que impulsan el cambio de uso del suelo mediante regresion geograficamente ponderada (GWR) y proyectar escenarios de cambio LULC al 2030 y 2040 utilizando modelos CA-Markov, evaluando sus implicaciones para los servicios ecosistemicos y la planificacion territorial.

**Variables (modelo de drivers):**
| Variable | Tipo | Fuente |
|----------|------|--------|
| Tasa de deforestacion | Dependiente | OE1 |
| Distancia a carreteras | Independiente/continua | OSM / IGAC |
| Distancia a rios | Independiente/continua | GEE HydroSHEDS |
| Elevacion (DEM) | Independiente/continua | SRTM 30m |
| Pendiente | Independiente/continua | Derivada de SRTM |
| Precipitacion media | Independiente/continua | CHIRPS |
| Densidad poblacional | Independiente/continua | WorldPop |
| Distancia a centros urbanos | Independiente/continua | GEE settlements |
| Presencia historica FARC | Independiente/categorica | CERAC / datos conflicto |
| Aptitud agropecuaria | Independiente/categorica | UPRA |

**Metricas de exito:**
- GWR: R2 local > 0.60, significancia p < 0.05
- CA-Markov: Validacion con OA > 85% (simulacion 2020 vs. observado 2020)
- Escenarios: BAU (Business as Usual), Conservacion, Desarrollo planificado

**Indicadores SDG:** 15.3.1, 13.2.1

---

## 4. HIPOTESIS CIENTIFICAS

### H1: Hipotesis de aceleracion post-conflicto (OE1)
> **H1:** La tasa de conversion de bosque a pasturas y cultivos en el Magdalena Medio se incremento significativamente (>30%) en el periodo post-acuerdo (2017-2024) respecto al periodo pre-acuerdo (2012-2016), asociada a la expansion ganadera y la especulacion de tierras tras la retirada de las FARC.

**Hipotesis nula (H0):** No existe diferencia significativa en la tasa de conversion bosque-agropecuario entre los periodos pre y post-acuerdo.

**Prueba estadistica:** Test t pareado / Mann-Whitney U para tasas de cambio; Chi-cuadrado para matrices de transicion.

**Datos necesarios:** Mapas LULC 2012, 2016, 2020, 2024; matrices de transicion.

---

### H2: Hipotesis de clustering espacial (OE2)
> **H2:** Los cambios de uso del suelo post-acuerdo presentan un patron de clustering espacial significativo (Moran's I > 0, p < 0.01), con hotspots de deforestacion concentrados en areas de baja presencia institucional y alta accesibilidad vial, particularmente en los municipios de Yondo, Puerto Berrio y la zona sur del Magdalena Medio.

**Hipotesis nula (H0):** Los cambios de uso del suelo se distribuyen aleatoriamente en el espacio (Moran's I = 0).

**Prueba estadistica:** Moran's I global; Getis-Ord Gi* local; Prueba de significancia por permutaciones (999 iteraciones).

**Datos necesarios:** Mapas de cambio rasterizados; grilla de analisis; variables socioambientales.

---

### H3: Hipotesis de degradacion de servicios ecosistemicos (OE3)
> **H3:** La deforestacion y cambio de uso del suelo post-acuerdo han generado una perdida significativa (>15%) en el almacenamiento de carbono regional, una reduccion en la regulacion hidrica, y una disminucion en la calidad de habitat, afectando desproporcionadamente las areas de bosque ripario del rio Magdalena y sus tributarios.

**Hipotesis nula (H0):** No existe cambio significativo en los indicadores de servicios ecosistemicos entre los periodos pre y post-acuerdo.

**Prueba estadistica:** Test t / Wilcoxon para diferencias en indicadores ES por periodo; ANOVA espacial para comparacion por zonas.

**Datos necesarios:** Outputs InVEST para cada periodo; estadisticas zonales por municipio y subcuenca.

---

### H4: Hipotesis de heterogeneidad espacial de drivers (OE4)
> **H4:** Los factores que impulsan la deforestacion en el Magdalena Medio presentan variacion espacial significativa, con la accesibilidad vial y la densidad poblacional como predictores dominantes en la zona norte (mas urbanizada), y la presencia historica de actores armados y la aptitud ganadera como predictores dominantes en la zona sur (mas rural).

**Hipotesis nula (H0):** Los coeficientes de regresion de los drivers de deforestacion son espacialmente estacionarios (modelo global OLS suficiente).

**Prueba estadistica:** Comparacion AICc entre modelo OLS global y GWR; test de variabilidad de coeficientes; bandwidth optimization.

**Datos necesarios:** Variables socioambientales espacializadas; matrices de cambio; datos de conflicto armado.

---

## 5. JUSTIFICACION DEL VALOR CIENTIFICO Y SOCIAL

### 5.1 Valor Cientifico

**Novedad:**
1. Primer analisis LULCC multi-clase y multi-temporal del Magdalena Medio usando GEE (los estudios previos solo analizan deforestacion binaria con Hansen GFC)
2. Primera cuantificacion de perdida de servicios ecosistemicos vinculada al post-conflicto en esta region
3. Integracion metodologica innovadora: RF-GEE + LandTrendr + InVEST + GWR + CA-Markov en un solo marco analitico
4. Primera aplicacion de GWR para analisis de drivers de deforestacion en Magdalena Medio

**Reproducibilidad:**
- Codigo abierto en GEE (JavaScript) y Python
- Datos satelitales de acceso libre (Landsat, Sentinel, CHIRPS, SRTM)
- Protocolo metodologico detallado y replicable

**Contribucion teorica:**
- Aporta al debate sobre el "paradox of peace" (como la paz puede generar degradacion ambiental)
- Contribuye a la comprension del nexo conflicto-medio ambiente en paisajes tropicales
- Proporciona evidencia para el marco teorico de "land use transition" en contextos post-conflicto

### 5.2 Valor Social y Politico

**Relevancia para politica publica:**
1. Insumo para los PDET (Planes de Desarrollo con Enfoque Territorial) del Magdalena Medio
2. Linea base para monitoreo de compromisos ambientales del acuerdo de paz
3. Identificacion de areas prioritarias para conservacion y restauracion
4. Cuantificacion de perdidas de servicios ecosistemicos para valoracion economica

**Relevancia para la Camara de Comercio de Medellin:**
- Datos para orientar inversion sostenible en la region
- Informacion para proyectos de desarrollo territorial
- Evidencia para iniciativas de economia verde y mercados de carbono

**Beneficiarios:**
- Agencia de Renovacion del Territorio (ART)
- IDEAM (monitoreo forestal)
- Corporaciones Autonomas Regionales (CARs)
- Comunidades rurales del Magdalena Medio
- Sector privado con inversion en la region
- Academia y comunidad cientifica internacional

---

## 6. MARCO CONCEPTUAL

```
                    ACUERDO DE PAZ 2016
                          |
              +-----------+-----------+
              |                       |
    Retirada FARC              Gobernanza debil
              |                       |
    +----v----+----+         +---v---+---+
    | Acceso a     |         | Vacios de |
    | territorios  |         | control   |
    +----+----+----+         +---+---+---+
         |    |                  |   |
         v    v                  v   v
    Ganaderia  Coca     Especulacion  Mineria
    extensiva  cultivos de tierras    ilegal
         |    |              |        |
         +----+----+---------+--------+
                   |
            CAMBIO LULC (OE1)
                   |
        +----------+----------+
        |          |          |
    Deforestacion  Fragmentacion  Degradacion
        |          |              |
        v          v              v
    Perdida    Perdida de     Perdida de
    carbono    conectividad   regulacion
    (OE3)      (OE2)         hidrica (OE3)
        |          |              |
        +----------+----------+---+
                   |
        SERVICIOS ECOSISTEMICOS
                   |
        +----------+----------+
        |          |          |
    Comunidades  Biodiversidad  Clima
        |          |              |
        v          v              v
    Bienestar   Extincion     Cambio
    humano      local         climatico
                   |
            ESCENARIOS FUTUROS (OE4)
            BAU / Conservacion / PDET
```

---

## 7. PREGUNTAS DE INVESTIGACION

Derivadas de los objetivos e hipotesis:

**P1 (OE1):** Cuales son las principales transiciones de uso del suelo en el Magdalena Medio entre 2012 y 2024, y como difieren las tasas de cambio entre los periodos pre y post-acuerdo de paz?

**P2 (OE2):** Donde se concentran espacialmente los cambios de uso del suelo y que patrones de fragmentacion del paisaje han emergido en el periodo post-acuerdo?

**P3 (OE3):** Cual es la magnitud de la perdida de servicios ecosistemicos (carbono, agua, habitat) atribuible a los cambios de uso del suelo post-acuerdo?

**P4 (OE4):** Que factores socioambientales explican la variacion espacial de la deforestacion y cuales son los escenarios probables de cambio al 2030 y 2040?

---

## 8. KEYWORDS PROPUESTOS

**Ingles:** Land use/land cover change; Post-conflict deforestation; Ecosystem services; Google Earth Engine; Magdalena Medio; Peace agreement

**Espanol:** Cambio de uso del suelo; Deforestacion post-conflicto; Servicios ecosistemicos; Google Earth Engine; Magdalena Medio; Acuerdo de paz
