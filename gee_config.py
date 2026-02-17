import ee
import os
from dotenv import load_dotenv

load_dotenv()

# Inicializar Earth Engine
try:
    ee.Initialize(project=os.getenv('GEE_PROJECT_ID'))
    print(f"GEE inicializado: {os.getenv('GEE_PROJECT_ID')}")
except Exception as e:
    print(f"Error GEE: {e}")
    raise

# ============================================================
# AREA DE ESTUDIO: Magdalena Medio (30 municipios)
# ============================================================

# Bounding box ajustado para la region del Magdalena Medio
# Cubre: Santander, Antioquia, Bolivar, Cesar (zona Magdalena Medio)
STUDY_AREA_BBOX = ee.Geometry.Rectangle([-75.0, 6.0, -73.5, 8.0])

# Municipios del Magdalena Medio por departamento
MUNICIPIOS = {
    'santander': [
        'Barrancabermeja', 'Puerto Wilches', 'Sabana de Torres',
        'Cimitarra', 'Puerto Parra', 'Rionegro',
        'San Vicente de Chucuri', 'El Carmen de Chucuri',
        'Betulia', 'Bolivar', 'Landazuri', 'Simacota',
        'Santa Helena del Opon', 'El Penon'
    ],
    'antioquia': [
        'Puerto Berrio', 'Yondo', 'Puerto Nare',
        'Puerto Triunfo', 'Caracoli', 'Maceo'
    ],
    'bolivar': [
        'San Pablo', 'Simiti', 'Santa Rosa del Sur',
        'Cantagallo', 'Morales', 'Arenal'
    ],
    'cesar': [
        'Aguachica', 'Gamarra', 'San Martin', 'San Alberto'
    ]
}

# Para cargar limites municipales reales desde DANE:
# admin2 = ee.FeatureCollection("FAO/GAUL/2015/level2")
# colombia_mun = admin2.filter(ee.Filter.eq('ADM0_NAME', 'Colombia'))
# Filtrar por nombres de municipios y hacer dissolve

# ============================================================
# PERIODOS DE ANALISIS
# ============================================================

PERIODS = {
    'pre_acuerdo': {
        'label': 'T1: Pre-acuerdo (conflicto activo)',
        'map_year': 2013,
        'start': '2012-01-01',
        'end': '2014-12-31',
        'context': 'Negociaciones La Habana; FARC controlan territorios'
    },
    'transicion': {
        'label': 'T2: Transicion (post-cese al fuego)',
        'map_year': 2016,
        'start': '2015-01-01',
        'end': '2017-06-30',
        'context': 'Firma acuerdo nov 2016; desmovilizacion; vacios de gobernanza'
    },
    'post_acuerdo_1': {
        'label': 'T3: Post-acuerdo temprano',
        'map_year': 2020,
        'start': '2019-01-01',
        'end': '2021-06-30',
        'context': 'Implementacion PDET; COVID-19; expansion ganadera'
    },
    'post_acuerdo_2': {
        'label': 'T4: Post-acuerdo reciente',
        'map_year': 2024,
        'start': '2023-01-01',
        'end': '2024-12-31',
        'context': 'Gobierno Petro; Paz Total; reduccion deforestacion nacional'
    }
}

# ============================================================
# COLECCIONES SATELITALES GEE
# ============================================================

COLLECTIONS = {
    'landsat8': 'LANDSAT/LC08/C02/T1_L2',
    'landsat9': 'LANDSAT/LC09/C02/T1_L2',
    'sentinel2': 'COPERNICUS/S2_SR_HARMONIZED',
    'modis_ndvi': 'MODIS/061/MOD13Q1',
    'modis_lst': 'MODIS/061/MOD11A2',
    'chirps': 'UCSB-CHG/CHIRPS/DAILY',
    'era5': 'ECMWF/ERA5_LAND/MONTHLY_AGGR',
    'srtm': 'USGS/SRTMGL1_003',
    'hansen': 'UMD/hansen/global_forest_change_2024_v1_12',
    'jrc_water': 'JRC/GSW1_4/GlobalSurfaceWater',
    'worldpop': 'WorldPop/GP/100m/pop',
    'ghsl': 'JRC/GHSL/P2023A/GHS_SMOD_V2-0',
    'soilgrids': 'projects/soilgrids-isric/assets/BDOD/BDOD_0-5cm_mean',
}

# ============================================================
# CLASES LULC (7 clases)
# ============================================================

LULC_CLASSES = {
    1: {'name': 'Bosque denso', 'color': '#006400', 'description': 'Cobertura arborea >60%'},
    2: {'name': 'Bosque secundario', 'color': '#32CD32', 'description': 'Cobertura arborea 30-60%, sucesion'},
    3: {'name': 'Pasturas', 'color': '#FFD700', 'description': 'Ganaderia, pastos naturales e introducidos'},
    4: {'name': 'Cultivos', 'color': '#FF8C00', 'description': 'Palma, arroz, coca, otros'},
    5: {'name': 'Agua', 'color': '#0000FF', 'description': 'Rios, cienagas, humedales'},
    6: {'name': 'Urbano', 'color': '#FF0000', 'description': 'Centros poblados, infraestructura'},
    7: {'name': 'Suelo desnudo', 'color': '#8B4513', 'description': 'Areas expuestas, mineria'},
}

# ============================================================
# PARAMETROS RANDOM FOREST
# ============================================================

RF_PARAMS = {
    'numberOfTrees': 500,
    'minLeafPopulation': 5,
    'bagFraction': 0.632,
    'seed': 42,
}

# ============================================================
# BANDAS E INDICES ESPECTRALES
# ============================================================

# Mapeo de bandas Landsat 8/9
LANDSAT_BANDS = {
    'blue': 'SR_B2',
    'green': 'SR_B3',
    'red': 'SR_B4',
    'nir': 'SR_B5',
    'swir1': 'SR_B6',
    'swir2': 'SR_B7',
}

# Mapeo de bandas Sentinel-2
SENTINEL_BANDS = {
    'blue': 'B2',
    'green': 'B3',
    'red': 'B4',
    'nir': 'B8',
    'swir1': 'B11',
    'swir2': 'B12',
}

# Indices espectrales a calcular
SPECTRAL_INDICES = ['NDVI', 'EVI', 'NDWI', 'NDBI', 'BSI', 'NBR', 'SAVI', 'MNDWI']

# ============================================================
# PARAMETROS LANDTRENDR
# ============================================================

LANDTRENDR_PARAMS = {
    'maxSegments': 6,
    'spikeThreshold': 0.9,
    'vertexCountOvershoot': 3,
    'preventOneYearRecovery': True,
    'recoveryThreshold': 0.25,
    'pvalThreshold': 0.05,
    'bestModelProportion': 0.75,
    'minObservationsNeeded': 6,
}

# ============================================================
# POOLS DE CARBONO (Mg C/ha) - Tier 2 Colombia (Alvarez et al. 2012)
# ============================================================

CARBON_POOLS = {
    1: {'c_above': 125, 'c_above_se': 15, 'c_below': 31, 'c_below_se': 8, 'c_soil': 57, 'c_soil_se': 12, 'c_dead': 18, 'c_dead_se': 5},  # Bosque denso (Alvarez 2012 + IFN Colombia)
    2: {'c_above': 55, 'c_above_se': 12, 'c_below': 14, 'c_below_se': 4, 'c_soil': 50, 'c_soil_se': 10, 'c_dead': 8, 'c_dead_se': 3},    # Bosque secundario (Alvarez 2012)
    3: {'c_above': 5, 'c_above_se': 2, 'c_below': 3, 'c_below_se': 1, 'c_soil': 35, 'c_soil_se': 8, 'c_dead': 0.5, 'c_dead_se': 0.3},    # Pasturas (FAO/IPCC)
    4: {'c_above': 8, 'c_above_se': 3, 'c_below': 2, 'c_below_se': 1, 'c_soil': 35, 'c_soil_se': 8, 'c_dead': 0.5, 'c_dead_se': 0.3},    # Cultivos
    5: {'c_above': 0, 'c_above_se': 0, 'c_below': 0, 'c_below_se': 0, 'c_soil': 0, 'c_soil_se': 0, 'c_dead': 0, 'c_dead_se': 0},          # Agua
    6: {'c_above': 2, 'c_above_se': 1, 'c_below': 0, 'c_below_se': 0, 'c_soil': 18, 'c_soil_se': 5, 'c_dead': 0, 'c_dead_se': 0},          # Urbano
    7: {'c_above': 0, 'c_above_se': 0, 'c_below': 0, 'c_below_se': 0, 'c_soil': 15, 'c_soil_se': 4, 'c_dead': 0, 'c_dead_se': 0},          # Suelo desnudo
}

print("Configuracion completa cargada")
print(f"  Area de estudio: Magdalena Medio ({len(sum(MUNICIPIOS.values(), []))} municipios)")
print(f"  Periodos: {len(PERIODS)}")
print(f"  Clases LULC: {len(LULC_CLASSES)}")
print(f"  Colecciones GEE: {len(COLLECTIONS)}")
