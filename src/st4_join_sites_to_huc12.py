import geopandas as gpd
import requests
from pathlib import Path
import matplotlib.pyplot as plt
from io import BytesIO
import pandas as pd
import numpy as np
import zipfile
import rasterio
from rasterstats import zonal_stats

from st0_config import data_dir, crs_wgs84, crs_projected

# spatially join SSC points to huc12 polygons and count membership
# 1) join stations to HUC12 polygons
## read in HUC12 and station gpkgs
stations_fp = f'{data_dir}/wqp_ssc_stations.gpkg'
stations_ks = gpd.read_file(stations_fp, layer="kansas").to_crs(crs_wgs84)

huc12_fp = f'{data_dir}/huc12_ks.gpkg'
huc12_ks = gpd.read_file(stations_fp, layer="kansas").to_crs(crs_wgs84)

## spatial join
stations_huc12 = gpd.sjoin(
    stations_ks,
    huc12_ks[['huc12','geometry']],
    predicate='within',
    how='left'
).drop(columns=['index_right'], errors='ignore')

print('Stations with HUC12 assigned:', 
      stations_huc12['huc12'].notna().sum(), 'of', len(stations_huc12))
print('Stations missing HUC12:', stations_huc12['huc12'].isna().sum())

# 2) count SSC sites per HUC12
sites_per_huc = (
    stations_huc12.groupby('huc12')
    .size()
    .rename('n_ssc_sites')
    .reset_index()
)

# 3) merge counts back to polygons
huc12_ks = huc12_ks.merge(sites_per_huc, on='huc12', how='left')
huc12_ks['n_ssc_sites'] = huc12_ks['n_ssc_sites'].fillna(0).astype(int)

print('HUC12 with >=1 SSC site:', 
      (huc12_ks['n_ssc_sites'] > 0).sum(), 'of', len(huc12_ks))




# compute geometry metrics
# 1) define analysis CRS for meters
analysis_crs = crs_projected # NAD83 

# 2) reproject polygons to projected CRS for area/length
huc_proj = huc12_ks.to_crs(analysis_crs).copy()

# 3) geometry metrics
huc_proj['area_km2'] = huc_proj.geometry.area / 1e6
huc_proj['perimeter_km'] = huc_proj.geometry.length / 1000.0

# 4) compactness (1 is circle, smaller is elongated/irregular)
huc_proj['compact'] = (4*np.pi*huc_proj.geometry.area) / (huc_proj.geometry.length**2)

# 5) bring metrics back to WGS84 layer (keep geo in EPSG:4326 for mapping)
huc12_ks = huc12_ks.merge(
    huc_proj[['huc12', 'area_km2', 'perimeter_km', 'compact']],
    on='huc12',
    how='left'
)

# quick check
print(huc12_ks[['area_km2', 'perimeter_km', 'compact']].describe())




# make projection-correctness stats
# 1) pick single huc12 with station
one_huc = huc12_ks.loc[huc12_ks['n_ssc_sites'] > 0].iloc[0:1].copy()

# 2) WRONG: area in EPSG:4326 (degrees units, not meaningful)
one_huc_wgs = one_huc.to_crs(crs_wgs84)
wrong_area = one_huc_wgs.geometry.area.iloc[0]

# 3) CORRECT: area in projected CRS (m^2 to km^2)
one_huc_proj = one_huc.to_crs(analysis_crs)
correct_area = (one_huc_proj.geometry.area.iloc[0]) / 1e6

print(f'Example HUC12: {one_huc['huc12'].iloc[0]}')
print(f'Area computed in {crs_wgs84}: {wrong_area} <- (degrees^2, meaningless)')
print(f'Area computed in {analysis_crs}: {correct_area} <- (km^2)')