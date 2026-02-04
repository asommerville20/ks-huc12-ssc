import geopandas as gpd
from pathlib import Path

from st0_config import data_dir, crs_wgs84

# spatially join SSC points to huc12 polygons and count membership
# 1) join stations to HUC12 polygons
## read in HUC12 and station gpkgs
stations_fp = data_dir / 'wqp_ssc_stations.gpkg'
stations_ks = gpd.read_file(stations_fp, layer="stations_ks").to_crs(crs_wgs84)

huc12_fp = data_dir / 'huc12_ks.gpkg'
huc12_ks = gpd.read_file(huc12_fp, layer="huc12_ks").to_crs(crs_wgs84)

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
huc12_ks_joined = huc12_ks.merge(sites_per_huc, on='huc12', how='left')
huc12_ks_joined['n_ssc_sites'] = huc12_ks_joined['n_ssc_sites'].fillna(0).astype(int)

print('HUC12 with >=1 SSC site:', 
      (huc12_ks_joined['n_ssc_sites'] > 0).sum(), 'of', len(huc12_ks))

# write to disk
huc12_ks_joined_fp = data_dir / 'huc12_ks_joined.gpkg'
huc12_ks_joined.to_file(huc12_ks_joined_fp, layer='huc12_ks_joined', driver='GPKG')
print(f'Saved: {huc12_ks_joined_fp}')
