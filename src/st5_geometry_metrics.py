import geopandas as gpd
from pathlib import Path
import numpy as np

from st0_config import data_dir, crs_wgs84, crs_projected

# compute geometry metrics
# 1) read in joined huc12 gpkg
in_fp = data_dir / 'huc12_ks_joined.gpkg'
huc12_joined = gpd.read_file(in_fp, layer='huc12_ks_joined').to_crs(crs_wgs84)

# 2) reproject polygons to projected CRS for area/length
huc_proj = huc12_joined.to_crs(crs_projected).copy()

# 2) geometry metrics
huc_proj['area_km2'] = huc_proj.geometry.area / 1e6
huc_proj['perimeter_km'] = huc_proj.geometry.length / 1000.0

# 3) compactness (1 is circle, smaller is elongated/irregular)
huc_proj['compact'] = (4*np.pi*huc_proj.geometry.area) / (huc_proj.geometry.length**2)

# 4) bring metrics back to WGS84 layer (keep geo in EPSG:4326 for mapping)
huc12_metrics = huc12_joined.merge(
    huc_proj[['huc12', 'area_km2', 'perimeter_km', 'compact']],
    on='huc12',
    how='left'
)

# quick check
print(huc12_metrics[['area_km2', 'perimeter_km', 'compact']].describe())




# make projection-correctness stats
# 1) pick single huc12 with station
one_huc = huc12_metrics.loc[huc12_metrics['n_ssc_sites'] > 0].iloc[0:1].copy()

# 2) WRONG: area in EPSG:4326 (degrees units, not meaningful)
one_huc_wgs = one_huc.to_crs(crs_wgs84)
wrong_area = one_huc_wgs.geometry.area.iloc[0]

# 3) CORRECT: area in projected CRS (m^2 to km^2)
one_huc_proj = one_huc.to_crs(crs_projected)
correct_area = (one_huc_proj.geometry.area.iloc[0]) / 1e6

# 4) compare against already computed metric
metric_area_km2 = float(one_huc['area_km2'].iloc[0])

print(f'Example HUC12: {one_huc['huc12'].iloc[0]}')
print(f'Area computed in {crs_wgs84}: {wrong_area} <- (degrees², meaningless)')
print(f'Area computed in {crs_projected}: {correct_area} <- (km², correct)')
print(f'Area from stored metric: {metric_area_km2} <- (km²)')

# write to disk
out_fp = data_dir / 'huc12_metrics.gpkg'
huc12_metrics.to_file(out_fp, layer='huc12_ks_metrics', driver='GPKG')
print(f'Saved: {out_fp}')