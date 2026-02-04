import geopandas as gpd
from pathlib import Path
import rasterio
from rasterstats import zonal_stats

from st0_config import data_dir, yr, name_resolution

# read in PRISM raster
prism_rast = next(data_dir.glob(f'*ppt*{name_resolution}*{yr}*.tif'), None)
if prism_rast is None:
    prism_rast = next(data_dir.glob(f'*ppt*{name_resolution}*{yr}*.bil'), None)
if prism_rast is None:
    raise FileNotFoundError('PRISM raster not found for given yr/resolution')

# read in polygons from st5 output (metrics)
huc_fp = data_dir / 'huc12_metrics.gpkg'
huc12_ks = gpd.read_file(huc_fp, layer='huc12_ks_metrics')

# open raster and inspect CRS
with rasterio.open(prism_rast) as src:
    rast_crs = src.crs
    print(f'PRISM CRS: {rast_crs}')

# reproject HUC12 polygons to PRISM CRS
huc12_rast_crs = huc12_ks.to_crs(rast_crs)

# zonal CRS
zs = zonal_stats(
    huc12_rast_crs.geometry,
    prism_rast,
    stats=['mean'],
    nodata=None
)

huc12_ks['ppt_mean'] = [z['mean'] for z in zs]
print(huc12_ks['ppt_mean'].describe())

# save to disk
out_fp = data_dir / 'huc12_ks_ppt.gpkg'
huc12_ks.to_file(out_fp, layer='huc12_ks_ppt', driver='GPKG')
print(f'Saved HUC12 + ppt_mean to : {out_fp}')