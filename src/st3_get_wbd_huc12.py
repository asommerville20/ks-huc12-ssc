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

from st0_config import data_dir, crs_wgs84, crs_int

def query_arcgis_features(layer_url,
                          where='1=1',
                          envelope=None,
                          out_fields='*',
                          chunk=200):
    '''
    Query ArcGIS REST MapServer layer and return all features as GeoDataFrame.
    Uses resultOffset pagination. Requests GeoJSON output.
    '''
    features = []
    offset = 0

    while True:
        params = {
            'f': 'geojson',
            'where': where,
            'outFields': out_fields,
            'returnGeometry': 'true',
            'outSR': crs_int,
            'resultRecordCount': chunk,
            'resultOffset': offset
        }

        if envelope is not None:
            geom = f"{envelope['xmin']},{envelope['ymin']},{envelope['xmax']},{envelope['ymax']}"
            params.update({
                'geometry': geom,
                'geometryType': 'esriGeometryEnvelope',
                'spatialRel': 'esriSpatialRelIntersects',
                'inSR': crs_int
            })

        r = requests.get(f'{layer_url}/query', params=params, timeout=120)
        r.raise_for_status()
        gj = r.json()

        batch = gj.get('features', [])
        if not batch:
            break

        # keep only features with geometry
        batch_valid = [f for f in batch if f.get('geometry') is not None]
        features.extend(batch_valid)

        # advance by server batch size, not valid count
        offset += len(batch)

        if not gj.get('exceededTransferLimit', False):
            break

    # if nothing returned, return empty GeoDataFrame with geometry + CRS
    if len(features) == 0:
        return gpd.GeoDataFrame(columns=['geometry'], 
                                geometry='geometry', 
                                crs=crs_wgs84)

    gdf = gpd.GeoDataFrame.from_features(features)

    # ensure geometry is set
    if 'geometry' in gdf.columns:
        gdf = gdf.set_geometry('geometry')
    else:
        gdf = gpd.GeoDataFrame(gdf, geometry='geometry')

    # ensure CRS is EPSG:4326
    if gdf.crs is None:
        gdf = gdf.set_crs(crs_wgs84)
    else:
        gdf = gdf.to_crs(crs_wgs84)

    return gdf




def get_wbd_huc12_ks(kansas_gdf):
    '''
    Pull WBD HU12 polygons that intersect Kanas, then clip to Kansas
    '''
    # WBD map server layer 6 is HU12 polygons
    huc12_layer = 'https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6'

    minx,miny,maxx,maxy = kansas_gdf.total_bounds
    env = {
        'xmin': float(minx),
        'ymin': float(miny),
        'xmax': float(maxx),
        'ymax': float(maxy)
    }

    huc12 = query_arcgis_features(
        layer_url = huc12_layer,
        where = '1=1',
        envelope = env,
        out_fields = 'huc12,name,states,areasqkm,tnmid'
    )

    # clip to kansas by intersection
    huc12_clip = gpd.overlay(huc12, kansas_gdf[['geometry']], how='intersection')
    huc12_clip = huc12_clip.to_crs(crs_wgs84)

    return huc12_clip

# attempt to run huc12_clip generation function
## read in boundary
ks_fp = f'{data_dir}/ks_boundary.gpkg'
kansas = gpd.read_file(ks_fp, layer="kansas").to_crs(crs_wgs84)

huc12_ks = get_wbd_huc12_ks(kansas)

## save to disk
huc12_fp = f'{data_dir}/huc12_ks.gpkg'
huc12_ks.to_file(huc12_fp, layer='huc12_ks', driver='gpkg')

print('HUC12 polygons returned:', len(huc12_ks))
print('HUC12 CRS:', huc12_ks.crs)
print('HUC12 columns:', list(huc12_ks.columns))

# map
## read in stations
stations_fp = f'{data_dir}/wqp_ssc_stations.gpkg'
stations_ks = gpd.read_file(stations_fp, layer="kansas").to_crs(crs_wgs84)

## plot
fig, ax = plt.subplots(figsize=(8,8))
kansas.boundary.plot(ax=ax, linewidth=1)
huc12_ks.boundary.plot(ax=ax, linewidth=0.4)
stations_ks.plot(ax=ax, markersize=5, color='black')

ax.set_title('Kansas: WBD HUC12 + SSC Stations')
plt.tight_layout()
plt.show()
#plt.savefig(out_dir / 'huc12_ssc_sites_KS.png', dpi=200)

# quick checks
## check for valid geometry
print('Valid geometry count:')
print(huc12_ks.geometry.is_valid.value_counts())
## check for unique geometry
print('Unique geometry check:')
print(huc12_ks.nunique())
print('Total geometry count:')
print(len(huc12_ks))