import geopandas as gpd
import requests
from pathlib import Path
import matplotlib.pyplot as plt
from io import BytesIO
import pandas as pd
import numpy as np

data_dir = Path('data')
out_dir = Path('outputs')

url = 'https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_state_20m.zip'
shp_zip = data_dir / 'cb_2023_us_state_20m.zip'

# download step
if not shp_zip.exists():
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    shp_zip.write_bytes(r.content)

# read in file and query for kansas 
states = gpd.read_file(f'zip://{shp_zip}')
kansas = states.loc[states['STUSPS'] == 'KS'].to_crs('EPSG:4326').copy()

print('ks crs', kansas.crs)
print('ks bounds', kansas.total_bounds)

#ax = kansas.boundary.plot(figsize=(6,6))
#ax.set_title('Kansas Boundary (EPSG:4326)')
#plt.tight_layout()
#plt.savefig(out_dir / 'kansas_boundary.png', dpi = 200)
# ----------------------------------------------------------------------------------- #



def get_wqp_ssc_stations_ks(
        pcode='80154', providers=None
):
    '''
    Download monitoring locations in Kansas that have SSC in WQP metadata.
    Returns a geodataframe (points).    
    '''
    base = 'https://www.waterqualitydata.us/data/Station/search'
    params = {
        'statecode': 'US:20',
        'pCode': str(pcode),
        'mimeType': 'geojson',
        'zip': 'no'
    }

    # for specifying providers
    if providers is not None:
        params['providers'] = str(providers)

    # build url that includes parameters
    import requests
    r = requests.Request('GET', base, params=params).prepare()
    url = r.url
    print('Request URL:\n', url, '\n')

    # download bytes directly
    resp = requests.get(url, timeout=120)
    if resp.status_code != 200:
        ## print small snippet of server message to see diagnose
        msg = resp.text[:1000] if resp.text else ""
        raise RuntimeError(f'WQP request failed: {resp.status_code}\n{msg}')

    gdf = gpd.read_file(BytesIO(resp.content))

    # standardize name fields
    def first_existing(cands):
        for c in cands:
            if c in gdf.columns:
                return c
        return None
    
    id_col = first_existing(["MonitoringLocationIdentifier", 
                             "monitoringLocationIdentifier", 
                             "stationId"])
    name_col = first_existing(["MonitoringLocationName", 
                               "monitoringLocationName", 
                               "stationName"])
    
    if id_col and id_col != 'site_id':
        gdf = gdf.rename(columns={id_col: 'site_id'})
    if name_col and name_col != 'site_name':
        gdf = gdf.rename(columns={name_col: 'site_name'})

    if 'site_id' not in gdf.columns:
        gdf['site_id'] = pd.Series(range(len(gdf))).astype(str)
    else:
        gdf['site_id'] = gdf['site_id'].astype(str)
    
    if 'site_name' not in gdf.columns:
        gdf['site_name'] = ''

    # ensure CRS is in WGS84
    if gdf.crs is None:
        gdf = gdf.set_crs('EPSG:4326')
    else:
        gdf = gdf.to_crs('EPSG:4326')

    return gdf[['site_id', 'site_name', 'geometry']].copy()

# run function
stations = get_wqp_ssc_stations_ks(pcode = '80154', providers=None)
print('Stations rows:', len(stations))
print('Stations CRS:', stations.crs)
print('Stations columns:', list(stations.columns))

# keep only points in Kansas polygon
## spatial join
stations_ks = gpd.sjoin(stations, kansas[['geometry']],
                        predicate='within',
                        how='inner').drop(columns=['index_right'])
print('Stations within Kansas polygon:', len(stations_ks))
# ----------------------------------------------------------------------------------- #




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
            'outSR': 4326,
            'resultRecordCount': chunk,
            'resultOffset': offset
        }

        if envelope is not None:
            geom = f"{envelope['xmin']},{envelope['ymin']},{envelope['xmax']},{envelope['ymax']}"
            params.update({
                'geometry': geom,
                'geometryType': 'esriGeometryEnvelope',
                'spatialRel': 'esriSpatialRelIntersects',
                'inSR': 4326
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
                                crs='EPSG:4326')

    gdf = gpd.GeoDataFrame.from_features(features)

    # ensure geometry is set
    if 'geometry' in gdf.columns:
        gdf = gdf.set_geometry('geometry')
    else:
        gdf = gpd.GeoDataFrame(gdf, geometry='geometry')

    # ensure CRS is EPSG:4326
    if gdf.crs is None:
        gdf = gdf.set_crs('EPSG:4326')
    else:
        gdf = gdf.to_crs('EPSG:4326')

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
    huc12_clip = huc12_clip.to_crs('EPSG:4326')

    return huc12_clip

# attempt to run huc12_clip generation function
huc12_ks = get_wbd_huc12_ks(kansas)

print('HUC12 polygons returned:', len(huc12_ks))
print('HUC12 CRS:', huc12_ks.crs)
print('HUC12 columns:', list(huc12_ks.columns))

# map
#fig, ax = plt.subplots(figsize=(8,8))
#kansas.boundary.plot(ax=ax, linewidth=1)
#huc12_ks.boundary.plot(ax=ax, linewidth=0.4)
#stations_ks.plot(ax=ax, markersize=5, color='black')

#ax.set_title('Kansas: WBD HUC12 + SSC Stations')
#plt.tight_layout()
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
# ----------------------------------------------------------------------------------- #




# spatially join SSC points to huc12 polygons and count membership
# 1) join stations to HUC12 polygons
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
# ----------------------------------------------------------------------------------- #




# compute geometry metrics
# 1) define analysis CRS for meters
analysis_crs = 'EPSG:26914' # NAD83 

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
print(huc12_ks[['area_km2', 'perimeter_km', 'compact']]).describe()