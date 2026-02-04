import geopandas as gpd
from pathlib import Path
from io import BytesIO
import pandas as pd

from st0_config import data_dir, crs_wgs84, pcode

def get_wqp_ssc_stations_ks(
        pcode=pcode, providers=None
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
        gdf = gdf.set_crs(crs_wgs84)
    else:
        gdf = gdf.to_crs(crs_wgs84)

    return gdf[['site_id', 'site_name', 'geometry']].copy()

# run function
stations = get_wqp_ssc_stations_ks(pcode = pcode, providers=None)
print('Stations rows:', len(stations))
print('Stations CRS:', stations.crs)
print('Stations columns:', list(stations.columns))

# keep only points in Kansas polygon
## read in ks boundary gpkg
ks_fp = data_dir / 'kansas_boundary.gpkg'
kansas = gpd.read_file(ks_fp, layer='kansas').to_crs(crs_wgs84)

## spatial join
stations_ks = gpd.sjoin(stations, kansas[['geometry']],
                        predicate='within',
                        how='inner').drop(columns=['index_right'])
print('Stations within Kansas polygon:', len(stations_ks))

# save sites to disk
stations_fp = data_dir / 'wqp_ssc_stations.gpkg'
stations_ks.to_file(stations_fp, layer="stations_ks", driver="GPKG")
print(f'Saved: {stations_fp}')
