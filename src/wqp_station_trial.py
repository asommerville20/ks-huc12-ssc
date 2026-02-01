import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from io import BytesIO

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
## pull kansas boundary shp again
data_dir = Path('data')
shp_zip = data_dir / 'cb_2023_us_state_20m.zip'
states = gpd.read_file(f'zip://{shp_zip}')
kansas = states.loc[states['STUSPS'] == 'KS'].to_crs('EPSG:4326').copy()

## spatial join
stations_ks = gpd.sjoin(stations, kansas[['geometry']],
                        predicate='within',
                        how='inner').drop(columns=['index_right'])
print('Stations within Kansas polygon:', len(stations_ks))

# basic map
fig, ax = plt.subplots(figsize=(7,7))
kansas.boundary.plot(ax=ax, linewidth=1)
stations_ks.plot(ax=ax, markersize=8)
ax.set_title('SSC Monitoring Locations KS (WQP)')
plt.tight_layout()

out_dir = Path('outputs')
plt.savefig(out_dir / 'ssc_sites_KS.png', dpi=200)


