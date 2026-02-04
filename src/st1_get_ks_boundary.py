import geopandas as gpd
import requests
from pathlib import Path
import matplotlib.pyplot as plt

from st0_config import data_dir, crs_wgs84

# set URL and zip destination
url = 'https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_state_20m.zip'
shp_zip = data_dir / 'cb_2023_us_state_20m.zip'
ks_fp = data_dir / 'kansas_boundary.gpkg'

# download step
if not shp_zip.exists():
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    shp_zip.write_bytes(r.content)

# read in file and query for kansas 
states = gpd.read_file(f'zip://{shp_zip}')
kansas = states.loc[states['STUSPS'] == 'KS'].to_crs(crs_wgs84).copy()

print('ks crs', kansas.crs)
print('ks bounds', kansas.total_bounds)

ax = kansas.boundary.plot(figsize=(6,6))
ax.set_title('Kansas Boundary (EPSG:4326)')
plt.tight_layout()
plt.show()

# save to disc
kansas.to_file(
    ks_fp,
    layer='kansas',
    driver='GPKG'
)
print(f'Saved Kansas boundary to: {ks_fp}')