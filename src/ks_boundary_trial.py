import geopandas as gpd
import requests
from pathlib import Path
import matplotlib.pyplot as plt

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

ax = kansas.boundary.plot(figsize=(6,6))
ax.set_title('Kansas Boundary (EPSG:4326)')
plt.tight_layout()
plt.savefig(out_dir / 'kansas_boundary.png', dpi = 200)
