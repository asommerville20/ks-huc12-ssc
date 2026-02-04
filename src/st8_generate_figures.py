import geopandas as gpd
from pathlib import Path
import matplotlib.pyplot as plt

from st0_config import data_dir, crs_wgs84

# input data from st7
in_fp = data_dir / 'huc12_ks_ppt.gpkg'
huc12_ks = gpd.read_file(in_fp, layer='huc12_ks_ppt').to_crs(crs_wgs84)

# 1) scatter: area vs ppt_mean
fig, ax = plt.subplots(figsize=(4,6))

ax.scatter(
    huc12_ks['area_km2'],
    huc12_ks['ppt_mean'],
    s=10,
    alpha=0.5
)
ax.set_xlabel('HUC12 Watershed Area (km²)')
ax.set_ylabel('PRISM PPT Mean (mm)')
ax.set_title('PRISM Precipitation vs Watershed Area (KS)')
plt.tight_layout()
plt.show()

# 2) map: precipitation gradient (choropleth)
fig, ax = plt.subplots(figsize=(6,6))

## polygons colored by ppt_mean
huc12_ks.plot(
    column='ppt_mean',
    ax=ax,
    legend=True
)

ax.set_title('PRISM Mean Annual Precipitation (mm) by HUC12')
plt.tight_layout()
plt.show()