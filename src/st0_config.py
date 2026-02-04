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

class Config:
    data_dir: Path = Path('data')
    out_dir: Path = Path('outputs')

    # CRS choices:
    crs_wgs84: str = 'EPSG:4326'
    crs_int = int = 4326
    crs_projected: str = 'EPSG:26914' # will be analysis CRS

    # WQP parameter code:
    pcode = '80154'

cfg = Config()
cfg.data_dr.mkdir(exist_ok=True)
cfg.out_dir.mkdir(exist_ok=True)

