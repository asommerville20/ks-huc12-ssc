from pathlib import Path

class Config:
    data_dir: Path = Path('data')
    out_dir: Path = Path('outputs')

    # CRS choices:
    crs_wgs84: str = 'EPSG:4326'
    crs_int = int = 4326
    crs_projected: str = 'EPSG:26914' # will be analysis CRS

    # WQP parameter code:
    pcode = '80154'

    # year:
    yr = 2023

    # resolution:
    resolution = '4km' # for quicker processing
    name_resolution = '25m'

cfg = Config()
cfg.data_dir.mkdir(exist_ok=True)
cfg.out_dir.mkdir(exist_ok=True)

# ---- module-level exports ----
data_dir = cfg.data_dir
out_dir = cfg.out_dir

crs_wgs84 = cfg.crs_wgs84
crs_int = cfg.crs_int
crs_projected = cfg.crs_projected

pcode = cfg.pcode
yr = cfg.yr
resolution = cfg.resolution
name_resolution = cfg.name_resolution
