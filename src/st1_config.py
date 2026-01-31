from pathlib import Path

class Config:
    data_dir: Path = Path('data')
    out_dir: Path = Path('outputs')

    # CRS choices:
    crs_wgs84: str = 'EPSG:4326'
    crs_projected: str = 'EPSG:26914'

cfg = Config()
cfg.data_dr.mkdir(exist_ok=True)
cfg.out_dir.mkdir(exist_ok=True)