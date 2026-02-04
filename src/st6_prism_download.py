import requests
from pathlib import Path
import zipfile

from st0_config import data_dir, yr

# 1) function to download PRISM zip
def download_prism_ppt(year, res='4km', region='us'):
    """
    Downloads PRISM annual precipitation raster via ZIP file

    params:
        year (int, optional): Defaults to 2023.
        res (str, optional): Defaults to '4km'.
    """
    url = f'https://services.nacse.org/prism/data/public/{res}/ppt/{year}'
    url = f'https://services.nacse.org/prism/data/get/{region}/{res}/ppt/{year}'
    out_zip = data_dir / f'ppt_{res}_{year}_prism.zip'

    if not out_zip.exists():
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        out_zip.write_bytes(r.content)

    return out_zip

# 2) function to extract PRISM raster from zip
def extract_raster_prism(zip_path):
    '''
    Extract raster file (.bil or .tif) from PRISM zip.
    :param zip_path: zip file path
    '''
    with zipfile.ZipFile(zip_path, 'r') as z:
        names = z.namelist()
        rast_candidates = [n for n in names if n.lower().endswith(('.bil','.tif'))]

        if not rast_candidates:
            raise RuntimeError(f'No raster found in {zip_path}')
        
        rast_name = rast_candidates[0]
        out_path = data_dir / Path(rast_name).name
        
        if not out_path.exists():
            z.extract(rast_name, path='data')
            extracted = data_dir / rast_name
            extracted = extracted.resolve()

            if extracted != out_path.resolve():
                out_path.parent.mkdir(parents=True,exist_ok=True)
                extracted.rename(out_path)

    return out_path

# 3) run 
prism_zip = download_prism_ppt(year=yr, res='4km')
prism_rast = extract_raster_prism(prism_zip)

print(f'Data saved: {prism_rast}')