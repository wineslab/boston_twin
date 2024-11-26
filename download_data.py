import zipfile
from pathlib import Path

import requests

NU_URL = "https://repository.library.northeastern.edu/downloads/neu:ms36tq790?datastream_id=content"

download_dir = Path(".")
try:
    print("Downloading the dataset from the Northeastern repository..")
    zip_dataset_path = download_dir.joinpath("BostonTwinDataset.zip")
    r = requests.get(NU_URL, stream=True, headers={"User-Agent": "'XYZ/3.0'"})
    if not r.status_code == 404:
        with open(zip_dataset_path, "wb") as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        print("Extracting..")
        with zipfile.ZipFile(zip_dataset_path, "r") as zip_ref:
            zip_ref.extractall(download_dir)
        zip_dataset_path.unlink()
except FileNotFoundError as e:
    print(
        f"Can't download from the Northeastern repository. Try the BPDA website. ({e})"
    )