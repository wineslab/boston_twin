import urllib.request
import zipfile
from pathlib import Path
import geopandas as gdp
import requests


class AppURLopener(urllib.request.FancyURLopener):
    version = "Mozilla/5.0"


BASE_MODEL_URL = "https://maps.bostonplans.org/3d/Bos3d_BldgModels_20230927_OBJ"
BASE_GROUND_URL = (
    "https://maps.bostonplans.org/3d/Bos3d_Terrain_2011_OBJ"
)

letters = ["A", "B", "C", "D", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
nums = range(1,13)

save_dir = Path("data")
if not save_dir.is_dir():
    save_dir.mkdir(parents=True, exist_ok=True)

# %% Project file
# source: https://www.cityschema.org/tile_scheme/index.htm
proj_url = "https://cityschema.github.io/repository-catalog/Bos3d_CityWide_Data/Bos3d_TIleGrid/Metro_Boston_3D_CRS.zip"

zip_proj_file_path = save_dir.joinpath("Metro_Boston_3D_CRS.zip")
proj_file_path = save_dir.joinpath("Metro_Boston_3D_CRS.prj")

if not proj_file_path.is_file():
    # req = urlrequest.Request(link, headers={'User-Agent': 'XYZ/3.0'})
    r = requests.get(proj_url, stream=True, headers={"User-Agent": "'XYZ/3.0'"})
    if not r.status_code==404:
        print("Downloading the 3D projection file...")

        with open(zip_proj_file_path, "wb") as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        print("Extracting..")
        with zipfile.ZipFile(zip_proj_file_path, "r") as zip_ref:
            zip_ref.extractall(save_dir)
        zip_proj_file_path.unlink()

# %% Model files
downloaded = []
for l in letters:
    for n in nums:
        filename = f"BOS_{l}_{n}_BldgModels_OBJ"
        downloaded.append(filename)
        out_tile_dir = save_dir.joinpath(filename)
        zip_file_path = save_dir.joinpath(filename + ".zip")
        if Path(out_tile_dir).is_dir():
            print(f"{out_tile_dir} already downloaded. Skipping.")
        else:
            if ~zip_file_path.is_file():
                try:
                    url = BASE_MODEL_URL + "/" + filename + ".zip"

                    # req = urlrequest.Request(link, headers={'User-Agent': 'XYZ/3.0'})
                    r = requests.get(
                        url, stream=True, headers={"User-Agent": "XYZ/3.0"}
                    )
                    if r.status_code==404:
                        continue
                    print("Downloading " + filename + "...")
                    # print(r.status_code)
                    # z = zipfile.ZipFile(io.BytesIO(r.content))
                    with open(zip_file_path, "wb") as fd:
                        for chunk in r.iter_content(chunk_size=128):
                            fd.write(chunk)

                    print("Extracting..")
                    with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
                        zip_ref.extractall(out_tile_dir)
                    zip_file_path.unlink()
                except FileNotFoundError as e:
                    print(url)
                    continue

        ground_name = f"BOS_{l}_{n}_TerrainMesh_2011_OBJ"
        out_ground_path = save_dir.joinpath(ground_name)
        zip_ground_path = save_dir.joinpath(ground_name + ".zip")
        if out_ground_path.is_dir():
            print(f"{ground_name} already downloaded. Skipping.")
            continue

        if ~out_ground_path.is_dir():
            if ~zip_ground_path.is_file():
                try:
                    url = BASE_GROUND_URL + "/" + ground_name + ".zip"
                    print(url)
                    r = requests.get(
                        url, headers={"User-Agent": "XYZ/3.0"}, stream=True
                    )
                    if r.status_code==404:
                        continue
                    print("Downloading " + ground_name + "...")
                    # z = zipfile.ZipFile(io.BytesIO(r.content))
                    with open(zip_ground_path, "wb") as fd:
                        for chunk in r.iter_content(chunk_size=128):
                            fd.write(chunk)

                    print("Extracting..")
                    with zipfile.ZipFile(zip_ground_path, "r") as zip_ref:
                        zip_ref.extractall(out_ground_path)
                    zip_ground_path.unlink()
                except FileNotFoundError as e:
                    print(url)
                    continue

# %% extract single OBJ models
extract_objs = False
if extract_objs:
    for f in downloaded:
        filepath = Path(f)
        for fzip in filepath.iterdir():
            print(f"Extracting individual models from {fzip}..")
            if "_OBJ" in fzip.stem and fzip.suffix == ".zip":
                with zipfile.ZipFile(str(fzip.resolve()), "r") as zip_ref:
                    zip_ref.extractall(str(fzip.parent.resolve()))
                    zip_ref.unlink()

# # %% extract and convert geodataframe
# extract_geodf = True
# if extract_geodf:
#     read_columns = [
#         "Centr_Lat",
#         "Centr_Lon",
#         "Gnd_El_Ft",
#         "Z_Max_Ft",
#         "Z_MIn_Ft",
#         "Height_Ft",
#         "geometry",
#         "StructType",
#     ]
#     write_columns = ["Z_m", "geometry", "StructType"]
#     new_rows = []
#     for f in downloaded:
#         filepath = Path(f, "ModelCatalog", "ModelCatalog.geojson")
#         gdf = gdp.read_file(filepath)
#         gdf = gdf[read_columns]
#         new_row = gdp.GeoDataFrame(columns=write_columns, geometry="geometry")
#         new_row.loc[:, "Z_m"] = gdf.apply(
#             lambda x: str2float(x["Z_Max_Ft"]) * 0.3048, axis=1
#         )
#         new_row.loc[:, "geometry"] = gdf["geometry"]
#         new_row.loc[:, "StructType"] = gdf["StructType"]
#         new_rows.append(new_row)
#     new_gdf = gdp.GeoDataFrame(pd.concat(new_rows), geometry="geometry")
#     new_gdf.to_file("complete_models.geojson")
