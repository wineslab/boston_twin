import mitsuba as mi
from pathlib import Path
from src.utils.BostonModelDownloader import BostonModel
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

# import osmnx as ox

model_rootdir = Path("data/data")
antenna_rootdir = Path("data/")
output_rootdir = Path("bostontwin")
output_3d_dir = output_rootdir.joinpath("boston3d")
output_antenna_dir = output_rootdir.joinpath("boston_antennas")

print(model_rootdir.resolve())
model = BostonModel(model_rootdir, output_3d_dir)

print(model.boston_center_m)
model.get_dataset()

# url = BASE_GROUND_URL + "/" + ground_name + ".zip"
# print(url)
# r = requests.get(
#     url, headers={"User-Agent": "XYZ/3.0"}, stream=True
# )
# if r.status_code==404:
#     continue
# print("Downloading " + ground_name + "...")
# # z = zipfile.ZipFile(io.BytesIO(r.content))
# with open(zip_ground_path, "wb") as fd:
#     for chunk in r.iter_content(chunk_size=128):
#         fd.write(chunk)

das_approv = gpd.read_file(
    antenna_rootdir.joinpath("DAS_Small_Cell_Approved_Locations.geojson")
)
das_approv.rename(
    {
        "attachment": "Attachment_Or_Replacement",
        "existing_p": "Original_Pole_Type",
        "installati": "Install_Date",
        "intended_c": "Intended_Commercial_Use",
        "latitude": "Lat",
        "longitude": "Long",
        "pole_id": "Pole_Identifying_Number",
        "replacemen": "New_Pole_Type",
        "spectrum": "Spectrum",
        "street_add": "Address",
        "vendor": "Vendor",
    },
    axis=1,
    inplace=True,
)
das_approv.drop(["FID", "fid_", "activation", "other_infr"], axis=1, inplace=True)
das_approv["Install_Date"] = das_approv["Install_Date"].apply(lambda x: pd.Timestamp(x))
print(das_approv.columns)
das_approv.reindex(
    das_approv.columns.tolist()
    + ["Neighborhood", "City_Reference", "Requester_Email_address"],
    axis=1,
)
das_approv = das_approv.sort_index(axis=1)

wireless_approv = gpd.read_file(
    antenna_rootdir.joinpath("Wireless_Antenna_Installation_Requests_Approved.geojson")
)
wireless_approv.rename({"Planned_Install_Date": "Install_Date"}, axis=1, inplace=True)
wireless_approv["Install_Date"] = wireless_approv["Install_Date"].apply(
    lambda x: pd.Timestamp(x)
)
wireless_approv = wireless_approv.sort_index(axis=1)
wireless_approv = wireless_approv.sort_values("Install_Date", axis=0)

print(
    f"{wireless_approv.shape[0]-wireless_approv['Pole_Identifying_Number'].unique().shape[0]} duplicate IDs"
)
wireless_approv.drop_duplicates(
    "Pole_Identifying_Number", keep="last", ignore_index=True, inplace=True
)
print(
    f"{wireless_approv.shape[0]-wireless_approv['Pole_Identifying_Number'].unique().shape[0]} duplicate IDs"
)
# print(wireless_approv.duplicate())

print(das_approv.columns)
print(wireless_approv.columns)
# results = pd.concat([das_approv, wireless_approv])
# results = pd.merge(
#     left=wireless_approv,
#     right=das_approv,
# )

new_antennas_ids = wireless_approv["Pole_Identifying_Number"].unique()
drop_row_ids = []
for row_idx, pole_id in das_approv["Pole_Identifying_Number"].items():
    if pole_id in new_antennas_ids:
        drop_row_ids.append(row_idx)
print(f"Dropping {len(drop_row_ids)} rows")
das_approv.drop(drop_row_ids, inplace=True)

results = pd.concat(
    [das_approv, wireless_approv],
    ignore_index=True,
    # verify_integrity=True,
    axis=0,
)
print(
    f"{das_approv.shape[0]} antennas before 2017.\n{wireless_approv.shape[0]} antennas after 2017.\n{results.shape[0]} current antennas."
)
if not output_antenna_dir.is_dir():
    output_antenna_dir.mkdir(parents=True,exist_ok=True)
results.to_file(output_antenna_dir.joinpath("antennas.geojson"), driver="GeoJSON")

# roads = ox.graph_from_bbox(
#     west=west, south=south, east=east, north=north, network_type="all"
# )
# roads = ox.projection.project_graph(roads, to_crs="epsg:6348")
