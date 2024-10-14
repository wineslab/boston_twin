import json
import time
import zipfile
from pathlib import Path
from typing import Union

import geopandas as gpd
import mitsuba as mi
import numpy as np
import pandas as pd
import pyproj
import requests
import shapely as shp
from pyproj.transformer import Transformer
from shapely.geometry import box

from bostontwin.utils.geo_utils import check_area_of_use, gdf2crs, get_crs
from bostontwin.utils.obj_utils import dir_obj2ply
from bostontwin.utils.utils import generate_mi_xml, print_eta, truncate_utf8_chars

mi.set_variant("scalar_rgb")

PUBLIC_URL = "https://www.bostonplans.org/3d-data-maps/3d-smart-model/3d-data-download"
BASE_MODEL_URL = "https://maps.bostonplans.org/3d/Bos3d_BldgModels_20230927_OBJ"
BASE_GROUND_URL = "https://maps.bostonplans.org/3d/Bos3d_Terrain_2011_OBJ"


def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2) + 1):
        yield chr(c)

letters = char_range("A", "O")
nums = range(1,13)  # range(1, 13)

class BostonModelDownloader:
    def __init__(
        self,
        in_model_dir: Union[Path, str],
        out_dataset_dir: Union[Path, str] = Path("bostontwin", "boston3d"),
    ) -> None:
        # Initialize the model from existing local directories
        self.in_model_dir = in_model_dir
        self.out_dataset_dir = out_dataset_dir

        self.flat = (
            True  # for now, we don't support ground elevation different from zero
        )

        self.set_local_projections()

        print(f"Data will be downloaded from {PUBLIC_URL}.")

        self.tiles_dict_path = self.in_model_dir.joinpath("tiles_info.json")
        if self.tiles_dict_path.is_file():
            print(f"Tile dict found: {self.tiles_dict_path}\nLoading..")
            with open(self.tiles_dict_path, "r") as f:
                self.tiles_dict = json.load(f)
                self.tiles_dict = {
                    k: {
                        kk: (
                            Path(vv) if (isinstance(vv, str) and ("path" in kk)) else vv
                        )
                        for kk, vv in v.items()
                    }
                    for k, v in self.tiles_dict.items()
                }
        else:
            self.tiles_dict = self._enumerate_tiles()
        self.n_tiles = len(self.tiles_dict)

    def download_data(self, save_dir: Union[Path, str], extract_objs=True) -> None:
        if self.tiles_dict_path.is_file():
            print(
                f"Tile dict already exists in {self.tiles_dict_path}. Delete it if you want to download the dataset again."
            )
            return

        if isinstance(save_dir, str):
            save_dir = Path(save_dir)
        if not save_dir.is_dir():
            save_dir.mkdir(parents=True, exist_ok=True)

        # %% projection file
        # source: https://www.cityschema.org/tile_scheme/index.htm
        proj_url = "https://cityschema.github.io/repository-catalog/Bos3d_CityWide_Data/Bos3d_TIleGrid/Metro_Boston_3D_CRS.zip"

        zip_proj_file_path = save_dir.joinpath("Metro_Boston_3D_CRS.zip")
        proj_file_path = save_dir.joinpath("Metro_Boston_3D_CRS.prj")

        if not proj_file_path.is_file():
            r = requests.get(proj_url, stream=True, headers={"User-Agent": "'XYZ/3.0'"})
            if not r.status_code == 404:
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
        for let in letters:
            for n in nums:
                filename = f"BOS_{let}_{n}_BldgModels_OBJ"
                out_tile_dir = save_dir.joinpath(filename)
                zip_file_path = save_dir.joinpath(filename + ".zip")
                if Path(out_tile_dir).is_dir():
                    print(f"{out_tile_dir} already downloaded. Skipping.")
                else:
                    if not zip_file_path.is_file():
                        try:
                            url = BASE_MODEL_URL + "/" + filename + ".zip"

                            r = requests.get(
                                url, stream=True, headers={"User-Agent": "XYZ/3.0"}
                            )
                            if r.status_code == 404:
                                continue
                            print("Downloading " + filename + "...")

                            with open(zip_file_path, "wb") as fd:
                                for chunk in r.iter_content(chunk_size=128):
                                    fd.write(chunk)

                            print("Extracting..")
                            with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
                                zip_ref.extractall(out_tile_dir)
                            zip_file_path.unlink()
                        except FileNotFoundError:
                            print(url)
                            continue
                downloaded.append(out_tile_dir)

                ground_name = f"BOS_{let}_{n}_TerrainMesh_2011_OBJ"
                out_ground_path = save_dir.joinpath(ground_name)
                zip_ground_path = save_dir.joinpath(ground_name + ".zip")
                if out_ground_path.is_dir():
                    print(f"{ground_name} already downloaded. Skipping.")
                else:
                    if not zip_ground_path.is_file():
                        try:
                            url = BASE_GROUND_URL + "/" + ground_name + ".zip"
                            print(url)
                            r = requests.get(
                                url, headers={"User-Agent": "XYZ/3.0"}, stream=True
                            )
                            if r.status_code == 404:
                                continue
                            print("Downloading " + ground_name + "...")
                            with open(zip_ground_path, "wb") as fd:
                                for chunk in r.iter_content(chunk_size=128):
                                    fd.write(chunk)

                            print("Extracting..")
                            with zipfile.ZipFile(zip_ground_path, "r") as zip_ref:
                                zip_ref.extractall(out_ground_path)
                            zip_ground_path.unlink()
                        except FileNotFoundError:
                            print(url)
                            continue

        # %% extract single OBJ models
        if extract_objs:
            for filepath in downloaded:
                # extract building models
                model_zip_folder = filepath.joinpath("objz")
                if model_zip_folder.is_dir():
                    print(f"Extracting individual models from {model_zip_folder}..")
                    for model_zip in model_zip_folder.iterdir():
                        model_name = model_zip.stem.replace("_OBJ", "")
                        if model_zip_folder.joinpath(model_name + ".obj").is_file():
                            continue

                        if model_zip.suffix == ".zip":
                            with zipfile.ZipFile(
                                str(model_zip.resolve()), "r"
                            ) as zip_ref:
                                zip_ref.extractall(str(filepath.resolve()))

        print("Preparing scene export..")
        self.set_local_projections()

        self.tiles_dict = self._enumerate_tiles()
        self.tiles = self.tiles_dict.keys()

        self.update_tiles_dict_json()

        print("Done.")

    def _enumerate_tiles(self) -> dict:
        centers_x_m = []
        centers_y_m = []
        total_n_models = 0
        tiles_dict = {}
        for tile_dir in self.in_model_dir.iterdir():
            if tile_dir.is_dir() and "_BldgModels_OBJ" in tile_dir.stem:
                tile_name = tile_dir.stem.replace("_BldgModels_OBJ", "")

                # check if the tile information geojson exists
                orig_model_catalog_path = tile_dir.joinpath("catalog_jsonp.js")
                out_model_catalog_path = orig_model_catalog_path.with_name(
                    "catalog.geojson"
                )
                with open(out_model_catalog_path, "w") as output:
                    with open(orig_model_catalog_path, "r") as input:
                        output.write(input.read()[7:])
                truncate_utf8_chars(out_model_catalog_path, 1)
                tile_catalog = gpd.GeoDataFrame.from_file(
                    out_model_catalog_path, crs="epsg:4326"
                )

                tile_info_path = tile_dir.joinpath("scene_bounds.geojson")

                # copy tile information geojson
                tile_catalog_bounds = tile_catalog.total_bounds
                tile_bounding_box = box(*tile_catalog_bounds)
                tile_info_epsg4326 = gpd.GeoDataFrame(
                    geometry=[tile_bounding_box], columns=["geometry"], crs="epsg:4326"
                )
                tile_info_epsg4326.to_file(tile_info_path, driver="GeoJSON")

                tile_catalog_bounds_list = [
                    (tile_catalog_bounds[0], tile_catalog_bounds[1]),
                    (tile_catalog_bounds[2], tile_catalog_bounds[3]),
                    (tile_catalog_bounds[0], tile_catalog_bounds[3]),
                    (tile_catalog_bounds[2], tile_catalog_bounds[1]),
                    ]
                # convert to projection CRS
                assert check_area_of_use(
                    tile_info_epsg4326.crs, self.local_crs, tile_catalog_bounds_list
                )
                tile_info_local = gdf2crs(tile_info_epsg4326, self.lonlat2local_transformer)
                tile_info_center = tile_info_local.unary_union.centroid
                tile_center_lonlat = self.local2lonlat_transformer.transform(tile_info_center.x, tile_info_center.y)

                model_list = [
                    model_path.stem
                    for model_path in tile_dir.iterdir()
                    if model_path.suffix.lower() == ".obj"
                    and "frame" not in model_path.stem.lower()
                ]

                n_models = len(model_list)
                total_n_models = total_n_models + n_models

                tile_dict = {
                    "name": tile_name,
                    "center_lon": tile_center_lonlat[0],
                    "center_lat": tile_center_lonlat[1],
                    "model_list": model_list,
                    "n_models": n_models,
                    "tile_info_path": tile_info_path,
                    "tile_catalog_path": out_model_catalog_path,
                }
                tiles_dict[tile_name] = tile_dict
        n_tiles = len(tiles_dict)
        self.n_tiles = n_tiles
        print(f"{n_tiles} imported.")

        tiles_dict["boston"] = {
            "center_x_m": np.mean(centers_x_m),
            "center_y_m": np.mean(centers_y_m),
            "n_models": total_n_models,
        }

        return tiles_dict

    def generate_dataset(self, create_xml) -> None:
        model_dir = self.out_dataset_dir.joinpath("meshes")
        if not model_dir.is_dir():
            model_dir.mkdir(exist_ok=True, parents=True)
            
        # convert all the obj files to ply
        t0 = time.perf_counter()
        dir_obj2ply(self.in_model_dir, model_dir, recursive=True, ft2m=True, center=True)
        t1 = time.perf_counter()
        print(f"Converted all OBJ files to PLY in {t1-t0:.2f} s.")

        # scenes are imported and converted to lon lat crs (epsg:4326) by default
        times = []
        boston_mitsuba_scene_dict = {
            "type": "scene",
            "integrator": {
                "type": "path",
            },
            "light": {"type": "constant"},
            "mat-itu_brick": {
                "type": "twosided",
                "bsdf": {
                    "type": "diffuse",
                    "reflectance": {
                        "type": "rgb",
                        "value": [0.401968, 0.111874, 0.086764],
                    },
                },
            },
            "mat-itu_concrete": {
                "type": "twosided",
                "bsdf": {
                    "type": "diffuse",
                    "reflectance": {
                        "type": "rgb",
                        "value": [0.539479, 0.539479, 0.539480],
                    },
                },
            },
            "mat-itu_medium_dry_ground": {
                "type": "twosided",
                "bsdf": {
                    "type": "diffuse",
                    "reflectance": {
                        "type": "rgb",
                        "value": [65 / 255, 60 / 255, 60 / 255],
                    },
                },
            },
        }
        for tile_idx, (tile_name, tile_dict) in enumerate(self.tiles_dict.items()):
            if tile_name == "boston":
                continue
            t0 = time.perf_counter()
            output_tile_scene_path = self.out_dataset_dir.joinpath(tile_name + ".xml")
            output_tile_info_path = self.out_dataset_dir.joinpath(
                tile_name + "_tileinfo" + ".geojson"
            )
            if (
                output_tile_scene_path.is_file()
                and output_tile_scene_path.with_suffix(".geojson").is_file()
                and output_tile_info_path.is_file()
            ):
                continue

            tile_model_catalog_path = tile_dict["tile_catalog_path"]
            tile_model_catalog_gdf = gpd.GeoDataFrame.from_file(
                tile_model_catalog_path
            )
            tile_model_catalog_gdf = tile_model_catalog_gdf.to_crs("epsg:4326")

            # drop the z coordinate from the geodataframe for faster processing 
            tile_model_catalog_gdf.geometry = tile_model_catalog_gdf.geometry.apply(
                lambda x: shp.force_2d(x)
            )
            if create_xml:
                # define CRS for the tile
                tile_crs = get_crs(scene_name=tile_name,
                                   scene_center_lon_lat=[tile_dict["center_lon"], tile_dict["center_lat"]])
                tile_transformer = Transformer.from_crs(self.local_crs_lonlat, tile_crs, always_xy=True)
                with open(model_dir.parent.joinpath(f"{tile_name}.wkt"), "w") as f:
                    f.write(tile_crs.to_wkt(output_axis_rule=True))

                # enumerate the models in the tile, prepare the structure for the XML, and check the data consistency
                model_list = tile_dict["model_list"]
                models_materials = []
                models_centers = []
                n_models_tile = 0
                for model_name in model_list:
                    ## There are two sources of model information: the geojson catalog and the info.json file
                    # The catalog is a geojson file with the model information for all models in the tile
                    # The info.json file is a json file with the model information for a single model
                    # We need to make sure they match. If not, there is something wrong with the data, and we skip the model

                    # read model info from geojson catalog
                    model_info_from_catalog = tile_model_catalog_gdf[
                        tile_model_catalog_gdf["Model_ID"] == model_name
                    ]

                    # read model info from info.json
                    model_info_path = tile_model_catalog_path.parent.joinpath(
                        model_name + ".json"
                    )
                    with open(model_info_path, "r") as f:
                        model_info_from_json = json.load(f)[0]

                    # check if the information matches
                    if not self.check_model_info(
                        model_info_from_catalog=model_info_from_catalog,
                        model_info_from_json=model_info_from_json,
                        model_name=model_name,
                    ):
                        continue
                    
                    ## Choose the model material. For now, we only have two materials: brick, for walls, and concrete, for everything else
                    model_struct_type = model_info_from_catalog["StructType"].values[0]
                    if model_struct_type == "Wall":
                        model_material = "mat-itu_brick"
                    else:
                        model_material = "mat-itu_concrete"
                    models_materials.append(model_material)
                        
                    # get the model center in local coordinates
                    model_translation_local = tile_transformer.transform(
                        model_info_from_json["Centr_Lon"],
                        model_info_from_json["Centr_Lat"],
                    )

                    models_centers.append(model_translation_local)
                    n_models_tile = n_models_tile + 1

                generate_mi_xml(
                    scene_name=tile_name,
                    models_list=model_list,
                    models_dir=model_dir,
                    out_dir=self.out_dataset_dir,
                    models_materials=models_materials,
                    models_center=models_centers,
                    create_ground=True,
                )

                self.tiles_dict[tile_name]["n_models"] = n_models_tile

                self.update_tiles_dict_json()

            tile_model_catalog_gdf.to_file(
                output_tile_scene_path.with_suffix(".geojson"), driver="GeoJSON"
            )

            tile_info = gpd.GeoDataFrame.from_file(
                self.tiles_dict[tile_name]["tile_info_path"]
            )
            # tile_info["Centr_X_Ft"] = self.tiles_dict[tile_name]["center_x_ft"]
            # tile_info["Centr_Y_Ft"] = self.tiles_dict[tile_name]["center_y_ft"]
            # tile_info["Centr_X_m"] = self.tiles_dict[tile_name]["center_x_m"]
            # tile_info["Centr_Y_m"] = self.tiles_dict[tile_name]["center_y_m"]
            tile_info["center_lon"] = self.tiles_dict[tile_name]["center_lon"]
            tile_info["center_lon"] = self.tiles_dict[tile_name]["center_lat"]
            tile_info["n_models"] = n_models_tile
            # tile_info["n_triangles"] = tile_n_tri
            tile_info.to_file(output_tile_info_path, driver="GeoJSON")

            print(
                f"Tile {tile_name} imported. There were {n_models_tile} models."
            )
            t1 = time.perf_counter()
            print_eta(t0, t1, times, tile_idx, self.n_tiles)
        output_boston_scene_path = self.out_dataset_dir.joinpath("boston" + ".xml")
        mi.xml.dict_to_xml(
            boston_mitsuba_scene_dict, str(output_boston_scene_path.resolve())
        )

        output_boston_gdf_path = self.out_dataset_dir.joinpath("boston" + ".geojson")
        self._aggregate_geojson(output_boston_gdf_path)

    def _aggregate_geojson(self, out_path):
        scene_gdf_list = []
        for tile_name, tile_dict in self.tiles_dict.items():
            if tile_name=="boston":
                continue
            tile_model_catalog_path = tile_dict["tile_catalog_path"]
            tile_model_catalog_gdf = gpd.GeoDataFrame.from_file(tile_model_catalog_path)
            tile_model_catalog_gdf.to_crs("epsg:4326", inplace=True)
            tile_model_catalog_gdf.geometry = tile_model_catalog_gdf.geometry.apply(lambda x: shp.force_2d(x))
            scene_gdf_list.append(tile_model_catalog_gdf)
        boston_gdf = gpd.GeoDataFrame(
            pd.concat(scene_gdf_list, ignore_index=True), crs="epsg:4326"
        )
        boston_gdf.to_file(out_path, driver="GeoJSON")

    def set_local_projections(self):
        original_prj_path = self.in_model_dir.joinpath("Metro_Boston_3D_CRS.prj")
        if not original_prj_path.is_file():
            print("Original file not found. Please download it with download_data().")
            self.original_crs = None
        else:
            with open(original_prj_path, "r") as f:
                prj_str = f.readline()
                self.original_crs = pyproj.CRS.from_wkt(prj_str)

        self.local_crs = pyproj.CRS.from_user_input("EPSG:26986")
        self.local_crs_lonlat = pyproj.CRS.from_user_input("EPSG:4326")

        self.lonlat2local_transformer = pyproj.Transformer.from_crs(self.local_crs_lonlat, self.local_crs, always_xy=True)
        self.local2lonlat_transformer = pyproj.Transformer.from_crs(
            self.local_crs, self.local_crs_lonlat, always_xy=True)

    def _get_tile_info(self):
        self.model_rootdir

    def update_tiles_dict_json(
        self,
    ):
        tmp_tiles_dict = self.tiles_dict.copy()
        tmp_tiles_dict = {
            k: {kk: str(vv) if isinstance(vv, Path) else vv for kk, vv in v.items()}
            for k, v in tmp_tiles_dict.items()
        }
        for v in tmp_tiles_dict.values():
            v.pop("tile_info_gdf", "")

        with open(self.tiles_dict_path, "w") as f:
            json.dump(tmp_tiles_dict, f, indent=4)

    @staticmethod
    def check_model_info(model_info_from_catalog, model_info_from_json, model_name):
        ## Read model info
        if model_info_from_json["Status"] != "Current":
            return False

        ## Check elevation and center of the model
        # get ground elevation [ft] from catalog and info.json
        ground_el_from_catalog = model_info_from_catalog["Gnd_El_Ft"].values[0]
        ground_el_from_info = model_info_from_json["Gnd_El_Ft"]

        # if the ground elevation is not available, skip the model
        if np.isnan(ground_el_from_catalog) or not ground_el_from_info:
            return False
        assert (
            ground_el_from_info == ground_el_from_catalog
        ), f"Mismatch between catalog (height: ({ground_el_from_catalog})) and info.json (height: ({ground_el_from_info})).\n\tModel: {model_name}.\n\tStatus:{model_info_from_json['Status']}"

        # get center of the model [ft] from catalog and info.json
        center_x_ft_from_catalog = model_info_from_catalog["Centr_X_Ft"].values[0]
        center_y_ft_from_catalog = model_info_from_catalog["Centr_Y_Ft"].values[0]
        center_x_ft_from_info = model_info_from_json["Centr_X_Ft"]
        center_y_ft_from_info = model_info_from_json["Centr_Y_Ft"]
        z_min_ft_from_info = model_info_from_json["Z_MIn_Ft"]
        assert (center_x_ft_from_catalog == center_x_ft_from_info) and (
            center_y_ft_from_catalog == center_y_ft_from_info
        ), f"Mismatch between catalog (center: ({center_x_ft_from_catalog}, {center_y_ft_from_catalog})) and info.json (center: ({center_x_ft_from_info}, {center_y_ft_from_info}))"

        return True
