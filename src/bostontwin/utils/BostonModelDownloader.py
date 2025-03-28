import json
import time
import zipfile
import logging
from pathlib import Path
from typing import Union

from copy import deepcopy

import geopandas as gpd
import mitsuba as mi
import numpy as np
import pandas as pd
import pyproj
import requests
import shapely as shp
from shapely.geometry import box

from src.bostontwin.utils.geo_utils import check_area_of_use, gdf2crs, get_crs
from src.bostontwin.utils.obj_utils import (
    create_ground_dict,
    dir_obj2ply,
    obj2ply_crs_conversion,
    read_mesh,
    merge_meshes,
    save_mesh,
)
from src.bostontwin.utils.utils import print_eta, truncate_utf8_chars
from src.bostontwin.utils.mi_utils import create_mitsuba_xml, get_transformation_matrix
from typing import List, Tuple

mi.set_variant("scalar_rgb")

logging.basicConfig(level=logging.INFO)

materials_str2dict = {
    "wet_ground": {
        "type": "wet_ground",
        "id": "wet_ground",
        "thickness": 100,
    },
    "medium_dry_ground": {
        "type": "medium_dry_ground",
        "id": "medium_dry_ground",
        "thickness": 100,
    },
    "very_dry_ground": {
        "type": "very_dry_ground",
        "id": "very_dry_ground",
        "thickness": 100,
    },
    "brick": {
        "type": "brick",
        "id": "brick",
        "thickness": 0.4,
    },
    "concrete": {
        "type": "concrete",
        "id": "concrete",
        "thickness": 1.0,
    },
}

PUBLIC_URL = (
    "https://www.bostonplans.org/3d-data-maps/3d-smart-model/3d-data-download"
)
BASE_MODEL_URL = "https://maps.bostonplans.org/3d/Bos3d_BldgModels_20250128_OBJ"
BASE_GROUND_URL = "https://maps.bostonplans.org/3d/Bos3d_Terrain_2011_OBJ"

NU_URL = "https://repository.library.northeastern.edu/downloads/neu:ms36tq790?datastream_id=content"

data_path = Path(__file__).parents[3].joinpath("data")
logging.info("Data path: %s", data_path)


def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2) + 1):
        yield chr(c)


letters = char_range("A", "O")
nums = range(1, 13)  # range(1, 13)


class BostonModelDownloader:
    def __init__(
        self,
        in_model_dir: Union[Path, str],
        out_dataset_dir: Union[Path, str] = Path("dataset", "scenes"),
    ) -> None:
        # Initialize the model from existing local directories
        self.in_model_dir = in_model_dir
        self.out_dataset_dir = out_dataset_dir

        self.bostontwin_center = (
            -71.08765495983191,
            42.337479190130736,
        )  # center of the local CRS
        self.boston_bounds = [-71.187459, 42.240483, -70.927061, 42.390985]

        self.set_local_projections()

        logging.info(f"Data will be downloaded from {PUBLIC_URL}.")

        self.tiles_dict_path = self.in_model_dir.joinpath("tiles_info.json")
        if self.tiles_dict_path.is_file():
            logging.info(f"Tile dict found: {self.tiles_dict_path.resolve()}\nLoading..")
            with open(self.tiles_dict_path, "r") as f:
                self.tiles_dict = json.load(f)
                self.tiles_dict = {
                    k: {
                        kk: (
                            Path(vv)
                            if (isinstance(vv, str) and ("path" in kk))
                            else vv
                        )
                        for kk, vv in v.items()
                    }
                    for k, v in self.tiles_dict.items()
                }
        else:
            self.tiles_dict = self._enumerate_tiles()
        self.n_tiles = len(self.tiles_dict)

    def download_data(
        self, save_dir: Union[Path, str], extract_objs=True
    ) -> None:
        # try:
        # TODO: change back to try
        if False:
            logging.info(
                "Downloading the dataset from the Northeastern repository.."
            )
            zip_dataset_path = save_dir.joinpath("BostonTwinDataset.zip")
            r = requests.get(
                NU_URL, stream=True, headers={"User-Agent": "'XYZ/3.0'"}
            )
            if not r.status_code == 404:
                with open(zip_dataset_path, "wb") as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)

                logging.info("Extracting..")
                with zipfile.ZipFile(zip_dataset_path, "r") as zip_ref:
                    zip_ref.extractall(self.out_dataset_dir)
                zip_dataset_path.unlink()

                return
        else:
            # except FileNotFoundError as e:
            # logging.info(f"Can't download from the Northeastern repository. Trying the BPDA website. ({e})")
            pass

        if self.tiles_dict_path.is_file():
            logging.info(
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
            r = requests.get(
                proj_url, stream=True, headers={"User-Agent": "'XYZ/3.0'"}
            )
            if not r.status_code == 404:
                logging.info("Downloading the 3D projection file...")

                with open(zip_proj_file_path, "wb") as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)

                logging.info("Extracting..")
                with zipfile.ZipFile(zip_proj_file_path, "r") as zip_ref:
                    zip_ref.extractall(save_dir)
                zip_proj_file_path.unlink()

        # %% Model files
        downloaded = []
        for let in letters:
            for n in nums:
                if "G" in let and n == 3:
                    continue
                filename = f"BOS_{let}_{n}_BldgModels_OBJ"
                out_tile_dir = save_dir.joinpath(filename)
                zip_file_path = save_dir.joinpath(filename + ".zip")
                url = BASE_MODEL_URL + "/" + filename + ".zip"
                if Path(out_tile_dir).is_dir():
                    logging.info(
                        f"{out_tile_dir} already downloaded. Skipping."
                    )
                else:
                    if not zip_file_path.is_file():
                        try:
                            r = requests.get(
                                url,
                                stream=True,
                                headers={"User-Agent": "XYZ/3.0"},
                            )
                            if r.status_code == 404:
                                continue
                            logging.info("Downloading " + filename + "...")

                            with open(zip_file_path, "wb") as fd:
                                for chunk in r.iter_content(chunk_size=128):
                                    fd.write(chunk)

                            logging.info("Extracting..")
                            with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
                                zip_ref.extractall(out_tile_dir)
                            zip_file_path.unlink()
                        except FileNotFoundError:
                            logging.info(url)
                            continue
                downloaded.append(out_tile_dir)

                ground_name = f"BOS_{let}_{n}_TerrainMesh_2011_OBJ"
                out_ground_path = save_dir.joinpath(ground_name)
                zip_ground_path = save_dir.joinpath(ground_name + ".zip")
                if out_ground_path.is_dir():
                    logging.info(f"{ground_name} already downloaded. Skipping.")
                else:
                    if not zip_ground_path.is_file():
                        try:
                            url = BASE_GROUND_URL + "/" + ground_name + ".zip"
                            logging.info(url)
                            r = requests.get(
                                url,
                                headers={"User-Agent": "XYZ/3.0"},
                                stream=True,
                            )
                            if r.status_code == 404:
                                continue
                            logging.info("Downloading " + ground_name + "...")
                            with open(zip_ground_path, "wb") as fd:
                                for chunk in r.iter_content(chunk_size=128):
                                    fd.write(chunk)

                            logging.info("Extracting..")
                            with zipfile.ZipFile(
                                zip_ground_path, "r"
                            ) as zip_ref:
                                zip_ref.extractall(out_ground_path)
                            zip_ground_path.unlink()
                        except FileNotFoundError:
                            logging.info(url)
                            continue

        # %% extract single OBJ models
        if extract_objs:
            for filepath in downloaded:
                # extract building models
                model_zip_folder = filepath.joinpath("objz")
                if model_zip_folder.is_dir():
                    logging.info(
                        f"Extracting individual models from {model_zip_folder}.."
                    )
                    for model_zip in model_zip_folder.iterdir():
                        model_name = model_zip.stem.replace("_OBJ", "")
                        if model_zip_folder.joinpath(
                            model_name + ".obj"
                        ).is_file():
                            continue

                        if model_zip.suffix == ".zip" and model_zip.stat().st_size>0:
                            try:
                                with zipfile.ZipFile(
                                    str(model_zip.resolve()), "r"
                                ) as zip_ref:
                                    zip_ref.extractall(str(filepath.resolve()))
                            except zipfile.BadZipFile as e:
                                logging.error(
                                    f"Bad zip file: {model_zip}. {e}"
                                )
                                continue

        logging.info("Preparing scene export..")
        self.set_local_projections()

        self.tiles_dict = self._enumerate_tiles()
        self.tiles = self.tiles_dict.keys()

        self.update_tiles_dict_json()

        logging.info("Done.")

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
                    geometry=[tile_bounding_box],
                    columns=["geometry"],
                    crs="epsg:4326",
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
                    tile_info_epsg4326.crs,
                    self.local_crs,
                    tile_catalog_bounds_list,
                )
                tile_info_local = gdf2crs(
                    tile_info_epsg4326, self.lonlat2local_transformer
                )
                tile_info_center = tile_info_local.unary_union.centroid
                tile_center_lonlat = self.local2lonlat_transformer.transform(
                    tile_info_center.x, tile_info_center.y
                )

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
        logging.info(f"{n_tiles} imported.")

        tiles_dict["boston"] = {
            "center_x_m": np.mean(centers_x_m),
            "center_y_m": np.mean(centers_y_m),
            "n_models": total_n_models,
        }

        return tiles_dict

    def generate_dataset(self) -> None:
        logging.info("Starting the scene generation..")

        model_dir = self.out_dataset_dir.joinpath("meshes")
        if not model_dir.is_dir():
            model_dir.mkdir(exist_ok=True, parents=True)

        # convert all the obj files to ply
        t0 = time.perf_counter()
        meshes_info = dir_obj2ply(
            self.in_model_dir,
            model_dir,
            recursive=True,
            transformer=self.original2local_transformer,
        )
        t1 = time.perf_counter()
        logging.info(f"Converted all OBJ files to PLY in {t1 - t0:.2f} s.")

        # scenes are imported and converted to lon lat crs (epsg:4326) by default
        boston_n_triangles = []
        boston_n_models = []
        boston_models = []
        boston_terrains = []
        boston_terrains_flat = []
        boston_materials = set()
        times = []
        ground_material = "medium_dry_ground"
        for tile_idx, (tile_name, tile_dict) in enumerate(
            self.tiles_dict.items()
        ):
            if tile_name == "boston":
                continue
            t0 = time.perf_counter()

            logging.info(f"Importing scene {tile_name}..")

            tile_model_catalog_path = tile_dict["tile_catalog_path"]
            tile_model_catalog_gdf = gpd.GeoDataFrame.from_file(
                tile_model_catalog_path
            )

            tile_model_catalog_gdf, tile_center = self.preprocess_geojson(
                lonlat2local_transformer=self.lonlat2local_transformer,
                in_geojson=tile_model_catalog_gdf,
            )

            scene_path = self.out_dataset_dir.joinpath(tile_name + ".xml")

            # enumerate the models in the tile, prepare the structure for the XML, and check the data consistency
            model_list = tile_model_catalog_gdf["Model_ID"].tolist()

            models_centers = []
            triangles_list = []
            scene_mesh_info_list = []
            n_models_tile = 0
            models_materials = set()
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
                    model_info_from_json = json.load(f)
                
                if len(model_info_from_json)==1:
                    model_info_from_json = model_info_from_json[0]
                else:
                    for m_info in model_info_from_json:
                        if m_info["Status"]==model_info_from_catalog["Status"].values[0]:
                            model_info_from_json = m_info
                            break

                # check if the information matches
                if not self.check_model_info(
                    model_info_from_catalog=model_info_from_catalog,
                    model_info_from_json=model_info_from_json,
                    model_name=model_name,
                ):
                    raise ValueError(
                        f"Model {model_name} has inconsistent information."
                    )

                model_mesh_info = meshes_info[model_name]
                model_center = model_mesh_info["center"]
                model_n_tri = model_mesh_info["n_tri"]

                ## Choose the model material. For now, we only have two materials: brick, for walls, and concrete, for everything else
                model_struct_type = model_info_from_catalog[
                    "StructType"
                ].values[0]
                if model_struct_type == "Wall":
                    model_material = "brick"
                else:
                    model_material = "concrete"
                models_materials.add(model_material)
                model_mesh_info["material_id"] = model_material

                tile_model_catalog_gdf.loc[
                    tile_model_catalog_gdf["Model_ID"] == model_name,
                    "Center_X_local",
                ] = model_center[0]
                tile_model_catalog_gdf.loc[
                    tile_model_catalog_gdf["Model_ID"] == model_name,
                    "Center_Y_local",
                ] = model_center[1]
                tile_model_catalog_gdf.loc[
                    tile_model_catalog_gdf["Model_ID"] == model_name,
                    "Center_Z_local",
                ] = model_center[2]
                tile_model_catalog_gdf.loc[
                    tile_model_catalog_gdf["Model_ID"] == model_name,
                    "n_triangles",
                ] = model_n_tri

                scene_mesh_info_list.append(model_mesh_info)
                triangles_list.append(model_n_tri)
                models_centers.append(model_center)
                n_models_tile = n_models_tile + 1

            # add information to the global Boston scene before creating the ground model
            boston_n_triangles.append(sum(triangles_list))
            boston_n_models.append(n_models_tile)
            boston_models.extend(scene_mesh_info_list)
            boston_materials.update(models_materials)

            for flat in [False, True]:
                scene_name = tile_name if not flat else tile_name + "_flat"
                # create the ground model
                # TODO: create both flat and terrain ground models
                ground_dict = self.add_ground_model(
                    flat=flat,
                    ground_material=ground_material,
                    tile_name=tile_name,
                    model_dir=model_dir,
                    obj_path=self.in_model_dir,
                    tile_center=tile_center,
                    transformer=self.original2local_transformer,
                )
                ground_dict["material_id"] = ground_material

                if flat:
                    boston_terrains_flat.append(ground_dict)
                else:
                    boston_terrains.append(ground_dict)

                scene_mesh_info_list.append(ground_dict)
                triangles_list.append(ground_dict["n_tri"])
                n_models_tile = n_models_tile + 1
                models_materials.add(
                    ground_material
                )  # add the ground material to the list of materials

                # create the mitsuba XML
                scene_path = self.out_dataset_dir.joinpath(scene_name + ".xml")
                models_materials_list = [
                    materials_str2dict[k] for k in models_materials
                ]
                if flat:
                    scene_mesh_info_list_flat = deepcopy(scene_mesh_info_list)
                    for mesh_info in scene_mesh_info_list_flat:
                        mesh_info.update(
                            {
                                "to_world": get_transformation_matrix(
                                    translation=[0, 0, -mesh_info["z_min"]]
                                )
                            }
                        )

                    create_mitsuba_xml(
                        scene_path,
                        scene_mesh_info_list_flat,
                        models_materials_list,
                    )
                else:
                    create_mitsuba_xml(
                        scene_path, scene_mesh_info_list, models_materials_list
                    )

                self.tiles_dict[tile_name]["n_models"] = n_models_tile

                # update the cached tile information
                self.update_tiles_dict_json()

                # save the catalog geojson in lonlat coordinates
                output_tile_scene_path = self.out_dataset_dir.joinpath(
                    scene_name
                )
                tile_model_catalog_gdf.to_file(
                    output_tile_scene_path.with_suffix(".geojson"),
                    driver="GeoJSON",
                )

                # save the tile information in a dedicated geojson for faster access
                tile_info = gpd.GeoDataFrame.from_file(
                    self.tiles_dict[tile_name]["tile_info_path"]
                )
                tile_info["center_lon"] = self.tiles_dict[tile_name][
                    "center_lon"
                ]
                tile_info["center_lat"] = self.tiles_dict[tile_name][
                    "center_lat"
                ]
                tile_info["center_x_local"] = tile_center[0]
                tile_info["center_y_local"] = tile_center[1]
                tile_info["n_models"] = n_models_tile
                tile_info["n_triangles"] = sum(triangles_list)
                output_tile_info_path = self.out_dataset_dir.joinpath(
                    scene_name + "_tileinfo" + ".geojson"
                )
                tile_info.to_file(output_tile_info_path, driver="GeoJSON")
                
                del scene_mesh_info_list[-1]
                del triangles_list[-1]
                n_models_tile = n_models_tile - 1
                del models_materials_list[-1]

                logging.info(
                    f"Tile {tile_name} imported. ({n_models_tile} models, {sum(triangles_list)} triangles)"
                )

                t1 = time.perf_counter()
                print_eta(t0, t1, times, tile_idx, self.n_tiles)
                
        # create the global Boston terrain model
        for flat in [False, True]:
            if flat:
                suffix = "_flat.ply"
                boston_terrains_path = [model_dir.parent.joinpath(Path(terrain["filename"]).with_name(Path(terrain["filename"]).stem + suffix)) for terrain in boston_terrains_flat]
            else:
                suffix = ".ply"
                boston_terrains_path = [model_dir.parent.joinpath(Path(terrain["filename"]).with_name(Path(terrain["filename"]).stem + suffix)) for terrain in boston_terrains]
            global_terrain_model_path = model_dir.joinpath(f"boston_terrain{suffix}")

            self.create_global_terrain_model(
                boston_terrains_path, global_terrain_model_path
            )
            global_ground_dict = create_ground_dict(
                model_material=ground_material,
                x_shift=0,
                y_shift=0,
                z_shift=0,
                scale_factor=1,
                base_rect_path=global_terrain_model_path,
                ply_path=global_terrain_model_path,
            )
            global_ground_dict["id"] = global_terrain_model_path.stem
            boston_models.append(global_ground_dict)

            # create the Mitsuba XML for the global Boston
            boston_materials = [materials_str2dict[k] for k in boston_materials]
            output_boston_scene_path = self.out_dataset_dir.joinpath(
                "boston" + ".xml"
            )
            create_mitsuba_xml(
                output_boston_scene_path, boston_models, boston_materials
            )

        # save the tile information in a dedicated geojson for faster access
        boston_info = gpd.GeoDataFrame(
            geometry=[box(*self.boston_bounds)],
            columns=["geometry"],
            crs="epsg:4326",
        )
        boston_info["center_lon"] = self.bostontwin_center[0]
        boston_info["center_lat"] = self.bostontwin_center[1]
        boston_info["n_models"] = sum(boston_n_models)
        boston_info["n_triangles"] = sum(boston_n_triangles)
        boston_info.to_file(
            self.out_dataset_dir.joinpath("boston_tileinfo.geojson"),
            driver="GeoJSON",
        )

        # save the catalog geojson in lonlat coordinates
        output_boston_gdf_path = self.out_dataset_dir.joinpath(
            "boston" + ".geojson"
        )
        self._aggregate_geojson(output_boston_gdf_path, boston_models)
        logging.info("Done. You can now use the BostonTwin.")

    def create_global_terrain_model(
        self, boston_terrains: List[Union[Path,str]], out_path: Union[Path, str]
    ) -> None:
        terrain_meshes = []
        for terrain_path in boston_terrains:
            terrain_meshes.append(read_mesh(terrain_path))
        terrain_mesh = merge_meshes(terrain_meshes)
        save_mesh(terrain_mesh, out_path)

    def _aggregate_geojson(
        self, out_path: Union[Path, str], model_list: List = []
    ):
        scene_gdf_list = []
        for tile_name, tile_dict in self.tiles_dict.items():
            if tile_name == "boston":
                continue
            tile_model_catalog_path = tile_dict["tile_catalog_path"]
            tile_model_catalog_gdf = gpd.GeoDataFrame.from_file(
                tile_model_catalog_path
            )
            if model_list:
                tile_model_catalog_gdf = tile_model_catalog_gdf[
                    tile_model_catalog_gdf["Model_ID"].isin(model_list)
                ]
            tile_model_catalog_gdf.to_crs("epsg:4326", inplace=True)
            tile_model_catalog_gdf.geometry = (
                tile_model_catalog_gdf.geometry.apply(lambda x: shp.force_2d(x))
            )
            scene_gdf_list.append(tile_model_catalog_gdf)
        boston_gdf = gpd.GeoDataFrame(
            pd.concat(scene_gdf_list, ignore_index=True), crs="epsg:4326"
        )
        boston_gdf.to_file(out_path, driver="GeoJSON")

    def set_local_projections(self):
        original_prj_path = (
            Path(__file__)
            .parents[3]
            .joinpath("data")
            .joinpath("Metro_Boston_3D_CRS.prj")
        )
        if not original_prj_path.is_file():
            raise FileNotFoundError(
                f"Projection file not found: {original_prj_path}"
            )
        else:
            with open(original_prj_path, "r") as f:
                prj_str = f.readline()
                self.original_crs = pyproj.CRS.from_wkt(prj_str)

        self.local_crs = get_crs(
            scene_name="BostonTwin", scene_center_lon_lat=self.bostontwin_center
        ).to_3d()
        with open(
            self.out_dataset_dir.parent.joinpath("BostonTwin.wkt"), "w"
        ) as f:
            f.write(self.local_crs.to_wkt(output_axis_rule=True))

        self.local_crs_lonlat = pyproj.CRS.from_user_input("EPSG:4326")

        self.original2local_transformer = pyproj.Transformer.from_crs(
            self.original_crs,
            self.local_crs,
            always_xy=True,
            allow_ballpark=False,
        )

        self.lonlat2local_transformer = pyproj.Transformer.from_crs(
            self.local_crs_lonlat, self.local_crs, always_xy=True
        )
        self.local2lonlat_transformer = pyproj.Transformer.from_crs(
            self.local_crs, self.local_crs_lonlat, always_xy=True
        )

    def update_tiles_dict_json(
        self,
    ):
        tmp_tiles_dict = self.tiles_dict.copy()
        tmp_tiles_dict = {
            k: {
                kk: str(vv) if isinstance(vv, Path) else vv
                for kk, vv in v.items()
            }
            for k, v in tmp_tiles_dict.items()
        }
        for v in tmp_tiles_dict.values():
            v.pop("tile_info_gdf", "")

        with open(self.tiles_dict_path, "w") as f:
            json.dump(tmp_tiles_dict, f, indent=4)

    @staticmethod
    def check_model_info(
        model_info_from_catalog, model_info_from_json, model_name
    ):
        ## Read model info
        if (model_info_from_json["Status"] == "Approved Demo") or (
            model_info_from_json["Status"] == "History"
        ):
            logging.info(
                f"Model {model_name} is not current ({model_info_from_json['Status']}). Skipping."
            )
            return False

        # ## Check elevation and center of the model
        # # get ground elevation [ft] from catalog and info.json
        # ground_el_from_catalog = model_info_from_catalog['Gnd_El_Ft'].values[0]
        # ground_el_from_info = model_info_from_json['Gnd_El_Ft']
        #
        # # if the ground elevation is not available, skip the model
        # if np.isnan(ground_el_from_catalog) or not ground_el_from_info:
        #     logging.info(f'Model {model_name} has no ground elevation. Skipping.')
        #     return False
        # assert (
        #     abs(ground_el_from_info - ground_el_from_catalog)<1
        # ), f'Mismatch between catalog (height: ({ground_el_from_catalog})) and info.json (height: ({ground_el_from_info})).\n\tModel: {model_name}.\n\tStatus:{model_info_from_json["Status"]}'

        # get center of the model [ft] from catalog and info.json
        center_x_ft_from_catalog = model_info_from_catalog["Centr_X_Ft"].values[
            0
        ]
        center_y_ft_from_catalog = model_info_from_catalog["Centr_Y_Ft"].values[
            0
        ]
        center_x_ft_from_info = model_info_from_json["Centr_X_Ft"]
        center_y_ft_from_info = model_info_from_json["Centr_Y_Ft"]
        # z_min_ft_from_info = model_info_from_json['Z_MIn_Ft']
        if not (
            (abs(center_x_ft_from_catalog - center_x_ft_from_info) < 1)
            and (abs(center_y_ft_from_catalog - center_y_ft_from_info) < 1)
        ):
            logging.error(
                f"Mismatch between catalog (center: ({center_x_ft_from_catalog}, {center_y_ft_from_catalog})) and info.json (center: ({center_x_ft_from_info}, {center_y_ft_from_info}))"
            )
            return False

        return True

    @staticmethod
    def preprocess_geojson(
        lonlat2local_transformer: pyproj.Transformer,
        in_geojson: gpd.GeoDataFrame,
    ) -> Tuple[gpd.GeoDataFrame, List]:
        tile_bounds = in_geojson.total_bounds
        tile_bounds = [
            lonlat2local_transformer.transform(
                tile_bounds[0],
                tile_bounds[1],
            ),
            lonlat2local_transformer.transform(
                tile_bounds[2],
                tile_bounds[3],
            ),
        ]
        tile_center = [
            (tile_bounds[0][0] + tile_bounds[1][0]) / 2,
            (tile_bounds[0][1] + tile_bounds[1][1]) / 2,
        ]

        in_geojson = in_geojson.to_crs("epsg:4326")
        in_geojson = in_geojson.loc[(in_geojson["Status"] != "History") & (in_geojson["Status"] != "Approved Demo"),:]
        return in_geojson, tile_center

    @staticmethod
    def add_ground_model(
        flat: bool,
        ground_material: str,
        tile_name: str,
        model_dir: Path,
        obj_path: Union[Path, str] = None,
        tile_center: Union[Tuple, List] = None,
        transformer: pyproj.Transformer = None,
    ) -> dict:
        if flat:  # if the flat option is enabled, we create a rectangle model as ground
            if tile_center is None:
                raise ValueError(
                    "tile_center must be provided if flat is True."
                )
            print("Tile center:", tile_center)
            ply_path = model_dir.joinpath(tile_name + "_terrain_flat.ply")
            ground_dict = create_ground_dict(
                model_material=ground_material,
                x_shift=tile_center[0],
                y_shift=tile_center[1],
                z_shift=0.0,
                scale_factor=1.0,
                base_rect_path=data_path.joinpath("rectangle.ply"),
                ply_path=ply_path
            )
        else:  # if the flat option is disabled, we use the terrain mesh
            # convert terrain mesh from obj to ply
            if obj_path is None:
                raise ValueError(
                    "obj_path must be provided when flat is False to locate the terrain mesh."
                )
            if transformer is None:
                raise ValueError(
                    "transformer must be provided when flat is False to convert coordinates."
                )
            terrain_name = tile_name + "_TerrainMesh_2011_OBJ"
            terrain_obj_path = obj_path.joinpath(
                terrain_name, terrain_name.replace("_OBJ", "") + ".obj"
            )
            ply_path = model_dir.joinpath(tile_name + "_terrain.ply")
            mesh_center, ground_n_tri, z_min = obj2ply_crs_conversion(
                terrain_obj_path,
                ply_path,
                transformer=transformer,
                flat=flat,
            )
            # create the ground dict
            ground_dict = {
                "filename": str(ply_path.relative_to(model_dir.parent)),
                "center": mesh_center,
                "n_tri": ground_n_tri,
                "material_id": ground_material,
                "z_min": 0.0,
            }
        id = ply_path.stem
        ground_dict["id"] = id
        return ground_dict
