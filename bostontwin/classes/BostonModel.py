import time
from pathlib import Path
from typing import Union

import geopandas as gpd
import mitsuba as mi

try:
    import open3d as o3d
except:
    print("open3d not available")
import shapely as shp

from ..utils.constants import LOCAL_CRS_ORIGIN_GEO, LOCAL_CRS_ORIGIN_STATE
from ..utils.geo_utils import gdf2crs
from ..utils.obj_utils import create_ground, get_mi_dict

try:
    mi.set_variant("llvm_ad_rgb")
except:
    mi.set_variant("scalar_rgb")


class BostonModel:
    def __init__(
        self,
        dataset_dir: Union[Path, str] = Path("dataset"),
    ) -> None:
        # Initialize the model from existing local directories
        if isinstance(dataset_dir, str):
            dataset_dir = Path(dataset_dir)
        self.dataset_dir = dataset_dir

        self.mesh_dir = dataset_dir.joinpath("meshes")

        self.tiles_dict = self._enumerate_scenes()
        self.tile_names = list(self.tiles_dict.keys())

        self.origin_epsg4326 = LOCAL_CRS_ORIGIN_GEO
        self.origin_epsg2249 = LOCAL_CRS_ORIGIN_STATE

        self.flat = (
            True  # for now, we don't support ground elevation different from zero
        )

    # def _check_in_model_dir(self):

    def _enumerate_scenes(self) -> dict:
        scenes_dict = {}
        for scene_file_path in self.dataset_dir.iterdir():
            if (scene_file_path.is_file()) and (scene_file_path.suffix == ".geojson") and ("tileinfo" not in scene_file_path.stem):
                scene_name = scene_file_path.stem

                mi_scene_path = scene_file_path.resolve().with_suffix(".xml")
                tile_info_path = scene_file_path.resolve().with_name(
                    f"{scene_name}_tileinfo.geojson"
                )
                # scene = mi.load_file(str(mi_scene_path))
                # tile_info = gpd.read_file(tile_scene_file)
                scenes_dict[scene_name] = {
                    "tileinfo_path": tile_info_path,
                    "geo_scene_path": scene_file_path,
                    "mi_scene_path": mi_scene_path,
                    # "models": [model.id() for model in scene.shapes()],
                }

        n_scenes = len(scenes_dict)
        self.n_tiles = n_scenes
        print(f"{n_scenes} scenes imported.")

        return scenes_dict

    def convert_to_ascii(self, out_mesh_dir: Union[str, Path] = None):
        if not out_mesh_dir:
            out_mesh_dir = self.mesh_dir
        if not out_mesh_dir.is_dir():
            out_mesh_dir.mkdir(parents=True, exist_ok=True)

        for model_in_path in self.mesh_dir.iterdir():
            model_in = o3d.io.read_triangle_mesh(str(model_in_path))

            model_out_path = out_mesh_dir.joinpath(model_in_path.stem + "_ascii.ply")
            o3d.io.write_triangle_mesh(str(model_out_path), model_in, write_ascii=True)

    def generate_simplified_dataset(
        self,
        scene_in_name: str,
        out_dir: Union[str, Path],
        precision=0.5,
        smooth=False,
        n_smooth_iterations=1,
    ):
        if scene_in_name not in self.tile_names:
            raise KeyError(
                f"scene_name must correspond to one of the scenes name in {self.dataset_dir}: {self.tile_names}\nInstead: {scene_in_name}"
            )
        scene_in_dict = self.tiles_dict[scene_in_name]

        if self.dataset_dir == out_dir:
            raise FileExistsError(
                "Specify a different folder for the simplified dataset.\nOverwriting the original dataset is not supported."
            )

        out_mesh_dir = out_dir.joinpath("meshes")
        if not out_mesh_dir.is_dir():
            out_mesh_dir.mkdir(exist_ok=True, parents=True)

        n_triangles_in = 0
        n_vert_in = 0
        n_triangles_out = 0
        n_vert_out = 0
        for model_name in scene_in_dict["models"]:
            model_in_path = self.mesh_dir.joinpath(model_name + ".ply")
            model_in = o3d.io.read_triangle_mesh(str(model_in_path.resolve()))
            n_triangles_in = n_triangles_in + len(model_in.triangles)
            n_vert_in = n_vert_in + len(model_in.vertices)

            if smooth:
                model_in = model_in.filter_smooth_simple(
                    number_of_iterations=n_smooth_iterations
                )

            model_out = model_in.simplify_vertex_clustering(voxel_size=precision)
            n_triangles_out = n_triangles_out + len(model_out.triangles)
            n_vert_out = n_vert_out + len(model_out.vertices)

            model_out_path = out_mesh_dir.joinpath(model_name + ".ply")
            o3d.io.write_triangle_mesh(str(model_out_path.resolve()), model_out)

        # copy the scene file
        scene_in_path = self.dataset_dir.joinpath(scene_in_name + ".xml")

        scene_out_path = out_dir.joinpath(scene_in_path.name)
        scene_out_path.write_text(scene_in_path.read_text())

        scene_out_path = out_dir.joinpath(scene_in_path.with_suffix(".geojson").name)
        scene_out_path.write_text(scene_in_path.with_suffix(".geojson").read_text())

        print(f"Input scene has {n_vert_in} vertices and {n_triangles_in} triangles.")
        print(
            f"Output scene has {n_vert_out} vertices and {n_triangles_out} triangles."
        )
        print(f"Saved in {scene_out_path}")

    def generate_scene_from_model_gdf(self, model_gdf, scene_center, scene_name) -> None:
        out_scene_path = self.dataset_dir
        meshes_path = self.mesh_dir

        output_scene_path = self.dataset_dir.joinpath(scene_name + ".xml")
        gdf_out_path = output_scene_path.with_suffix(".geojson")
        if not gdf_out_path.is_file():
            model_gdf.to_file(gdf_out_path, driver="GeoJSON")

        output_scene_info_path = self.dataset_dir.joinpath(
            scene_name + "_tileinfo" + ".geojson"
        )

        xmin, ymin, xmax, ymax = model_gdf.total_bounds
        tile_info = gpd.GeoDataFrame(
            geometry=[
                shp.Polygon(
                    [
                        (xmin, ymin),
                        (xmax, ymin),
                        (xmax, ymax),
                        (xmin, ymax),
                        (xmin, ymin),
                    ]
                )
            ],
            columns=["geometry"],
            crs="epsg:4326"
        )

        scene_crs = get_crs(scene_name=scene_name,
                            scene_center_lon_lat=[scene_center["center_lon"], scene_center["center_lat"]])
        scene_transformer = Transformer.from_crs("EPSG:4326", scene_crs, always_xy=True)
        with open(model_dir.parent.joinpath(f"{scene_name}.wkt"), "w") as f:
            f.write(scene_crs.to_wkt(output_axis_rule=True))

        model_gdf = gdf2crs(model_gdf)
        xmin, ymin, xmax, ymax = model_gdf.total_bounds
        scene_size_x = xmax-xmin
        scene_size_y = ymax-ymin
        print(f"Scene size is {scene_size:.2f}x{scene_size:.2f} m.")

        model_dict = out_model_dict.copy()

        scene_dict_mi["ground"] = model_dict

        model_gdf["center"] = model_gdf.centroid
        
        n_models_scene = model_gdf.shape[0]
        
        model_list = model_gdf["Model_ID"].tolist()
        model_materials = model_gdf["StructType"].apply(lambda x: "mat-itu_brick" if "Wall" in x else "mat-itu_concrete").tolist()
        model_centers = model_gdf.centroid.apply(lambda x:scene_transformer.transform(x)).tolist()

        generate_mi_xml(
            scene_name=scene_name,
            models_list=model_list,
            models_dir=model_dir,
            out_dir=self.out_dataset_dir,
            models_materials=models_materials,
            models_center=models_centers,
            create_ground=True,
            ground_size= max(scene_size_x,scene_size_y)
        )
        tile_info["Centr_X_m"] = scene_center_local[0]
        tile_info["Centr_Y_m"] = scene_center_local[1]
        tile_info["Centr_lon"] = scene_center[0]
        tile_info["Centr_lat"] = scene_center[1]
        tile_info["n_models"] = n_models_scene
        tile_info.to_file(output_scene_info_path, driver="GeoJSON")

        print(
            f"Scene {scene_name} imported. There were {n_models_scene} models ({scene_n_tri} triangles)."
        )
        self.tiles_dict = self._enumerate_scenes()
        self.tile_names = list(self.tiles_dict.keys())
