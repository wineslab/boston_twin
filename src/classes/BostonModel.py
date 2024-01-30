from pathlib import Path
from typing import Union, List, Tuple
from ..functions.constants import LOCAL_CRS_ORIGIN_GEO, LOCAL_CRS_ORIGIN_STATE
import geopandas as gpd
import mitsuba as mi
import shapely as shp
import open3d as o3d


class BostonModel:
    def __init__(
        self,
        dataset_dir: Union[Path, str] = Path("dataset"),
    ) -> None:
        # Initialize the model from existing local directories
        if isinstance(dataset_dir, Path):
            self.dataset_dir = dataset_dir
        else:
            self.dataset_dir = Path(dataset_dir)

        self.mesh_dir = dataset_dir.joinpath("meshes")

        self.tiles_dict = self._enumerate_scenes()
        self.tile_names = list(self.tiles_dict.keys())

        self.origin_epsg4326 = LOCAL_CRS_ORIGIN_GEO
        # origin_pnt = shp.Point(LOCAL_CRS_ORIGIN_GEO)
        # origin_gdf = gpd.GeoDataFrame(
        #     index=[0], geometry=[origin_pnt], crs="epsg:4326"
        # ).to_crs("epsg:2249")
        # self.origin_epsg2249 = list(origin_gdf.values[0][0].coords)[0]

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

    def generate_area_scene(
        self, origin_lonlat: Union[Tuple, List], scene_out_name: str
    ):
        mi.xml.dict_to_xml()
        # https://mitsuba.readthedocs.io/en/latest/src/api_reference.html#properties
