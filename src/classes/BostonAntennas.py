import geopandas as gpd
from typing import Union
from pathlib import Path
from src.utils.geo_utils import gdf2localcrs


class BostonAntennas:
    def __init__(self, dataset_dir: Union[str, Path]) -> None:
        if isinstance(dataset_dir,str):
            dataset_dir = Path(dataset_dir)
        self.dataset_dir = dataset_dir
        self.antennas_path = dataset_dir.joinpath("antennas.geojson")

        self.antenna_gdf_epsg4326 = gpd.read_file(self.antennas_path)  # crs="epsg:4326"
        # translate to the BDPA reference system
        self.antenna_gdf_local_crs = gdf2localcrs(self.antenna_gdf_epsg4326)

    def _get_antenna_from_local_gdf(self, scene_gdf_local):
        xmin, ymin, xmax, ymax = scene_gdf_local.total_bounds
        return self.get_antenna_location_from_local_bb(
            xmin, ymin, xmax, ymax
        )

    def get_antenna_location_from_gdf(self, scene_gdf):
        scene_gdf_local_crs = gdf2localcrs(scene_gdf)
        return self._get_antenna_from_local_gdf(scene_gdf_local_crs).reset_index()

    def get_antenna_location_from_local_bb(self, xmin, ymin, xmax, ymax):
        return self.antenna_gdf_local_crs.cx[xmin:xmax, ymin:ymax].reset_index(
            drop=True
        )
