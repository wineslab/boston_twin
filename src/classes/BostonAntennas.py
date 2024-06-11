from pathlib import Path
from typing import Union

import geopandas as gpd

from src.utils.geo_utils import gdf2localcrs


class BostonAntennas:

    def __init__(
        self, dataset_dir: Union[str, Path], antennas_filename: str = "antennas.geojson"
    ) -> None:
        if isinstance(dataset_dir,str):
            dataset_dir = Path(dataset_dir)
        self.dataset_dir = dataset_dir
        self.antennas_path = dataset_dir.joinpath(antennas_filename)

        self.antenna_gdf_epsg4326 = gpd.read_file(self.antennas_path)  # crs="epsg:4326"
        # translate to the BDPA reference system
        self.antenna_gdf_local_crs = gdf2localcrs(self.antenna_gdf_epsg4326)

    def get_antenna_location_from_gdf(self, scene_gdf):
        if "WGS 84" not in scene_gdf.crs.name:
            raise ValueError(f"Input GeoDataFrame must be in longitude,latitude coordinates. Instead {scene_gdf.crs.name}")
        xmin, ymin, xmax, ymax = scene_gdf.total_bounds
        out_gdf = self.get_antenna_location_from_bb(xmin, ymin, xmax, ymax)
        return out_gdf

    def get_antenna_location_from_bb(self, xmin, ymin, xmax, ymax):
        return self.antenna_gdf_epsg4326.cx[xmin:xmax, ymin:ymax].reset_index(drop=True)
