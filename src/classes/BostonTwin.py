from .BostonModel import BostonModel
from .BostonAntennas import BostonAntennas
from typing import Union
from pathlib import Path
import geopandas as gpd
import mitsuba as mi
import numpy as np
import pyproj
import time

from sionna.rt import load_scene, Transmitter, Receiver
from src.utils.geo_utils import gdf2localcrs, plot_geodf


class BostonTwin:
    def __init__(
        self,
        dataset_dir: Union[Path, str] = Path("dataset"),
    ):
        self.boston_model_path = dataset_dir.joinpath("boston3d")
        self.boston_model = BostonModel(self.boston_model_path)
        self.boston_antennas = BostonAntennas(dataset_dir.joinpath("boston_antennas"))

        self.current_scene_name = ""
        self.current_scene_gdf = None
        self.current_sionna_scene = None
        self.current_mi_scene = None
        self.current_scene_antennas = None

        self.node_height = 10
        self.txs = []
        self.rxs = []

    def _check_scene(self):
        if self.current_scene_name is None:
            raise ValueError("Scene not set! Run the set_scene(<scene_name>) method specifying the scene name.")

    def get_scene_names(self):
        return self.boston_model.tile_names

    def set_scene(self, scene_name:str):
        self.current_scene_name = scene_name

        self.tile_info_path = self.boston_model.tiles_dict[self.current_scene_name][
            "tileinfo_path"
        ]

        self.current_scene_info_gdf = gpd.GeoDataFrame.from_file(self.tile_info_path)
        self.current_scene_center = (
            self.current_scene_info_gdf["Centr_X_m"].values[0],
            self.current_scene_info_gdf["Centr_Y_m"].values[0],
        )

        self.mi_scene_path = self.boston_model.tiles_dict[self.current_scene_name][
            "mi_scene_path"
        ]

    def load_mi_scene(self):
        self._check_scene()
        self.current_mi_scene = mi.load_file(str(self.mi_scene_path.resolve()))

    def load_scene_geodf(self):
        self._check_scene()
        self.geo_scene_path = self.boston_model.tiles_dict[self.current_scene_name][
            "geo_scene_path"
        ]
        scene_gdf = gpd.GeoDataFrame.from_file(self.geo_scene_path)
        self.current_scene_gdf_lonlat = scene_gdf

        scene_gdf_localcrs = gdf2localcrs(scene_gdf)
        self.current_scene_gdf = self.translate_gdf(
            scene_gdf_localcrs,
            xoff=-self.current_scene_center[0],
            yoff=-self.current_scene_center[1],
        )

    def load_antennas(self):
        self._check_scene()
        self.scene_antennas_gdf_lonlat = self.boston_antennas.get_antenna_location_from_gdf(
            self.current_scene_info_gdf
        )
        self.current_scene_antennas = self.translate_gdf(
            self.scene_antennas_gdf_lonlat,
            xoff=-self.current_scene_center[0],
            yoff=-self.current_scene_center[1],
        )

    def load_bostontwin(
        self, scene_name, load_sionna=True, load_mi_scene=False, load_geodf=False
    ):
        self.set_scene(scene_name)

        self.load_antennas()

        if load_sionna:
            self.current_sionna_scene = load_scene(str(self.mi_scene_path))

        if load_mi_scene:
            self.load_mi_scene()

        if load_geodf:
            self.current_scene_gdf = self.load_scene_geodf(scene_name)

        return self.current_sionna_scene, self.current_scene_antennas

    def get_mi_scene(self):
        return self.current_mi_scene

    def plot_buildings(self, basemap:bool=False, **plot_kwargs):
        self.load_scene_geodf()

        ax = plot_geodf(
            self.current_scene_gdf_lonlat,
            basemap=basemap,
            title=self.current_scene_name,
            **plot_kwargs,
        )
        return ax

    def plot_antennas(self, basemap:bool=False, **plot_kwargs):
        return plot_geodf(
            self.scene_antennas_gdf_lonlat,
            basemap=basemap,
            title=self.current_scene_name,
            **plot_kwargs,
        )

    def plot_twin(self, basemap:bool=False):
        ax = self.plot_buildings(basemap=basemap)
        ax = self.plot_antennas(basemap=False, ax=ax)
        return ax

    def add_scene_antennas(
        self,
        tx_antenna_ids,
        rx_antenna_ids,
        tx_names=[],
        rx_names=[],
        tx_params=[],
        rx_params=[],
    ):
        if not tx_names:
            tx_names = [f"TX_{i}" for i in range(len(tx_antenna_ids))]

        for tx_idx, (tx_name, tx_antenna_idx) in enumerate(
            zip(tx_names, tx_antenna_ids)
        ):
            antenna_coords = list(
                self.current_scene_antennas.loc[tx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            if len(tx_params) > 0:
                tx_par = tx_params[tx_idx]
            else:
                tx_par = {}
            tx = Transmitter(tx_name, position=antenna_coords, **tx_par)
            self.current_sionna_scene.add(tx)
            self.txs.append(tx)

        if not rx_names:
            rx_names = [f"RX_{i}" for i in range(len(rx_antenna_ids))]

        for rx_idx, (rx_name, rx_antenna_idx) in enumerate(
            zip(rx_names, rx_antenna_ids)
        ):
            antenna_coords = list(
                self.current_scene_antennas.loc[rx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            if len(rx_params) > 0:
                rx_par = rx_params[rx_idx]
            else:
                rx_par = {}
            rx = Receiver(rx_name, position=antenna_coords, **rx_par)
            self.current_sionna_scene.add(rx)
            self.rxs.append(rx)

        return dict(zip(tx_names + rx_names, self.txs + self.rxs))

    def generate_scene_from_radius(
        self, scene_name, center_lon, center_lat, side_m, load=False
    ):
        radius = np.sqrt(2) * side_m  # m
        azimuths = [45, 225]

        geod = pyproj.Geod(ellps="WGS84")
        lon1, lat1, _ = geod.fwd(center_lon, center_lat, azimuths[0], radius)
        lon2, lat2, _ = geod.fwd(center_lon, center_lat, azimuths[1], radius)
        bbox = [lon1, lat1, lon2, lat2]
        print("Selecting models within the area...")
        t0 = time.time()
        boston_gdf = gpd.GeoDataFrame.from_file(
            self.boston_model_path.joinpath("boston.geojson"), bbox=bbox
        )
        t1 = time.time()
        print(f"Done. ({t1-t0:.2f} s)")
        self.boston_model.generate_scene_from_model_gdf(
            boston_gdf, (center_lon, center_lat), scene_name
        )

    @staticmethod
    def translate_gdf(gdf, xoff, yoff):
        gdf["geometry"] = gdf.translate(xoff=xoff, yoff=yoff)
        return gdf
