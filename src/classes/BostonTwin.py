from .BostonModel import BostonModel
from .BostonAntennas import BostonAntennas
from typing import Union
from pathlib import Path
import geopandas as gpd
import mitsuba as mi

from sionna.rt import load_scene, Transmitter, Receiver
from src.functions.geo_utils import gdf2localcrs


class BostonTwin:
    def __init__(
        self,
        dataset_dir: Union[Path, str] = Path("dataset"),
    ):
        self.boston_model = BostonModel(dataset_dir.joinpath("boston3d"))
        self.boston_antennas = BostonAntennas(dataset_dir.joinpath("boston_antennas"))

        self.current_scene_name = ""
        self.current_scene_gdf = None
        self.current_sionna_scene = None
        self.current_mi_scene = None
        self.current_scene_antennas = None

        self.node_height = 10
        self.txs = []
        self.rxs = []

    def get_scene_names(self):
        return self.boston_model.tile_names

    def set_scene(self, scene_name):
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
        self.current_mi_scene = mi.load_file(str(self.mi_scene_path.resolve()))

    def load_scene_geodf(self):
        self.geo_scene_path = self.boston_model.tiles_dict[self.current_scene_name][
            "geo_scene_path"
        ]
        scene_gdf = gpd.GeoDataFrame.from_file(self.geo_scene_path)
        scene_gdf = gdf2localcrs(scene_gdf)
        self.current_scene_gdf = self.translate_gdf(
            scene_gdf,
            xoff=-self.current_scene_center[0],
            yoff=-self.current_scene_center[1],
        )

    def load_antennas(self):
        scene_antennas_gdf = self.boston_antennas.get_antenna_location_from_gdf(
            self.current_scene_info_gdf
        )
        self.current_scene_antennas = self.translate_gdf(
            scene_antennas_gdf,
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

    def add_scene_antennas(
        self,
        tx_antenna_ids,
        rx_antenna_ids,
        tx_names=[],
        rx_names=[],
        tx_params=[{}],
        rx_params=[{}],
    ):
        if not tx_names:
            tx_names = [f"TX_{i}" for i in range(len(tx_antenna_ids))]

        
        for tx_idx, (tx_name, tx_antenna_idx) in enumerate(zip(tx_names, tx_antenna_ids)):
            antenna_coords = list(
                self.current_scene_antennas.loc[tx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            tx_par = tx_params[tx_idx]
            tx = Transmitter(tx_name, position=antenna_coords, **tx_par)
            self.current_sionna_scene.add(tx)
            self.txs.append(tx)

        if not rx_names:
            rx_names = [f"RX_{i}" for i in range(len(rx_antenna_ids))]

        for rx_idx, (rx_name, rx_antenna_idx) in enumerate(zip(rx_names, rx_antenna_ids)):
            antenna_coords = list(
                self.current_scene_antennas.loc[rx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            print(antenna_coords)
            rx_par = rx_params[rx_idx]
            rx = Receiver(rx_name, position=antenna_coords, **rx_par)
            self.current_sionna_scene.add(rx)
            self.rxs.append(rx)

        return dict(zip(tx_names + rx_names, self.txs + self.rxs))

    @staticmethod
    def translate_gdf(gdf, xoff, yoff):
        gdf["geometry"] = gdf.translate(xoff=xoff, yoff=yoff)
        return gdf
