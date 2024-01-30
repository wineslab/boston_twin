from .BostonModel import BostonModel
from .BostonAntennas import BostonAntennas
from typing import Union
from pathlib import Path
import geopandas as gpd

from sionna.rt import load_scene, Transmitter, Receiver, PlanarArray, Camera


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

    def load_bostontwin(self, scene_name):
        self.tile_info_path = self.boston_model.tiles_dict[scene_name]["tileinfo_path"]
        self.geo_scene_path = self.boston_model.tiles_dict[scene_name]["geo_scene_path"]
        self.mi_scene_path = self.boston_model.tiles_dict[scene_name]["mi_scene_path"]

        self.current_scene_name = scene_name
        self.current_scene_info_gdf = gpd.GeoDataFrame.from_file(self.tile_info_path)
        self.current_scene_center = (
            self.current_scene_info_gdf["Centr_X_m"].values[0],
            self.current_scene_info_gdf["Centr_Y_m"].values[0],
        )

        self.current_sionna_scene = load_scene(str(self.mi_scene_path))
        self.current_mi_scene = self.current_sionna_scene.mi_scene

        self.current_scene_antennas = (
            self.boston_antennas.get_antenna_location_from_gdf(
                self.current_scene_info_gdf
            )
        )
        self.current_scene_antennas["geometry"] = self.current_scene_antennas.translate(
            xoff=-self.current_scene_center[0], yoff=-self.current_scene_center[1]
        )
        return self.current_sionna_scene, self.current_scene_antennas

    def get_mi_scene(self):
        return self.current_mi_scene

    def add_scene_antennas(
        self, tx_antenna_ids, rx_antenna_ids, tx_names=[], rx_names=[]
    ):
        if not tx_names:
            tx_names = [f"TX_{i}" for i in range(len(tx_antenna_ids))]

        for tx_name, tx_antenna_idx in zip(tx_names, tx_antenna_ids):
            antenna_coords = list(
                self.current_scene_antennas.loc[tx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            tx = Transmitter(tx_name, position=antenna_coords)
            self.current_sionna_scene.add(tx)
            self.txs.append(tx)

        if not rx_names:
            rx_names = [f"RX_{i}" for i in range(len(rx_antenna_ids))]

        for rx_name, rx_antenna_idx in zip(rx_names, rx_antenna_ids):
            antenna_coords = list(
                self.current_scene_antennas.loc[rx_antenna_idx, "geometry"].coords[0]
            )
            antenna_coords.append(self.node_height)
            print(antenna_coords)
            rx = Receiver(rx_name, position=antenna_coords)
            self.current_sionna_scene.add(rx)
            self.rxs.append(rx)

        return dict(zip(tx_names + rx_names, self.txs + self.rxs))
