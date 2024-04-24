import json
import time
from pathlib import Path
from typing import List, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import mitsuba as mi
import numpy as np
import pyproj
from matplotlib.axes import Axes
from sionna.rt import Receiver, Transmitter, load_scene

from src.utils.geo_utils import gdf2localcrs, plot_geodf

from .BostonAntennas import BostonAntennas
from .BostonModel import BostonModel


class BostonTwin:
    """BostonTwin Class

    The BostonTwin class implements BostonTwin, the Boston Digital Twin for wireless communications.
    The class contains two main variables:
        - `boston_model`, that includes a number of methods to access the 3D model of the structures in Boston (buildings, bridges, walls, etc..)
        - `boston_antennas`, that includes the georeferenced locations to the wireless antennas in Boston

    Attributes
    ----------
    dataset_dir : Path
        Path to the Boston Twin.
    boston_model_path : Path
        Path to the Boston Model folder (dataset_dir/boston3d).
    boston_model : BostonModel
        BostonModel instance, containing the information on the 3D model of the structures in Boston.
    boston_antennas : BostonAntennas
        BostonAntennas instance, that includes the georeferenced locations to the wireless antennas in Boston.
    current_scene_name : str
        Name of the current scene.
    current_scene_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the information on the structures of the current scene.
    current_sionna_scene : sionna.rt.Scene
        Sionna Scene instance of the current scene.
    current_mi_scene : mitsuba.Scene
        Mitsuba Scene instance of the current scene.
    current_scene_antennas_localcrs : geopandas.GeoDataFrame
        GeoDataFrame containing the information on the antennas of the current scene.
    current_scene_txrx_localcrs : geopandas.GeoDataFrame
        GeoDataFrame containing the information on the transmitters and receivers of the current scene.
    node_height : float
        Height of the antennas of the scene. For now, all the antennas have the same height.
    txs: list
        List of transmitters in the current Sionna scene.
    rxs: list
        List of receivers in the current Sionna scene.
    _node_pos_dict: dict
        Dictionary containing the position (in the local CRS) of the nodes in the current scene.
        The keys correspond to the node index in the current_scene_antennas_localcrs GeoDataFrame.
        The values are [x,y] vectors (local CRS).

    Methods
    -------
    get_scene_names()
        Return the list of scene names currently present in BostonTwin.
    load_boston_twin(scene_name, load_sionna=True, load_mi_scene=False, load_geodf=False)
        Load `scene_name` as the current scene. Choose among those listed by `get_scene_names()`,
        or add a new one with the `generat_scene_from_radius()` method.
    plot_buildings()
        Plot the 2D footprint of the buildings in the current scene.
    plot_antennas()
        Plot the location of the wirelessa antennas in the current scene.
    plot_twin()
        Combine `plot_antennas()` and `plot_buildings()`
    add_scene_antennas(tx_antenna_ids, rx_antenna_ids, tx_names=[], rx_names=[], tx_params=[], rx_params=[])
        Add the specified Transmitters and Receivers to the Sionna scene.
    generate_scene_from_radius()
        Generate a new scene specifying its center and radius.
    """

    def __init__(
        self,
        dataset_dir: Union[Path, str] = Path("dataset"),
    ):
        if isinstance(dataset_dir, str):
            dataset_dir = Path(dataset_dir)
        self.dataset_dir = dataset_dir
        self.boston_model_path = dataset_dir.joinpath("boston3d")
        self.boston_model = BostonModel(self.boston_model_path)
        self.boston_antennas = BostonAntennas(dataset_dir.joinpath("boston_antennas"))

        self.current_scene_name = ""
        self.current_scene_gdf = None
        self.current_sionna_scene = None
        self.current_mi_scene = None
        self.current_scene_antennas_localcrs = None
        self.current_scene_txrx_localcrs = None

        self.node_height = 10
        self.txs = []
        self.rxs = []
        self._node_pos_dict = {}

    def _check_scene(self):
        if self.current_scene_name is None:
            raise ValueError(
                "Scene not set! Run the set_scene(<scene_name>) method specifying the scene name."
            )

    def _generate_node_pos_dict(self):
        self._node_pos_dict = dict(
            zip(
                self.current_scene_antennas_localcrs.index,
                [
                    {"x": c.coords[0][0], "y": c.coords[0][1]}
                    for c in self.current_scene_antennas_localcrs.geometry
                ],
            )
        )

    def get_scene_names(self) -> List[str]:
        """Return the list of scene names currently present in BostonTwin. The files describing are found in the `self.dataset_dir` directory.

        Returns
        -------
        list
            List of the names of scenes available in BostonTwin.
        """
        return self.boston_model.tile_names

    def set_scene(self, scene_name: str):
        """Set the current scene to `scene_name`.

        Parameters
        ----------
        scene_name : str
            Name of the scene to be set. Must be among those returned by `get_scene_names()`.
        """
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
        self.scene_antennas_gdf_lonlat = (
            self.boston_antennas.get_antenna_location_from_gdf(
                self.current_scene_info_gdf
            )
        )
        scene_antennas_gdf_localcrs = gdf2localcrs(self.scene_antennas_gdf_lonlat)
        self.current_scene_antennas_localcrs = self.translate_gdf(
            scene_antennas_gdf_localcrs,
            xoff=-self.current_scene_center[0],
            yoff=-self.current_scene_center[1],
        )

    def load_bostontwin(
        self, scene_name: str, load_sionna=True, load_mi_scene=False, load_geodf=False
    ):
        """Load `scene_name` as the current scene.

        Parameters
        ----------
        scene_name : str
            Name of the scene to load. The scene files must be present in the `self.dataset_dir` directory.
        load_sionna : bool, optional
            Load the sionna scene. Defaults to True.
        load_mi_scene : bool, optional
            Load the Mitsuba scene. Defaults to False.
        load_geodf : bool, optional
            Load the GeoDataframe of the scene. Defaults to False.

        Returns
        -------
        current_sionna_scene : sionna.rt.scene
            The sionna.rt.Scene representing the current scene.
        current_scene_antennas : gpd.GeoDataFrame
            A Geopandas GeoDataFrame containing the information and location of the antennas present in the current scene.
        """
        self.set_scene(scene_name)

        self.load_antennas()

        self._generate_node_pos_dict()

        if load_sionna:
            self.current_sionna_scene = load_scene(str(self.mi_scene_path))

        if load_mi_scene:
            self.load_mi_scene()

        if load_geodf:
            self.current_scene_gdf = self.load_scene_geodf(scene_name)

        return self.current_sionna_scene, self.current_scene_antennas_localcrs

    def get_mi_scene(self):
        return self.current_mi_scene

    def plot_buildings(
        self, basemap: bool = False, local_crs: bool = False, **plot_kwargs
    ) -> Axes:
        """Plot the buildings 2D footprint.

        Parameters
        ----------
        basemap : bool, optional
            Add a map as background. Defaults to False.
        local_crs : bool, optional
            Use the local Coordinate Reference System. Incompatible with basemap. Defaults to False.

        Returns
        -------
        ax : Axes
            Building footprint plot.
        """

        if basemap and local_crs:
            raise ValueError(
                "'basemap' and 'local_crs' are currently incompatible. Please choose one."
            )

        self.load_scene_geodf()

        if local_crs:
            plot_gdf = self.current_scene_gdf
        else:
            plot_gdf = self.current_scene_gdf_lonlat

        ax = plot_geodf(
            plot_gdf,
            basemap=basemap,
            title=self.current_scene_name,
            **plot_kwargs,
        )
        return ax

    def plot_antennas(
        self,
        basemap: bool = False,
        local_crs: bool = False,
        annotate: bool = False,
        **plot_kwargs,
    ) -> Axes:
        """Plot the location of the antennas in the current scene.

        Parameters
        ----------
        basemap : bool, optional
            Add a map as background. Defaults to False.
        local_crs : bool, optional
            Use the local Coordinate Reference System. Incompatible with basemap. Defaults to False.
        annotate : bool, optional
            Add annotation to map. Defaults to False.

        Returns
        -------
        ax : Axes
            Antenna location plot.
        """
        if basemap and local_crs:
            raise ValueError(
                "'basemap' and 'local_crs' are currently incompatible. Please choose one."
            )

        if local_crs:
            plot_gdf = self.current_scene_antennas_localcrs
        else:
            plot_gdf = self.scene_antennas_gdf_lonlat

        ax = plot_geodf(
            plot_gdf,
            basemap=basemap,
            title=self.current_scene_name,
            **plot_kwargs,
        )
        if annotate:
            for idx, row in plot_gdf.iterrows():
                xy = row.geometry.centroid.coords[0]
                format_coords = [f"{x:.1f}" for x in row.geometry.centroid.coords[0]]
                ax.annotate(
                    text=format_coords,
                    xy=xy,
                    horizontalalignment="center",
                    verticalalignment="top",
                    size=4,
                    bbox=dict(
                        facecolor=[1, 1, 1, 0.3],
                        edgecolor="black",
                        boxstyle="round,pad=2",
                    ),
                )
                if self.current_scene_txrx_localcrs is not None:
                    pole_ID_lonlat = row["Pole_Identifying_Number"]
                    pole_ID_localcrs = self.current_scene_txrx_localcrs[
                        "Pole_Identifying_Number"
                    ]
                    match_id = pole_ID_localcrs == pole_ID_lonlat
                    name_text = self.current_scene_txrx_localcrs.loc[
                        match_id, "Name"
                    ].values[0]
                else:
                    name_text = f"idx-{idx}"

                ax.annotate(
                    text=name_text,
                    xy=xy,
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    size=8,
                )

        return ax

    def plot_twin(self, basemap: bool = False, local_crs : bool = False, annotate : bool = False, ax : plt.Axes = None) -> Axes:
        """Plot the location of the antennas and the building footprint.

        Parameters
        ----------
        basemap : bool, optional
            Add a map as background. Defaults to False.
        local_crs : bool, optional
            Use the local Coordinate Reference System. Incompatible with basemap. Defaults to False.
        annotate : bool, optional
            Add antenna location annotation to map. Defaults to False.

        Returns
        -------
        ax : Axes
            Plot of the antennas among the building 2D footprint.
        """

        if basemap and local_crs:
            raise ValueError(
                "'basemap' and 'local_crs' are currently incompatible. Please choose one."
            )

        ax = self.plot_buildings(basemap=basemap, color="k", local_crs=local_crs, ax=ax)
        ax = self.plot_antennas(basemap=False, ax=ax, color="r", local_crs=local_crs, annotate=annotate)
        return ax

    def add_scene_antennas(
        self,
        tx_antenna_ids: Union[List[int], np.typing.ArrayLike],
        rx_antenna_ids: Union[List[int], np.typing.ArrayLike],
        tx_names: List[str] = [],
        rx_names: List[str] = [],
        tx_params: list = [],
        rx_params: list = [],
    ) -> dict:
        """Add antennas to the current scene.

        Parameters
        ----------
        tx_antenna_ids : List[int]
            List of indices of the antennas to be used as Transmitters.
        rx_antenna_ids : List[int]
            List of IDs to be used as Receiver names (RX_{id}).
        tx_names : List[str], optional
            List of Transmitter names (Default: TX_{id}). Defaults to [].
        rx_names : List[str], optional
            List of Receiver names (Default: RX_{id}). Defaults to [].
        tx_params : list, optional
            Parameters for the Transmitters. Refer to the Sionna Ray Tracer documentation. Defaults to [].
        rx_params : list, optional
            Parameters for the Receivers. Refer to the Sionna Ray Tracer documentation. Defaults to [].

        Returns
        -------
        nodes_dict : dict
            Dictionary with the names of the Transmitters/Receivers as keys and the corresponding Sionna object as values.
        """
        txrx_ids = np.concatenate([tx_antenna_ids, rx_antenna_ids])
        self.current_scene_txrx_localcrs = self.current_scene_antennas_localcrs.loc[
            txrx_ids, :
        ]
        self.current_scene_txrx_localcrs.loc[tx_antenna_ids, "TX/RX"] = "TX"
        self.current_scene_txrx_localcrs.loc[rx_antenna_ids, "TX/RX"] = "RX"

        if not tx_names:
            tx_names = [f"TX_{i}" for i in range(len(tx_antenna_ids))]
        self.current_scene_txrx_localcrs.loc[tx_antenna_ids, "Name"] = tx_names

        for tx_idx, (tx_name, tx_antenna_idx) in enumerate(
            zip(tx_names, tx_antenna_ids)
        ):
            antenna_coords = list(
                self.current_scene_antennas_localcrs.loc[
                    tx_antenna_idx, "geometry"
                ].coords[0]
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
        self.current_scene_txrx_localcrs.loc[rx_antenna_ids, "Name"] = rx_names

        for rx_idx, (rx_name, rx_antenna_idx) in enumerate(
            zip(rx_names, rx_antenna_ids)
        ):
            antenna_coords = list(
                self.current_scene_antennas_localcrs.loc[
                    rx_antenna_idx, "geometry"
                ].coords[0]
            )
            antenna_coords.append(self.node_height)
            if len(rx_params) > 0:
                rx_par = rx_params[rx_idx]
            else:
                rx_par = {}
            rx = Receiver(rx_name, position=antenna_coords, **rx_par)
            self.current_sionna_scene.add(rx)
            self.rxs.append(rx)

        nodes_dict = dict(zip(tx_names + rx_names, self.txs + self.rxs))
        return nodes_dict

    def generate_scene_from_radius(
        self,
        scene_name: str,
        center_lon: float,
        center_lat: float,
        side_m: float,
        load=False,
    ):
        """Generate a new scene specifying its center and radius.

        Parameters
        ----------
        scene_name : str
            Name of the new scene.
        center_lon : float
            Longitude of the center.
        center_lat : float
            Latitude of the center.
        side_m : float
            Radius of the scene.
        load : bool, optional
            Load the scene as current scene. Defaults to False.
        """
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

        if load:
            self.set_scene(scene_name)

    def export_scene_antennas(self, out_path: Union[Path, str]):
        """Export to file the location (in the local CRS) of the antennas in the current scene.

        Parameters
        ----------
        out_path : Union[Path,str]
            Path where to export the antenna location.
        """
        if not isinstance(out_path, Path):
            out_path = Path(out_path)

        if out_path.suffix.lower() != ".json":
            raise ValueError("Currently, only JSON is supported as export format.")

        if out_path.is_file():
            raise FileExistsError(f"File {out_path} already exists!")

        if not out_path.parent.is_dir():
            raise FileNotFoundError(f"Directory {out_path.parent} not found!")

        with open(out_path, "w") as f:
            json.dump(self._node_pos_dict, f, indent=4)

    # Static Methods
    @staticmethod
    def translate_gdf(gdf, xoff, yoff):
        gdf["geometry"] = gdf.translate(xoff=xoff, yoff=yoff)
        return gdf
