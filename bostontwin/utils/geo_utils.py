from functools import partial
from pathlib import Path
from typing import List, Tuple, Union

import contextily as cx
import geopandas as gpd
import pyproj
from pyproj.crs import CompoundCRS, GeographicCRS, ProjectedCRS, VerticalCRS
from pyproj.crs.coordinate_operation import TransverseMercatorConversion
from pyproj.crs.coordinate_system import Cartesian2DCS, VerticalCS
from shapely.geometry import Point, box
from shapely.ops import transform

data_dir = Path(__file__).parents[2].joinpath("data")

# def gdf2localcrs(in_coords_gdf):
#     # # translate to the BDPA reference system
#     # origin_coords_epsg4326 = Point(LOCAL_CRS_ORIGIN_GEO)
#     # origin_gdf_epsg2249 = gpd.GeoDataFrame(
#     #     index=[0], geometry=[origin_epsg4326], crs="epsg:4326"
#     # ).to_crs("epsg:2249")
#     # local_offset = list(origin_gdf_epsg2249["geometry"][0].coords[0])
#     # local_offset = LOCAL_CRS_ORIGIN_STATE

#     # in_coords_gdf_epsg2249 = in_coords_gdf.to_crs("epsg:2249")
#     # in_coords_gdf_local_crs = in_coords_gdf_epsg2249.translate(
#     #     xoff=-local_offset[0], yoff=-local_offset[1]
#     # )

#     prj_path = data_dir.joinpath("Metro_Boston_3D_CRS.prj")
#     with open(prj_path, "r") as f:
#         prj_str = f.readline()
#         crs = pyproj.CRS.from_wkt(prj_str)

#     # convert to meters
#     out_coords_gdf_local_crs = in_coords_gdf.to_crs(crs)
#     out_coords_gdf_local_crs["geometry"] = out_coords_gdf_local_crs.scale(
#         FT2M_FACTOR, FT2M_FACTOR, FT2M_FACTOR, origin=(0, 0)
#     )
#     return out_coords_gdf_local_crs

def gdf2crs(in_gdf:gpd.GeoDataFrame, crsTransformer:pyproj.Transformer)->gpd.GeoDataFrame:
    transformer = partial(transform, crsTransformer.transform)
    out_gdf = in_gdf.set_geometry(
        in_gdf.geometry.apply(transformer), inplace=False
    )
    out_gdf.crs = crsTransformer.target_crs
    return out_gdf

def plot_geodf(
    geodf, basemap: bool = False, title: str = "", **plot_kwargs
):
    ax = geodf.plot(**plot_kwargs)
    if basemap:
        cx.add_basemap(
            ax,
            crs=geodf.crs,
            source=cx.providers.CartoDB.VoyagerNoLabels,
            zoom=17,
        )
        cx.add_basemap(
            ax,
            crs=geodf.crs,
            source=cx.providers.CartoDB.VoyagerOnlyLabels,
            zoom=16,
        )
    if geodf.crs.is_geographic:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    else:
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")

    if title:
        ax.set_title(title)

    return ax

def get_crs(scene_name:str, scene_center_lon_lat:Tuple):
    center_lon, center_lat = scene_center_lon_lat
    vertcrs = VerticalCRS(
        name="NAVD88 height",
        datum="North American Vertical Datum 1988",
        vertical_cs=VerticalCS(),
        geoid_model="GEOID12B",
    )
    projcrs = ProjectedCRS(
        conversion=TransverseMercatorConversion(
            latitude_natural_origin=center_lat,
            longitude_natural_origin=center_lon,
            false_easting=0,
            false_northing=0,
            scale_factor_natural_origin=1,
        ),
        geodetic_crs=GeographicCRS(datum="North American Datum 1983"),
        cartesian_cs=Cartesian2DCS(),
    )
    compcrs = CompoundCRS(
        name=scene_name, components=[projcrs, vertcrs]
    )
    return compcrs

def check_point_in_area_of_use(crs, in_location):
    # assert in_crs.is_geographic, "Input CRS must be geographic"
    # get out_crs bounds
    bounding_box = box(*crs.area_of_use.bounds)
    out_polygon = Point(*in_location)
    return bounding_box.contains(out_polygon)

def check_area_of_use(in_crs:pyproj.CRS, out_crs:pyproj.CRS, in_coords:Union[Tuple,List])->bool:
    # get out_crs bounds
    if len(in_coords)>2:
        for point_coord in in_coords:
            in_location_out_crs = pyproj.Transformer.from_proj(
                in_crs, out_crs, always_xy=True
            ).transform(*point_coord)
            return check_point_in_area_of_use(out_crs, in_location_out_crs) & check_point_in_area_of_use(in_crs, point_coord)
    else:
        in_location_out_crs = pyproj.Transformer.from_proj(
            in_crs, out_crs, always_xy=True
        ).transform(*in_coords)
        return check_point_in_area_of_use(
            in_crs, out_crs, in_coords
        ) & check_point_in_area_of_use(out_crs, in_location_out_crs)
