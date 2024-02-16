import mitsuba as mi
from .constants import FT2M_FACTOR
import pyproj
from pathlib import Path

data_dir = Path(__file__).parents[2].joinpath("data")

def obj2ply_mi(
    x_shift,
    y_shift,
    z_shift,
    scale_factor,
    obj_path,
    ply_path,
    ply_path_relative,
):
    ## Read the OBJ file and convert it to PLY, changing the unit to meters and centering the model file
    # the PLY file is saved in relative coordinates, centered in [0,0]
    # note that the unit is converted from feet to meters
    to_relative = mi.ScalarTransform4f.scale(scale_factor).translate(
        [
            x_shift,
            y_shift,
            z_shift,
        ]
    )
    in_model_dict = {
        "type": "obj",
        "filename": str(obj_path.resolve()),
        "face_normals": True,
        "to_world": to_relative,
    }
    model = mi.load_dict(in_model_dict)  # load the model into mitsuba
    n_tri = model.face_count()

    # write the scaled and centered model to PLY file
    model.write_ply(str(ply_path.resolve()))

    # create the new mitsuba dict for the PLY file
    center_x_m_from_info = x_shift * scale_factor
    center_y_m_from_info = y_shift * scale_factor
    z_min_m_from_info = z_shift * scale_factor
    to_world = mi.ScalarTransform4f.translate(
        [
            -center_x_m_from_info,
            -center_y_m_from_info,
            -z_min_m_from_info,
        ]
    )
    out_model_dict = {
        "to_world": to_world,
        "type": "ply",
        "filename": str(ply_path.relative_to(ply_path_relative)),
        "face_normals": True,
    }
    return out_model_dict, n_tri


def gdf2localcrs(in_coords_gdf):
    # # translate to the BDPA reference system
    # origin_coords_epsg4326 = Point(LOCAL_CRS_ORIGIN_GEO)
    # origin_gdf_epsg2249 = gpd.GeoDataFrame(
    #     index=[0], geometry=[origin_epsg4326], crs="epsg:4326"
    # ).to_crs("epsg:2249")
    # local_offset = list(origin_gdf_epsg2249["geometry"][0].coords[0])
    # local_offset = LOCAL_CRS_ORIGIN_STATE

    # in_coords_gdf_epsg2249 = in_coords_gdf.to_crs("epsg:2249")
    # in_coords_gdf_local_crs = in_coords_gdf_epsg2249.translate(
    #     xoff=-local_offset[0], yoff=-local_offset[1]
    # )

    prj_path = data_dir.joinpath("Metro_Boston_3D_CRS.prj")
    with open(prj_path, "r") as f:
        prj_str = f.readline()
        crs = pyproj.CRS.from_wkt(prj_str)

    # convert to meters
    out_coords_gdf_local_crs = in_coords_gdf.to_crs(crs)
    out_coords_gdf_local_crs["geometry"] = out_coords_gdf_local_crs.scale(
        FT2M_FACTOR, FT2M_FACTOR, origin=(0, 0)
    )
    return out_coords_gdf_local_crs
