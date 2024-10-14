from pathlib import Path

import mitsuba as mi
import numpy as np
import open3d as o3d

from bostontwin.utils.constants import FT2M_FACTOR


def dir_obj2ply(obj_dir, ply_dir, recursive=True, **kwargs):
    obj_list = list(obj_dir.glob("**/*.obj")) if recursive else list(obj_dir.glob("*.obj"))
    for obj_path in obj_list:
        ply_path = ply_dir.joinpath(obj_path.stem + ".ply")
        if not ply_path.is_file():
            obj2ply(obj_path, ply_path, **kwargs)

def get_mi_dict(
    model_material,
    model_center_x,
    model_center_y,
    model_z,
    ply_path,
    ply_path_relative,
):
    # create the new mitsuba dict for the PLY file
    to_world = mi.ScalarTransform4f.translate(
        [
            model_center_x,
            model_center_y,
            model_z
        ]
    )
    out_model_dict = {
        "to_world": to_world,
        "type": "ply",
        "filename": str(ply_path.relative_to(ply_path_relative)),
        "face_normals": True,
        "bsdf": {"type": "ref", "id": model_material},
    }
    return out_model_dict


def get_model_center(model_path):
    mesh = o3d.io.read_triangle_mesh(str(model_path))
    return mesh.get_center()


def render_model(scene_dict, model_dict, model_name, sensor_name, center_view=True):
    scene_dict[model_name] = model_dict

    if center_view:
        center = get_model_center(model_dict["filename"])
        scene_dict[sensor_name]["to_world"] = mi.Transform4f.look_at(
            origin=[center[0] - 10, center[1] - 10, 30],  # camera position
            target=center,  # look at
            up=[0, 0, 1],  # up vector
        )

    scene = mi.load_dict(scene_dict)
    render = mi.render(scene)
    return render


def obj2ply(obj_path, ply_path, ft2m=True, center=True):
    ## Convert the OBJ file to PLY, changing the unit to meters and centering the model file
    # the PLY file is saved in relative coordinates, centered in [0,0]
    # note that the unit is converted from feet to meters

    # check if input file exists
    obj_path = Path(obj_path)
    if not obj_path.is_file():
        raise FileNotFoundError(f"File not found: {obj_path}")

    # check if output directory exists
    ply_path = Path(ply_path)
    if not ply_path.parent.is_dir():
        raise FileNotFoundError(f"Directory not found: {ply_path.parent}")

    # load the obj file
    obj_path = obj_path.resolve()
    mesh = o3d.io.read_triangle_mesh(str(obj_path))

    # if ft2m is True, convert the mesh from feet to meters
    if ft2m:
        mesh.scale(FT2M_FACTOR, center=mesh.get_center())

    # if center is True, center the mesh in (0,0,0)
    mesh_trans = [0, 0, 0]
    if center:
        mesh_center = mesh.get_center()
        mesh_min_z = np.asarray(mesh.vertices)[:, 2].min()
        mesh_trans = [-mesh_center[0], -mesh_center[1], -mesh_min_z]
        mesh.translate(mesh_trans)

    # compute the normals
    mesh.compute_vertex_normals()
    mesh.compute_triangle_normals()

    # save the mesh as a ply file
    o3d.io.write_triangle_mesh(str(ply_path), mesh)

    mesh_n_tri = len(mesh.triangles)

    return mesh_trans, mesh_n_tri


def create_ground_dict(
    model_material,
    x_shift,
    y_shift,
    z_shift,
    scale_factor,
    base_rect_path,
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
        "type": "ply",
        "filename": str(base_rect_path.resolve()),
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
        "bsdf": {"type": "ref", "id": model_material},
    }
    return out_model_dict, n_tri
