from pathlib import Path

import mitsuba as mi
import numpy as np
import open3d as o3d
from typing import Union
from src.bostontwin.utils.constants import FT2M_FACTOR
from src.bostontwin.utils.mi_utils import get_transformation_matrix


def read_mesh(mesh_path: Union[Path, str]) -> o3d.geometry.TriangleMesh:
    mesh_path = Path(mesh_path)
    if not mesh_path.is_file():
        raise FileNotFoundError(f"File not found: {mesh_path}")
    mesh = o3d.io.read_triangle_mesh(str(mesh_path))
    return mesh

def merge_meshes(meshes: list) -> o3d.geometry.TriangleMesh:
    merged_mesh = o3d.geometry.TriangleMesh()
    for mesh in meshes:
        merged_mesh += mesh
    return merged_mesh

def save_mesh(mesh: o3d.geometry.TriangleMesh, mesh_path: Union[Path, str]):
    mesh_path = Path(mesh_path)
    if not mesh_path.parent.is_dir():
        raise FileNotFoundError(f"Directory not found: {mesh_path.parent}")
    o3d.io.write_triangle_mesh(str(mesh_path), mesh)

def dir_obj2ply(obj_dir, ply_dir, recursive=True, **kwargs):
    meshes_info = {}
    obj_list = (
        list(obj_dir.glob("**/*.obj"))
        if recursive
        else list(obj_dir.glob("*.obj"))
    )
    for obj_path in obj_list:
        ply_path = ply_dir.joinpath(obj_path.stem + ".ply")
        if ply_path.is_file():
            # read the ply file to get the center and n_tri
            mesh = o3d.io.read_triangle_mesh(str(ply_path))
            mesh_center = mesh.get_center()
            mesh_n_tri = len(mesh.triangles)
            mesh_z_min = np.asarray(mesh.vertices)[:, 2].min()
        else:
            mesh_center, mesh_n_tri, mesh_z_min = obj2ply_crs_conversion(
                obj_path, ply_path, **kwargs
            )

        id = obj_path.stem
        meshes_info[id] = {
            "filename": str(ply_path.relative_to(ply_dir.parent)),
            "center": mesh_center,
            "n_tri": mesh_n_tri,
            "id": id,
            "z_min": float(mesh_z_min),
        }
    return meshes_info


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
        [model_center_x, model_center_y, model_z]
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


def render_model(
    scene_dict, model_dict, model_name, sensor_name, center_view=True
):
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


def obj2ply(obj_path, ply_path, ft2m=True, center=False):
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

    # if center is True, center the mesh in (0,0,0)
    mesh_trans = [0, 0, 0]
    if center:
        mesh_center = mesh.get_center()
        if isinstance(center, tuple) and len(center) == 2:
            mesh_center[0] = center[0]
            mesh_center[1] = center[1]

        mesh_min_z = np.asarray(mesh.vertices)[:, 2].min()
        mesh_trans = [-mesh_center[0], -mesh_center[1], -mesh_min_z]
    else:
        mesh_center = mesh.get_center()
    mesh.translate(mesh_trans)

    # if ft2m is True, convert the mesh from feet to meters
    if ft2m:
        mesh.scale(FT2M_FACTOR, center=[0, 0, 0])

    # compute the normals
    mesh.compute_vertex_normals()
    mesh.compute_triangle_normals()

    # save the mesh as a ply file
    o3d.io.write_triangle_mesh(str(ply_path), mesh)

    mesh_n_tri = len(mesh.triangles)

    return mesh_trans, mesh_n_tri


def obj2ply_crs_conversion(obj_path, ply_path, transformer, flat=True):
    ## Convert the OBJ file to PLY, using the CRS transformer

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

    vertices = np.asarray(mesh.vertices)

    new_vertices = transformer.transform(
        vertices[:, 0], vertices[:, 1], vertices[:, 2], errcheck=True
    )
    new_vertices = np.array(new_vertices).T
    if np.all(np.abs(new_vertices[:, 2] - vertices[:, 2]) < 1e-6):
        # print("No change in Z values. Applying manual conversion.")
        # TODO: fix this
        new_vertices[:, 2] = vertices[:, 2] * FT2M_FACTOR

    z_min = new_vertices[:, 2].min()
    if flat:
        new_vertices[:, 2] = new_vertices[:, 2] - z_min

    new_vertices_dummy = np.zeros(new_vertices.shape)
    for v in range(new_vertices.shape[0]):
        for c in range(new_vertices.shape[1]):
            new_vertices_dummy[v, c] = new_vertices[v, c]

    mesh.vertices = o3d.utility.Vector3dVector(new_vertices_dummy)

    mesh_center = mesh.get_center()

    # compute the normals
    mesh.compute_vertex_normals()
    mesh.compute_triangle_normals()

    # save the mesh as a ply file
    o3d.io.write_triangle_mesh(str(ply_path), mesh)

    mesh_n_tri = len(mesh.triangles)

    return mesh_center, mesh_n_tri, z_min


def create_ground_dict(
    model_material,
    x_shift,
    y_shift,
    z_shift,
    scale_factor,
    base_rect_path,
    ply_path,
):
    ## create xml dict for the flat ground model

    # check if the base mesh exists
    base_rect_path = Path(base_rect_path)
    if not base_rect_path.is_file():
        raise FileNotFoundError(f"File not found: {base_rect_path}")
    # check if the base mesh exists in `ply_path`
    base_rect_ply_path = ply_path.parent.joinpath(base_rect_path.name)
    # copy the base mesh to `ply_path` if it doesn't exist
    if not base_rect_ply_path.is_file():
        base_rect_ply_path.write_bytes(base_rect_path.read_bytes())

    # the PLY file is saved in relative coordinates, centered in [0,0]
    # note that the unit is converted from feet to meters

    to_world = get_transformation_matrix(
        translation=[x_shift, y_shift, z_shift], scale=scale_factor
    )
    model_dict = {
        "type": "ply",
        "filename": str(base_rect_ply_path.resolve()),
        "face_normals": True,
        "to_world": to_world,
        "material_id": model_material,
    }

    mesh = read_mesh(base_rect_ply_path)
    n_tri = len(mesh.triangles)
    model_dict["n_tri"] = n_tri
    model_dict["z_min"] = 0.0

    return model_dict
