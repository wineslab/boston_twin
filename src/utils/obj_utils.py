import mitsuba as mi

def obj2ply_mi(
    model_material,
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
            0,  # -z_min_m_from_info,
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


def create_ground(
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

def get_mi_dict(ply_path,x_shift,y_shift,z_shift,ply_path_relative,model_material):
    in_model_dict = {
        "type": "ply",
        "filename": str(ply_path.resolve()),
        "face_normals": True,
    }
    model = mi.load_dict(in_model_dict)  # load the model into mitsuba
    n_tri = model.face_count()

    # create the new mitsuba dict for the PLY file
    center_x_m_from_info = x_shift
    center_y_m_from_info = y_shift
    z_min_m_from_info = z_shift
    to_world = mi.ScalarTransform4f.translate(
        [
            -center_x_m_from_info,
            -center_y_m_from_info,
            0,  # -z_min_m_from_info,
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
