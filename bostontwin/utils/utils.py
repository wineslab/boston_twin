import math
import os
from pathlib import Path
from typing import Union

import mitsuba as mi

try:
    mi.set_variant("cuda_ad_rgb")
except:
    try:
        mi.set_variant("llvm_ad_rgb")
    except:
        mi.set_variant("scalar_rgb")

from .constants import FT2M_FACTOR
from .obj_utils import create_ground_dict, get_mi_dict


def generate_mi_xml(
    scene_name,
    models_list,
    models_dir,
    out_dir,
    models_materials,
    models_center,
    create_ground=True,
    ground_size = FT2M_FACTOR * 2650
):
    
    scene_dict = generate_mi_scene_dict(models_list, models_dir, out_dir, models_materials, models_center, create_ground)
    if create_ground:
        # load frame obj and use it as (flat) ground
        frame_name = f"{scene_name}_terrain_flat"
        base_rect_path = models_dir.joinpath("rectangle.ply")
        frame_material = "mat-itu_medium_dry_ground"
        ground_dict, _ = create_ground_dict(
            frame_material,
            0,
            0,
            0,
            ground_size,  # tile size is 5000 ft x 5000 ft, we increase it a bit
            base_rect_path.resolve(),
            models_dir.joinpath(frame_name + ".ply"),
            out_dir,
        )
        scene_dict["ground"] = ground_dict
        mi.xml.dict_to_xml(scene_dict, str(out_dir.joinpath(scene_name + ".xml").resolve()))

def generate_mi_scene_dict(models_list, models_dir, out_dir, models_materials, models_center, create_ground=True):
    mitsuba_scene_dict = {
        "type": "scene",
        "integrator": {
            "type": "path",
        },
        "light": {"type": "constant"},
        "mat-itu_brick": {
            "type": "twosided",
            "bsdf": {
                "type": "diffuse",
                "reflectance": {
                    "type": "rgb",
                    "value": [0.401968, 0.111874, 0.086764],
                },
            },
        },
        "mat-itu_concrete": {
            "type": "twosided",
            "bsdf": {
                "type": "diffuse",
                "reflectance": {
                    "type": "rgb",
                    "value": [0.539479, 0.539479, 0.539480],
                },
            },
        },
        "mat-itu_medium_dry_ground": {
            "type": "twosided",
            "bsdf": {
                "type": "diffuse",
                "reflectance": {
                    "type": "rgb",
                    "value": [65 / 255, 60 / 255, 60 / 255],
                },
            },
        },
    }
    for model_name, model_material, model_center in zip(
        models_list, models_materials, models_center
    ):
        ply_path = models_dir.joinpath(model_name + ".ply")
        
        # create the model dictionary for Mitsuba
        tile_model_dict = get_mi_dict(
            model_material,
            model_center_x=model_center[0],
            model_center_y=model_center[1],
            model_z=0,  # model_info_from_json["Z_MIn_Ft"]*FT2M_FACTOR,
            ply_path=ply_path,
            ply_path_relative=out_dir,
        )

        # add model to the scene
        mitsuba_scene_dict[model_name] = tile_model_dict
    
    return mitsuba_scene_dict



def str2float(x):
    try:
        return float(x)
    except ValueError:
        if not x:
            return math.nan
        else:
            raise ValueError


def check_file_exists(filepath: Union[Path, str]):
    if not filepath.is_file():
        raise FileNotFoundError(
            f"Missing file {filepath.relative_to(filepath.parents[-2])}"
        )


def time2str(t: float, precision: str = "") -> str:
    mins = int(t // 60)
    secs_float = t % 60
    secs = int(secs_float)
    time_str = f"{mins} min {secs} s"
    if precision == "ms":
        ms = round(secs_float * 1e3)
        time_str + f" {ms} ms"
    return time_str


def print_eta(
    t0: float, t1: float, times: float, iter_idx: int, n_iters: int, **kwargs
) -> None:
    t10 = t1 - t0
    t10_str = time2str(t10)
    times.append(t10)
    mean_time = sum(times) / (iter_idx + 1)
    iter_left = n_iters - (iter_idx + 1)
    t_left = iter_left * mean_time
    t_left_str = time2str(t_left)

    print(
        f"Iter {iter_idx+1}/{n_iters} completed in {t10_str}. {t_left_str} left.",
        **kwargs,
    )


def truncate_utf8_chars(filename, count, ignore_newlines=True):
    """
    Truncates last `count` characters of a text file encoded in UTF-8.
    :param filename: The path to the text file to read
    :param count: Number of UTF-8 characters to remove from the end of the file
    :param ignore_newlines: Set to true, if the newline character at the end of the file should be ignored
    """
    with open(filename, "rb+") as f:
        last_char = None

        size = os.fstat(f.fileno()).st_size

        offset = 1
        chars = 0
        while offset <= size:
            f.seek(-offset, os.SEEK_END)
            b = ord(f.read(1))

            if ignore_newlines:
                if b == 0x0D or b == 0x0A:
                    offset += 1
                    continue

            if b & 0b10000000 == 0 or b & 0b11000000 == 0b11000000:
                # This is the first byte of a UTF8 character
                chars += 1
                if chars == count:
                    # When `count` number of characters have been found, move current position back
                    # with one byte (to include the byte just checked) and truncate the file
                    f.seek(-1, os.SEEK_CUR)
                    f.truncate()
                    return
            offset += 1
