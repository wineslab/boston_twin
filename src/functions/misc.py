import math
from typing import Union
from pathlib import Path

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
