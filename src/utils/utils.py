import math
import os
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
