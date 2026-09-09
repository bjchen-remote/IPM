#!/usr/bin/env python3
"""Read native MATLAB v7.3 checkpoint telemetry without MATLAB or scipy.

Monitoring only: no signature/norm/trusted-prefix validation, no rho read,
no PDE, no checkpoint mutation. Use only atomically installed checkpoints.
Example on this Mac:
  python3 research/longtime_lab/read_checkpoint_scalars.py CHECKPOINT.mat
The MATLAB HDF5 library is discovered locally; --library can select another.
"""

import argparse
import ctypes as C
import glob
import json
import math
from pathlib import Path


class HDF5:
    def __init__(self, library, filename):
        # Apple's protected Python strips DYLD_LIBRARY_PATH. Preload the
        # MATLAB sibling dependency by its absolute path instead.
        sibling = Path(library).parent / "libsz.2.dylib"
        self.dependencies = [C.CDLL(str(sibling), mode=C.RTLD_GLOBAL)] if sibling.is_file() else []
        self.lib = C.CDLL(library)
        hid, size, ptr = C.c_int64, C.c_uint64, C.c_void_p
        specs = {
            "H5open": (C.c_int, []),
            "H5Eset_auto2": (C.c_int, [hid, ptr, ptr]),
            "H5Fopen": (hid, [C.c_char_p, C.c_uint, hid]),
            "H5Fclose": (C.c_int, [hid]),
            "H5Dopen2": (hid, [hid, C.c_char_p, hid]),
            "H5Dclose": (C.c_int, [hid]),
            "H5Dget_space": (hid, [hid]),
            "H5Dget_type": (hid, [hid]),
            "H5Tget_class": (C.c_int, [hid]),
            "H5Tclose": (C.c_int, [hid]),
            "H5Sget_simple_extent_ndims": (C.c_int, [hid]),
            "H5Sget_simple_extent_dims": (C.c_int, [hid, ptr, ptr]),
            "H5Sselect_hyperslab": (C.c_int, [hid, C.c_int, ptr, ptr, ptr, ptr]),
            "H5Screate_simple": (hid, [C.c_int, ptr, ptr]),
            "H5Sclose": (C.c_int, [hid]),
            "H5Dread": (C.c_int, [hid, hid, hid, hid, hid, ptr]),
            "H5Gopen2": (hid, [hid, C.c_char_p, hid]),
            "H5Gclose": (C.c_int, [hid]),
        }
        for name, (result, arguments) in specs.items():
            function = getattr(self.lib, name)
            function.restype, function.argtypes = result, arguments
        self.lib.H5open()
        self.lib.H5Eset_auto2(0, None, None)
        self.native_double = hid.in_dll(self.lib, "H5T_NATIVE_DOUBLE_g").value
        self.file = self.lib.H5Fopen(str(filename).encode(), 0, 0)
        if self.file < 0:
            raise ValueError("Not a readable HDF5/MAT v7.3 file")
        self.callback_type = C.CFUNCTYPE(C.c_int, hid, C.c_char_p, ptr, ptr)
        self.iterate = getattr(self.lib, "H5Literate2", None) or self.lib.H5Literate
        self.iterate.restype = C.c_int
        self.iterate.argtypes = [hid, C.c_int, C.c_int, C.POINTER(size), self.callback_type, ptr]

    def close(self):
        self.lib.H5Fclose(self.file)

    def members(self, path):
        group = self.lib.H5Gopen2(self.file, path.encode(), 0)
        if group < 0:
            raise ValueError("Missing group: " + path)
        names = []
        callback = self.callback_type(lambda _g, name, _i, _d: names.append(name.decode()) or 0)
        index = C.c_uint64(0)
        try:
            if self.iterate(group, 0, 0, C.byref(index), callback, None) < 0:
                raise ValueError("Could not list " + path)
        finally:
            self.lib.H5Gclose(group)
        return names

    def numeric(self, path, history_length=None, last_only=False):
        dataset = self.lib.H5Dopen2(self.file, path.encode(), 0)
        if dataset < 0:
            raise ValueError("Not a numeric dataset: " + path)
        space = datatype = memory = -1
        try:
            datatype = self.lib.H5Dget_type(dataset)
            if self.lib.H5Tget_class(datatype) not in (0, 1):
                raise ValueError("Unsupported nonnumeric dataset: " + path)
            space = self.lib.H5Dget_space(dataset)
            rank = self.lib.H5Sget_simple_extent_ndims(space)
            if rank < 1 or rank > 4:
                raise ValueError("Unsupported dataset rank: " + path)
            shape = (C.c_uint64 * rank)()
            self.lib.H5Sget_simple_extent_dims(space, shape, None)
            counts, start = list(shape), [0] * rank
            if last_only:
                start, counts = [length - 1 for length in shape], [1] * rank
            elif history_length is not None:
                axes = [i for i, length in enumerate(shape) if length == history_length]
                if not axes:
                    raise ValueError("History length mismatch: " + path)
                # MATLAB NxK histories are HDF5 KxN; time is the last axis.
                time_axis = axes[-1]
                start[time_axis] = history_length - 1
                counts[time_axis] = 1
            total = math.prod(counts)
            if total < 1 or total > 10000:
                raise ValueError("Refusing large/empty telemetry read: " + path)
            start_c, count_c = (C.c_uint64 * rank)(*start), (C.c_uint64 * rank)(*counts)
            if self.lib.H5Sselect_hyperslab(space, 0, start_c, None, count_c, None) < 0:
                raise ValueError("Could not select telemetry: " + path)
            memory = self.lib.H5Screate_simple(rank, count_c, None)
            buffer = (C.c_double * total)()
            if self.lib.H5Dread(dataset, self.native_double, memory, space, 0, buffer) < 0:
                raise ValueError("Could not read telemetry: " + path)
            values = list(buffer)
            return (values[0] if total == 1 else values), list(shape)
        finally:
            if memory >= 0:
                self.lib.H5Sclose(memory)
            if space >= 0:
                self.lib.H5Sclose(space)
            if datatype >= 0:
                self.lib.H5Tclose(datatype)
            self.lib.H5Dclose(dataset)


def read_telemetry(filename, library):
    source = Path(filename).resolve(strict=True)
    hdf = HDF5(library, source)
    try:
        root = "/checkpoint"
        if "checkpoint" not in hdf.members("/"):
            raise ValueError("Expected native variable 'checkpoint'")
        state, history = root + "/payload/state", root + "/payload/log/history"
        # Read one value and the shape, never the complete history vector.
        _, step_shape = hdf.numeric(history + "/common/acceptedStep", last_only=True)
        if sum(length > 1 for length in step_shape) > 1:
            raise ValueError("acceptedStep must be a scalar history vector")
        length = max(step_shape)
        report = {"checkpointFile": str(source), "monitoringOnly": True,
                  "signatureValidated": False, "trustedPrefixValidated": False,
                  "historyLength": length, "state": {}, "historyLast": {},
                  "skippedFields": []}
        for name in ("step", "normalizedTime", "remeshCount"):
            report["state"][name] = hdf.numeric(state + "/" + name)[0]
        report["state"]["scale"] = {name: hdf.numeric(state + "/scale/" + name)[0]
                                          for name in ("canonicalTime", "physicalTime", "logC_l", "logC_omega")}
        for group in ("common", "mesh", "gauge"):
            values = report["historyLast"][group] = {}
            for name in hdf.members(history + "/" + group):
                try:
                    values[name] = hdf.numeric(history + "/" + group + "/" + name, length)[0]
                except ValueError as error:
                    report["skippedFields"].append({"field": group + "." + name, "reason": str(error)})
        return report
    finally:
        hdf.close()


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite": str(value)}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint")
    parser.add_argument("--library", help="Path to an HDF5 shared library")
    arguments = parser.parse_args()
    candidates = sorted(glob.glob("/Applications/MATLAB_*.app/bin/maca64/libhdf5.*.dylib"), reverse=True)
    library = arguments.library or (candidates[0] if candidates else None)
    if library is None:
        parser.error("No MATLAB HDF5 library found; provide --library")
    try:
        result = read_telemetry(arguments.checkpoint, library)
    except (OSError, ValueError) as error:
        parser.exit(1, "Monitoring read failed: " + str(error) + "\n")
    print(json.dumps(json_safe(result), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
