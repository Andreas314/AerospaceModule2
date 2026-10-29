#!/usr/bin/env python3
"""
make_background_box.py

Reads the bounding box of a solid body (STL, ASCII or binary), builds a
background box around it for blockMesh, and computes a locationInMesh point
for snappyHexMesh that is inside the box but outside the body.

Flow is assumed along +x (inlet at x-min, outlet at x-max).

Example (Ahmed body, domain sizes in multiples of body length L):
    python make_background_box.py --stl constant/triSurface/ahmed.stl \
        --upstream 3 --downstream 8 --lateral 3 --top 3 --bottom 0.5 \
        --cell-size 0.05

Example (absolute sizes in metres):
    python make_background_box.py --stl ahmed.stl --absolute \
        --upstream 1 --downstream 4 --lateral 1.5 --top 1.5 --bottom 0.3 \
        --cell-size 0.05
"""
import argparse
import math
import re
import struct
from pathlib import Path


# ----------------------------------------------------------------- STL --
def stl_bounds(path):
    """Return (xmin, ymin, zmin, xmax, ymax, zmax) of an ASCII or binary STL."""
    data = Path(path).read_bytes()
    lo = [math.inf] * 3
    hi = [-math.inf] * 3

    def update(x, y, z):
        for i, v in enumerate((x, y, z)):
            lo[i] = min(lo[i], v)
            hi[i] = max(hi[i], v)

    # Binary STL: 80-byte header, uint32 n_tri, then 50 bytes per triangle
    is_binary = False
    if len(data) >= 84:
        n_tri = struct.unpack("<I", data[80:84])[0]
        is_binary = len(data) == 84 + 50 * n_tri

    if is_binary:
        for i in range(n_tri):
            off = 84 + 50 * i + 12          # skip the normal
            v = struct.unpack("<9f", data[off:off + 36])
            for k in range(0, 9, 3):
                update(*v[k:k + 3])
    else:
        text = data.decode("utf-8", errors="ignore")
        for m in re.finditer(
            r"vertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", text
        ):
            update(*map(float, m.groups()))

    if math.inf in lo:
        raise ValueError(f"No vertices found in {path}")
    return (*lo, *hi)


# --------------------------------------------------------- domain setup --
def build_domain(bb, up, down, lat, top, bot):
    xmin, ymin, zmin, xmax, ymax, zmax = bb
    return (
        xmin - up, ymin - lat, zmin - bot,
        xmax + down, ymax + lat, zmax + top,
    )


def inside_point(bb, dom):
    """
    Pick the middle of the largest gap between body bbox and domain wall.
    That point is inside the domain and outside the body bbox, so it is
    always in the fluid. A small offset keeps it off cell faces/edges.
    """
    xmin, ymin, zmin, xmax, ymax, zmax = bb
    dxmin, dymin, dzmin, dxmax, dymax, dzmax = dom
    cx, cy, cz = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2

    candidates = [  # (gap size, point)
        (xmin - dxmin, ((xmin + dxmin) / 2, cy, cz)),
        (dxmax - xmax, ((xmax + dxmax) / 2, cy, cz)),
        (ymin - dymin, (cx, (ymin + dymin) / 2, cz)),
        (dymax - ymax, (cx, (ymax + dymax) / 2, cz)),
        (zmin - dzmin, (cx, cy, (zmin + dzmin) / 2)),
        (dzmax - zmax, (cx, cy, (zmax + dzmax) / 2)),
    ]
    gap, pt = max(candidates, key=lambda c: c[0])
    if gap <= 0:
        raise ValueError("Domain does not extend beyond the body in any direction.")

    # tiny irrational-ish jitter so the point never sits exactly on a face
    jitter = 1.2345e-3 * gap
    return (pt[0] + jitter, pt[1] + 0.7 * jitter, pt[2] + 0.3 * jitter)


# --------------------------------------------------------------- output --
def cells(length, h):
    return max(1, round(length / h))


def write_blockmesh(path, dom, h):
    x0, y0, z0, x1, y1, z1 = dom
    nx, ny, nz = cells(x1 - x0, h), cells(y1 - y0, h), cells(z1 - z0, h)
    txt = f"""FoamFile
{{
    version 2.0;
    format  ascii;
    class   dictionary;
    object  blockMeshDict;
}}

convertToMeters 1;

vertices
(
    ({x0:.6g} {y0:.6g} {z0:.6g})  // 0
    ({x1:.6g} {y0:.6g} {z0:.6g})  // 1
    ({x1:.6g} {y1:.6g} {z0:.6g})  // 2
    ({x0:.6g} {y1:.6g} {z0:.6g})  // 3
    ({x0:.6g} {y0:.6g} {z1:.6g})  // 4
    ({x1:.6g} {y0:.6g} {z1:.6g})  // 5
    ({x1:.6g} {y1:.6g} {z1:.6g})  // 6
    ({x0:.6g} {y1:.6g} {z1:.6g})  // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    bottom {{ type patch; faces ((0 3 2 1)); }}
    top    {{ type patch; faces ((4 5 6 7)); }}
    sides  {{ type patch; faces ((0 1 5 4) (3 7 6 2)); }}
);

mergePatchPairs ();
"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(txt)
    return nx, ny, nz


def patch_snappy(path, pt):
    """Replace locationInMesh in an existing snappyHexMeshDict, if present."""
    p = Path(path)
    if not p.exists():
        return False
    s = p.read_text()
    new = f"locationInMesh ({pt[0]:.6g} {pt[1]:.6g} {pt[2]:.6g});"
    s2, n = re.subn(r"locationInMesh\s*\([^)]*\)\s*;", new, s)
    if n:
        p.write_text(s2)
    return bool(n)


# ----------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stl", required=True, help="STL of the solid body")
    ap.add_argument("--upstream", type=float, required=True, help="gap inlet -> body")
    ap.add_argument("--downstream", type=float, required=True, help="gap body -> outlet")
    ap.add_argument("--lateral", type=float, required=True, help="gap body -> side walls")
    ap.add_argument("--top", type=float, required=True, help="gap body -> top wall")
    ap.add_argument("--bottom", type=float, required=True, help="gap body -> bottom wall")
    ap.add_argument("--absolute", action="store_true",
                    help="gaps are in metres (default: multiples of body length L = x-extent)")
    ap.add_argument("--cell-size", type=float, required=True,
                    help="target background cell size (m)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="STL unit scale to metres (e.g. 0.001 for mm)")
    ap.add_argument("--case", default=".", help="OpenFOAM case directory")
    args = ap.parse_args()

    bb = tuple(v * args.scale for v in stl_bounds(args.stl))
    L = bb[3] - bb[0]
    f = 1.0 if args.absolute else L
    gaps = [g * f for g in (args.upstream, args.downstream, args.lateral,
                            args.top, args.bottom)]

    dom = build_domain(bb, *gaps)
    pt = inside_point(bb, dom)

    case = Path(args.case)
    nx, ny, nz = write_blockmesh(case / "system" / "blockMeshDict", dom, args.cell_size)
    patched = patch_snappy(case / "system" / "snappyHexMeshDict", pt)

    print(f"Body bbox  : {bb[:3]} -> {bb[3:]}  (L = {L:.4g} m)")
    print(f"Domain     : {dom[:3]} -> {dom[3:]}")
    print(f"Cells      : {nx} x {ny} x {nz} = {nx*ny*nz:,}")
    print(f"locationInMesh ({pt[0]:.6g} {pt[1]:.6g} {pt[2]:.6g});")
    print("snappyHexMeshDict updated." if patched else
          "Copy the locationInMesh line above into snappyHexMeshDict.")


if __name__ == "__main__":
    main()
