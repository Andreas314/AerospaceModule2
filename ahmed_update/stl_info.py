#!/usr/bin/env python3
"""stl_info.py <stl> [region]  ->  prints "L Aref"

L    = x-extent of the whole STL (reference length, flow along +x)
Aref = frontal area of <region> (default 'body') = 1/2 * sum(|n_x| * A) over its
       triangles. For a half model this is automatically the half-model area,
       consistent with the (half) forces that forceCoeffs integrates.
ASCII STL with named solids (solid body / solid legs); binary STL -> whole surface.
"""
import re, struct, sys
from pathlib import Path

def tris_ascii(text, region):
    parts = re.split(r"(?m)^\s*solid\s+(\S+)\s*$", text)
    out = {}
    for i in range(1, len(parts), 2):
        v = [tuple(map(float, m)) for m in
             re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", parts[i + 1])]
        out[parts[i]] = [v[j:j + 3] for j in range(0, len(v) - 2, 3)]
    return out

def main():
    path = sys.argv[1]
    region = sys.argv[2] if len(sys.argv) > 2 else "body"
    data = Path(path).read_bytes()
    n = struct.unpack("<I", data[80:84])[0] if len(data) >= 84 else -1
    if len(data) == 84 + 50 * n:                       # binary
        t = []
        for i in range(n):
            f = struct.unpack("<9f", data[84 + 50 * i + 12: 84 + 50 * i + 48])
            t.append([f[0:3], f[3:6], f[6:9]])
        regs = {"all": t}
    else:
        regs = tris_ascii(data.decode("utf-8", "ignore"), region)
    alltr = [tr for r in regs.values() for tr in r]
    xs = [p[0] for tr in alltr for p in tr]
    L = max(xs) - min(xs)
    sel = regs.get(region, alltr)
    s = 0.0
    for a, b, c in sel:
        uy, uz = b[1] - a[1], b[2] - a[2]
        wy, wz = c[1] - a[1], c[2] - a[2]
        s += abs(uy * wz - uz * wy) / 2          # |n_x| * A
    print(f"{L:.6g} {s / 2:.6g}")

if __name__ == "__main__":
    main()
