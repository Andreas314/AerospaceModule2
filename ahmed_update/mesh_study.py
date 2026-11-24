#!/usr/bin/env python3
"""mesh_study.py - mesh-independence study from summary.csv (S, M, L runs).

Writes mesh_study.txt and mesh_study.png and prints the table:
  * Cd, Cl, Cm per mesh and the deviation of Cd from the experiment (Lienhart & Becker 2003,
    25 deg slant, Re_H = 7.7e5: Cd ~ 0.285)
  * Grid Convergence Index (Celik et al., ASME JFE 2008) for Cd with the 3 meshes:
        r = (N_fine/N_coarse)^(1/3)           representative cell-size ratio
        p = |ln|e32/e21| + q| / ln r21        observed order (fixed-point iteration)
        Cd_ext = (r21^p Cd1 - Cd2) / (r21^p - 1)   Richardson extrapolation
        GCI21 = 1.25 |e21/Cd1| / (r21^p - 1)  uncertainty of the fine-mesh value
    Only meaningful if Cd changes monotonically S -> M -> L.
Usage: python3 mesh_study.py [summary.csv] [--ref 0.285]
"""
import argparse, csv, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    rows = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            r = {k: v for k, v in r.items() if k is not None}
            if r.get("mesh") in ("S", "M", "L") and r.get("cells") and r.get("Cd_mean"):
                rows[r["mesh"]] = r                      # last run of each mesh wins
    return rows


def gci(f, n):
    """f, n = [fine, medium, coarse]. Returns dict or a reason string."""
    f1, f2, f3 = f
    r21, r32 = (n[0] / n[1]) ** (1 / 3), (n[1] / n[2]) ** (1 / 3)
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return "differences are zero (fully converged or identical)"
    s = math.copysign(1.0, e32 / e21)
    if s < 0:
        return "oscillatory convergence (Cd not monotonic S->M->L): no GCI"
    p = 2.0
    for _ in range(100):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-6:
            break
        p = p_new
    ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    g21 = 1.25 * abs(e21 / f1) / (r21 ** p - 1)
    return dict(p=p, ext=ext, gci=g21, r21=r21, r32=r32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="summary.csv")
    ap.add_argument("--ref", type=float, default=0.285, help="experimental Cd")
    a = ap.parse_args()
    rows = load(a.csv)
    order = [m for m in ("L", "M", "S") if m in rows]           # fine -> coarse
    if len(order) < 2:
        raise SystemExit(f"need >= 2 of S/M/L in {a.csv}, found {list(rows)}")
    N = [float(rows[m]["cells"]) for m in order]
    Cd = [float(rows[m]["Cd_mean"]) for m in order]

    lines = [f"Mesh study ({', '.join(order)} = fine -> coarse), experiment Cd = {a.ref}", "",
             f"{'mesh':<5}{'cells':>11}{'Cd':>9}{'dev.exp':>9}{'Cl':>9}{'Cm':>9}{'y+ body':>9}{'time[s]':>9}{'Cd conv.':>10}"]
    for m in order:
        r = rows[m]
        g = lambda k, f="{:.4f}": f.format(float(r[k])) if r.get(k) not in (None, "") else "-"
        lines.append(f"{m:<5}{int(float(r['cells'])):>11,}{g('Cd_mean'):>9}"
                     f"{(float(r['Cd_mean'])/a.ref-1)*100:>+8.1f}%{g('Cl_mean'):>9}{g('Cm_mean'):>9}"
                     f"{g('yplus_avg_body','{:.2f}'):>9}{g('elapsed_s','{:.0f}'):>9}{r.get('Cd_converged',''):>10}")
    if len(order) == 3:
        res = gci(Cd, N)
        lines.append("")
        if isinstance(res, str):
            lines.append("GCI: " + res)
        else:
            lines += [f"r21 = {res['r21']:.2f}, r32 = {res['r32']:.2f}   observed order p = {res['p']:.2f}",
                      f"Richardson-extrapolated Cd = {res['ext']:.4f}",
                      f"GCI(fine) = {res['gci']*100:.2f} %   ->  Cd(L) = {Cd[0]:.4f} +- {res['gci']*100*Cd[0]/100:.4f}",
                      "(Note: the first-layer stack is fixed by y+, so the 3 meshes are not a purely "
                      "geometric refinement: read p and GCI as indicative.)"]
    txt = "\n".join(lines)
    print(txt)
    Path("mesh_study.txt").write_text(txt + "\n")

    h = [n ** (-1 / 3) for n in N]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(h, Cd, "o-", label="RANS k-omega SST")
    for hi, ci, m in zip(h, Cd, order):
        ax.annotate(f"{m} ({int(N[order.index(m)]/1000)}k)", (hi, ci), textcoords="offset points", xytext=(5, 6))
    ax.axhline(a.ref, color="g", ls="--", label=f"experiment {a.ref}")
    ax.axhspan(a.ref * 0.95, a.ref * 1.05, color="g", alpha=.1, label="+-5 %")
    ax.set(xlabel="cell-size proxy  N^(-1/3)  (-> 0 = finer)", ylabel="Cd [-]", title="Ahmed body 25 deg: mesh study")
    ax.set_xlim(left=0); ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig("mesh_study.png", dpi=150)


if __name__ == "__main__":
    main()
