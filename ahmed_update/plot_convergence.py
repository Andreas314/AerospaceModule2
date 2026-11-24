#!/usr/bin/env python3
"""Convergence analysis of the force coefficients (Cd, Cl, Cm) + residuals + y+.

Reads  postProcessing/forceCoeffs/*/forceCoeffs.dat  and  log.foamRun, writes
  convergence.png, convergence_report.txt, and appends a row to summary.csv.

Per coefficient (window = last 20 % of the iterations, tol = 1 %):
  mean, std, relative std, drift (mean of last window vs the window before),
  linear slope over the window, and the iteration after which the running value stays
  within +-tol of the final mean.  CONVERGED = relative std < tol AND drift < tol.
"""
import argparse, csv, glob, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = re.compile(r"Solving for (\S+), Initial residual = ([^,]+),")
YP = re.compile(r"patch (\S+) y\+ : min = (\S+), max = (\S+), average = (\S+)")
WIN, TOL = 0.20, 0.01


def parse_log(path):
    res, step, yp = defaultdict(dict), 0, {}
    if not Path(path).exists():
        return res, yp
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("Time = "):
            step += 1
        elif (m := RES.search(line)):
            res[m.group(1)][step] = float(m.group(2))
        elif (m := YP.search(line)):
            yp[m.group(1)] = tuple(float(x.strip(",")) for x in m.groups()[1:])
    return res, yp


def parse_forces(case):
    files = sorted(glob.glob(f"{case}/postProcessing/forceCoeffs/*/forceCoeffs.dat"),
                   key=lambda f: float(Path(f).parent.name) if Path(f).parent.name.replace('.', '', 1).isdigit() else 0)
    rows, cols = [], None
    for f in files:                                    # restarts -> several time dirs
        for line in Path(f).read_text().splitlines():
            if line.startswith("#"):
                t = line.lstrip("# ").split()
                if t and t[0] == "Time":
                    cols = t
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    if not rows or cols is None:
        return None, None
    n = min(len(cols), min(map(len, rows)))
    return cols[:n], np.array([r[:n] for r in rows])


def analyse(t, y):
    n = len(y)
    w = max(10, int(n * WIN))
    last, prev = y[-w:], y[-2 * w:-w] if n >= 2 * w else y[:w]
    mean, std = last.mean(), last.std()
    scale = abs(mean) if abs(mean) > 1e-12 else 1.0
    drift = abs(mean - prev.mean()) / scale
    slope = np.polyfit(t[-w:], last, 1)[0] * (t[-1] - t[-w]) / scale   # relative change over window
    out = np.nonzero(np.abs(y - mean) > TOL * scale)[0]
    settle = t[min(out[-1] + 1, n - 1)] if len(out) else t[0]
    conv = (std / scale < TOL) and (drift < TOL)
    return dict(mean=mean, std=std, rstd=std / scale, drift=drift, slope=slope,
                settle=settle, converged=conv, window=w)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    for a in ("cells", "subdomains", "elapsed", "mesh", "yplus", "U"):
        ap.add_argument(f"--{a}", default="")
    a = ap.parse_args()
    case = Path(a.case)

    res, yp = parse_log(case / "log.foamRun")
    # y+ comes from the yPlus function object: postProcessing/yPlus/<t0>/yPlus.dat
    # (columns: time patch min max average). Last time wins; the moving ground is skipped.
    ydat = sorted((case / "postProcessing" / "yPlus").glob("*/yPlus.dat"))
    if ydat:
        last = {}
        for line in ydat[-1].read_text().splitlines():
            f = line.split()
            if len(f) >= 5 and not line.lstrip().startswith("#"):
                try:
                    last[f[1]] = (float(f[2]), float(f[3]), float(f[4]))
                except ValueError:
                    pass
        yp = {k: v for k, v in last.items() if k != "bottom"} or yp
    cols, F = parse_forces(case)
    if F is None or len(F) < 20:
        print("No (or too few) forceCoeffs data found in postProcessing/forceCoeffs")
        return
    t = F[:, 0]
    coefs = {c: analyse(t, F[:, cols.index(c)]) for c in ("Cd", "Cl", "Cm") if c in cols}

    fig, ax = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    for ci, ca in zip(("Cd", "Cl"), ax[:2]):
        if ci not in coefs:
            continue
        y, s = F[:, cols.index(ci)], coefs[ci]
        ca.plot(t, y, lw=1, label=ci)
        ca.plot(t, np.cumsum(y) / np.arange(1, len(y) + 1), lw=1, ls="--", label="running mean")
        ca.axhline(s["mean"], color="k", lw=.7)
        ca.axhspan(s["mean"] * (1 - TOL), s["mean"] * (1 + TOL), color="g", alpha=.15, label="+-1 % of final mean")
        ca.axvspan(t[-s["window"]], t[-1], color="gray", alpha=.12, label="analysis window")
        ca.set_ylabel(ci + " [-]"); ca.grid(alpha=.3); ca.legend(ncol=4, fontsize=8)
        ca.set_title(f"{ci} = {s['mean']:.4f}  (std {s['rstd']*100:.2f} %, drift {s['drift']*100:.2f} %)"
                     f"  {'CONVERGED' if s['converged'] else 'NOT converged'}", fontsize=10)
    for name, d in sorted(res.items()):
        ax[2].semilogy(list(d), list(d.values()), lw=1, label=name)
    ax[2].set(xlabel="iteration", ylabel="initial residual [-]"); ax[2].legend(ncol=4, fontsize=8)
    ax[2].grid(True, which="both", alpha=.3)
    fig.tight_layout(); fig.savefig(case / "convergence.png", dpi=150)

    L = [f"Convergence report  mesh={a.mesh} cells={a.cells} U={a.U} y+_target={a.yplus}",
         f"{len(t)} iterations, analysis window = last {int(WIN*100)} %, tolerance = {TOL*100:.0f} %", "",
         f"{'coef':<5}{'mean':>10}{'std[%]':>9}{'drift[%]':>10}{'slope[%]':>10}{'settled@':>10}  verdict"]
    for c, s in coefs.items():
        L.append(f"{c:<5}{s['mean']:>10.4f}{s['rstd']*100:>9.3f}{s['drift']*100:>10.3f}"
                 f"{s['slope']*100:>10.3f}{s['settle']:>10.0f}  {'CONVERGED' if s['converged'] else 'NOT converged'}")
    if res:
        L += ["", "final initial residuals: " + "  ".join(f"{k}={list(v.values())[-1]:.2e}" for k, v in sorted(res.items()))]
    for p, (mn, mx, av) in yp.items():
        L.append(f"y+ on {p}: min {mn:.2f}  avg {av:.2f}  max {mx:.2f}   (target {a.yplus})")
    (case / "convergence_report.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))

    row = dict(mesh=a.mesh, cells=a.cells, subdomains=a.subdomains, yplus_target=a.yplus, U=a.U,
               elapsed_s=a.elapsed, iterations=len(t))
    for c, s in coefs.items():
        row.update({f"{c}_mean": s["mean"], f"{c}_std": s["std"], f"{c}_converged": s["converged"]})
    for p, (mn, mx, av) in yp.items():
        row[f"yplus_avg_{p}"] = av
    out = case / "summary.csv"                 
    old, names = [], []                        
    if out.exists():
        with out.open(newline="") as fh:
            rd = csv.DictReader(fh)
            names = list(rd.fieldnames or [])
            old = [{k: v for k, v in r.items() if k is not None} for r in rd]
    names += [k for k in row if k not in names]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=names, restval="")
        w.writeheader()
        w.writerows(old + [row])


if __name__ == "__main__":
    main()
