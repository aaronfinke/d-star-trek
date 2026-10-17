#!/usr/bin/env python3
"""Compare two d*TREK reflection lists (.ref), e.g. new output vs. the output
of the original PCS TOF-Laue d*TREK.

Modes
  hkl    join on (h, k, l); report count, rms/mean differences of numeric
         fields and the correlation of fIntensity (optionally for I/sigma > cut)
  spots  match unindexed spots by nearest (px0, px1, channel) within tolerances;
         report the fraction recovered and the mean/rms offsets

Usage
  refcmp.py new.ref legacy.ref --mode hkl --fields fCalc_pixel0,fCalc_pixel1,fCalwavelength
  refcmp.py new.ref legacy.ref --mode spots --tol 2,2,1
"""

import argparse
import sys

import numpy as np


def read_ref(path):
    """Return (names, columns) where columns maps field name -> numpy array
    (float for n/f fields, object for s fields)."""
    with open(path) as f:
        lines = f.read().split('\n')
    counts = [int(x) for x in lines[0].split()]
    n_int, n_flt, n_str = counts[:3]
    n_head = counts[3] if len(counts) > 3 else 0   # embedded header lines
    nf = n_int + n_flt + n_str
    start = 1 + n_head
    names = [l.strip() for l in lines[start:start + nf]]
    rows = [l.split() for l in lines[start + nf:] if l.strip()]
    rows = [r for r in rows if len(r) >= nf]
    cols = {}
    for i, name in enumerate(names):
        vals = [r[i] for r in rows]
        if i < n_int + n_flt:
            cols[name] = np.array(vals, dtype=float)
        else:
            cols[name] = np.array(vals, dtype=object)
    return names, cols


def stats(label, d):
    d = d[np.isfinite(d)]
    if d.size == 0:
        print(f'  {label:24s}  n=0')
        return
    print(f'  {label:24s}  n={d.size:6d}  mean={d.mean():+10.4f}  rms={np.sqrt(np.mean(d * d)):10.4f}'
          f'  median|d|={np.median(np.abs(d)):10.4f}')


def compare_hkl(a, b, fields, isig):
    def keys(c):
        return list(zip(c['nH'].astype(int), c['nK'].astype(int), c['nL'].astype(int)))
    ka, kb = keys(a), keys(b)
    ib = {}
    for j, k in enumerate(kb):
        ib.setdefault(k, j)
    pairs = [(i, ib[k]) for i, k in enumerate(ka) if k in ib]
    print(f'reflections: new={len(ka)} legacy={len(kb)} common hkl={len(pairs)}')
    if not pairs:
        return
    ia = np.array([p[0] for p in pairs])
    jb = np.array([p[1] for p in pairs])
    for f in fields:
        if f in a and f in b:
            va, vb = a[f][ia], b[f][jb]
            ok = (va > -998.0) & (vb > -998.0)
            stats(f, va[ok] - vb[ok])
        else:
            print(f'  {f:24s}  missing in {"new" if f not in a else "legacy"}')
    if 'fIntensity' in a and 'fIntensity' in b:
        Ia, Ib = a['fIntensity'][ia], b['fIntensity'][jb]
        sb = b['fSigmaI'][jb] if 'fSigmaI' in b else np.ones_like(Ib)
        ok = (sb > 0) & (Ib / np.where(sb > 0, sb, 1) > isig)
        if ok.sum() > 2:
            cc = np.corrcoef(Ia[ok], Ib[ok])[0, 1]
            ratio = np.median(Ia[ok] / np.where(Ib[ok] != 0, Ib[ok], np.nan))
            print(f'  fIntensity CC (legacy I/sig > {isig:g}) = {cc:.4f}  n={ok.sum()}  median ratio new/legacy={ratio:.4f}')


def compare_spots(a, b, tol, rotfield):
    pa = np.vstack([a['fObs_pixel0'], a['fObs_pixel1'], a[rotfield]]).T
    pb = np.vstack([b['fObs_pixel0'], b['fObs_pixel1'], b[rotfield]]).T
    t = np.array(tol, dtype=float)
    found = []
    for j, p in enumerate(pb):
        d = np.abs(pa - p) / t
        m = np.max(d, axis=1)
        i = int(np.argmin(m))
        if m[i] <= 1.0:
            found.append((i, j))
    print(f'spots: new={len(pa)} legacy={len(pb)} legacy recovered={len(found)} '
          f'({100.0 * len(found) / max(len(pb), 1):.1f}%) within tol {tol}')
    if found:
        ia = np.array([f[0] for f in found])
        jb = np.array([f[1] for f in found])
        for k, name in enumerate(['fObs_pixel0', 'fObs_pixel1', rotfield]):
            stats(name, pa[ia, k] - pb[jb, k])
        for name in ('fDetwavelength', 'fIntensity'):
            if name in a and name in b:
                stats(name, a[name][ia] - b[name][jb])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('new')
    ap.add_argument('legacy')
    ap.add_argument('--mode', choices=['hkl', 'spots'], default='hkl')
    ap.add_argument('--fields', default='fCalc_pixel0,fCalc_pixel1,fCalc_rot_mid,fCalwavelength,fObs_pixel0,fObs_pixel1,fObs_rot_mid')
    ap.add_argument('--tol', default='2,2,1', help='spots mode: px0,px1,channel tolerances')
    ap.add_argument('--isig', type=float, default=3.0, help='hkl mode: legacy I/sigma cut for the CC')
    args = ap.parse_args()

    _, a = read_ref(args.new)
    _, b = read_ref(args.legacy)
    if args.mode == 'hkl':
        compare_hkl(a, b, [f for f in args.fields.split(',') if f], args.isig)
    else:
        compare_spots(a, b, [float(x) for x in args.tol.split(',')], 'fObs_rot_mid')
    return 0


if __name__ == '__main__':
    sys.exit(main())
