#!/usr/bin/env python3
"""Validation: compare spot centroids measured in the stitched volume with the
positions in a legacy dtfind .ref list.  A correct stitch gives median offsets
near zero in px0, px1 and channel.

  check_stitch.py file.h5 legacy_find.ref [--n 200]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pcs_stitch import stitch  # noqa: E402
from refcmp import read_ref    # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('h5')
    ap.add_argument('ref')
    ap.add_argument('--n', type=int, default=200, help='use the N strongest legacy spots')
    ap.add_argument('--half', type=int, default=3, help='half-width of the centroid box (pixels)')
    args = ap.parse_args()

    vol, _ = stitch(args.h5)
    nchan, n1, n0 = vol.shape
    _, c = read_ref(args.ref)
    order = np.argsort(-c['fIntensity'])[:args.n]
    h = args.half
    offs = []
    for i in order:
        p0, p1, ch = c['fObs_pixel0'][i], c['fObs_pixel1'][i], c['fObs_rot_mid'][i]
        a0, a1, ac = int(p0), int(p1), int(ch)
        if not (h <= a0 < n0 - h and h <= a1 < n1 - h and 2 <= ac < nchan - 2):
            continue
        box = vol[ac - 2:ac + 3, a1 - h:a1 + h + 1, a0 - h:a0 + h + 1].astype(np.float64)
        bg = vol[ac - 2:ac + 3, max(a1 - 3 * h, 0):a1 + 3 * h + 1, max(a0 - 3 * h, 0):a0 + 3 * h + 1]
        box = box - np.median(bg)
        box[box < 0] = 0
        if box.sum() <= 0:
            continue
        gc, g1, g0 = np.mgrid[ac - 2:ac + 3, a1 - h:a1 + h + 1, a0 - h:a0 + h + 1]
        w = box.sum()
        offs.append(((g0 * box).sum() / w + 0.5 - p0, (g1 * box).sum() / w + 0.5 - p1,
                     (gc * box).sum() / w + 0.5 - ch))
    offs = np.array(offs)
    med = np.median(offs, axis=0)
    mad = np.median(np.abs(offs - med), axis=0)
    print(f'{os.path.basename(args.h5)}: n={len(offs)}  median offset px0={med[0]:+.2f} px1={med[1]:+.2f} '
          f'chan={med[2]:+.2f}   MAD px0={mad[0]:.2f} px1={mad[1]:.2f} chan={mad[2]:.2f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
