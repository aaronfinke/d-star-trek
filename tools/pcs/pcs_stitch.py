#!/usr/bin/env python3
"""Reference (numpy) implementation of the PCS NeXus segment stitching used by
d*TREK's nexuspcs.cpp, and a reader for d*TREK .img files.  Used to check the
C++ reader and to make quick projections.

  pcs_stitch.py file.h5 --compare image.img --chan 20
"""

import argparse
import sys

import h5py
import numpy as np


def stitch(path):
    """Return (volume[chan][px1][px0] float32, tof_edges) for a PCS HDF5 file."""
    f = h5py.File(path, 'r')
    run = [k for k in f if k.startswith('run_')][0]
    det = f[run]['NXinstrument/NXdetector1']
    segs = []
    k = 1
    while f'{run}/tof_detector_{k}/data' in f:
        segs.append(np.asarray(f[f'{run}/tof_detector_{k}/data'][()], dtype=np.float32))
        k += 1
    ntof, narc, nvert = segs[0].shape
    ratio = 1.0
    if 'horizontal_segment_size_in_pixels' in det:
        ratio = float(det['horizontal_segment_size_in_pixels'][()][0]) / narc
    overlap = (16 * narc) // 256
    if 'horizontal_segment_overlap_in_pixels' in det:
        overlap = int(det['horizontal_segment_overlap_in_pixels'][()][0] / ratio + 0.5)
    cutoff = (8 * narc) // 256
    if 'horizontal_border_cutoff_in_pixels' in det:
        cutoff = int(det['horizontal_border_cutoff_in_pixels'][()][0] / ratio + 0.5)
    stride = narc - overlap
    dim1 = stride * len(segs) + overlap
    if 'detector dimensions' in det:
        dim1 = max(dim1, int(det['detector dimensions'][()][1]))
    vol = np.zeros((ntof, dim1, nvert), dtype=np.float32)
    for s, d in enumerate(segs):
        a = np.arange(cutoff, narc - cutoff)
        p1 = s * stride + (a - cutoff)
        ok = (p1 >= 0) & (p1 < dim1)
        vol[:, p1[ok], :] = d[:, a[ok], :]
    tof = np.asarray(f[f'{run}/tof_detector_1/tof'][()], dtype=float)
    return vol, tof


def read_img(path):
    """Read a d*TREK .img with an uncompressed float or unsigned short body."""
    with open(path, 'rb') as fh:
        raw = fh.read()
    head = raw[:4096].split(b'}')[0].decode(errors='replace')
    kv = {}
    for item in head.replace('\n', '').strip('{').split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            kv[k.strip()] = v.strip()
    nbytes = int(kv['HEADER_BYTES'])
    n0, n1 = int(kv['SIZE1']), int(kv['SIZE2'])
    order = '<' if kv.get('BYTE_ORDER', 'little_endian').startswith('little') else '>'
    dt = order + ('f4' if kv['Data_type'].startswith('float') else 'u2')
    data = np.frombuffer(raw, dtype=dt, count=n0 * n1, offset=nbytes).reshape(n1, n0)
    return kv, data


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('h5')
    ap.add_argument('--compare', help='d*TREK image written by dtpcsproject')
    ap.add_argument('--chan', type=int, default=-1, help='channel; -1 = sum, -2 = max')
    args = ap.parse_args()
    vol, tof = stitch(args.h5)
    print(f'volume {vol.shape}, total counts {vol.sum():.6g}, tof edges {tof[0]}..{tof[-1]}')
    if args.compare:
        ref = vol[args.chan] if args.chan >= 0 else (vol.sum(0) if args.chan == -1 else vol.max(0))
        _, img = read_img(args.compare)
        diff = np.abs(img.astype(np.float64) - ref.astype(np.float64))
        print(f'compare {args.compare}: max|diff| = {diff.max():.6g}, image sum = {img.sum():.6g}, ref sum = {ref.sum():.6g}')
        return 0 if diff.max() <= 1e-3 * max(1.0, np.abs(ref).max()) else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
