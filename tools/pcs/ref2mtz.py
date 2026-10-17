#!/usr/bin/env python3
"""Convert d*TREK TOF-Laue reflection lists to an unmerged MTZ file with a
LAMBDA column (for AIMLESS, or for comparison with the LAUENORM/SCALA
output of the original PCS workflow).

Each input list is one batch, unless it has an sBatch field (dtlauenorm
output), whose values are then used as batch numbers.

  ref2mtz.py -o out.mtz [--cell a b c al be ga] [--spacegroup P21]
             [--head file.head] [--wavefield fCalwavelength] list1.ref ...

The cell and space group are taken from --cell/--spacegroup, else from
--head, else from the header embedded in the first list (lists written by
this d*TREK).  Needs gemmi and numpy.
"""

import argparse
import re
import sys

import gemmi
import numpy as np


def read_ref(path):
    """Return (header_dict, names, rows) of a d*TREK .ref file."""
    with open(path) as f:
        lines = f.read().split('\n')
    counts = [int(x) for x in lines[0].split()]
    n_head = counts[3] if len(counts) > 3 else 0
    nf = sum(counts[:3])
    header = {}
    for line in lines[1:1 + n_head]:
        m = re.match(r'\s*([A-Za-z0-9_]+)=(.*?);?\s*$', line)
        if m:
            header[m.group(1)] = m.group(2).rstrip(';').split()
    names = [l.strip() for l in lines[1 + n_head:1 + n_head + nf]]
    rows = [l.split() for l in lines[1 + n_head + nf:] if l.strip()]
    rows = [r for r in rows if len(r) >= nf]
    return header, names, rows


def read_head(path):
    header = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'\s*([A-Za-z0-9_]+)=(.*?);?\s*$', line)
            if m:
                header[m.group(1)] = m.group(2).rstrip(';').split()
    return header


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('lists', nargs='+')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--cell', nargs=6, type=float)
    ap.add_argument('--spacegroup')
    ap.add_argument('--head')
    ap.add_argument('--wavefield', default='fCalwavelength')
    ap.add_argument('--keep-negative-sigma', action='store_true',
                    help='keep reflections with sigma <= 0 (default: drop)')
    args = ap.parse_args()

    hdr0, _, _ = read_ref(args.lists[0])
    hdr = read_head(args.head) if args.head else hdr0
    cell = args.cell or [float(x) for x in hdr.get('CRYSTAL_UNIT_CELL', [])[:6]]
    sg = args.spacegroup or (hdr['CRYSTAL_SPACEGROUP'][0] if 'CRYSTAL_SPACEGROUP' in hdr else None)
    if len(cell) != 6 or sg is None:
        sys.exit('ERROR: give --cell and --spacegroup (or --head)')
    sg = gemmi.find_spacegroup_by_number(int(sg)) if sg.isdigit() else gemmi.SpaceGroup(sg)
    ops = sg.operations()
    asu = gemmi.ReciprocalAsu(sg)

    hkl, isym, batch, inten, sigma, lam = [], [], [], [], [], []
    for n, path in enumerate(args.lists, start=1):
        _, names, rows = read_ref(path)
        idx = {k: names.index(k) for k in ('nH', 'nK', 'nL', 'fIntensity', 'fSigmaI')}
        if args.wavefield not in names:
            sys.exit(f'ERROR: {path} has no field {args.wavefield}')
        iw = names.index(args.wavefield)
        ib = names.index('sBatch') if 'sBatch' in names else None
        for r in rows:
            s = float(r[idx['fSigmaI']])
            if s <= 0.0 and not args.keep_negative_sigma:
                continue
            h = [int(r[idx['nH']]), int(r[idx['nK']]), int(r[idx['nL']])]
            if sg.operations().is_systematically_absent(h):
                continue
            ha, isy = asu.to_asu(h, ops)
            hkl.append(ha)
            isym.append(isy)
            if ib is not None:
                try:
                    batch.append(int(r[ib]))
                except ValueError:
                    batch.append(n)
            else:
                batch.append(n)
            inten.append(float(r[idx['fIntensity']]))
            sigma.append(s)
            lam.append(float(r[iw]))

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = sg
    mtz.set_cell_for_all(gemmi.UnitCell(*cell))
    mtz.add_dataset('TOFLaue')
    for name, typ in (('M/ISYM', 'Y'), ('BATCH', 'B'), ('I', 'J'), ('SIGI', 'Q'),
                      ('LAMBDA', 'R')):
        mtz.add_column(name, typ)
    data = np.column_stack([np.array(hkl, dtype=float), np.array(isym, dtype=float),
                            np.array(batch, dtype=float), np.array(inten), np.array(sigma),
                            np.array(lam)])
    mtz.set_data(data)
    for b in sorted(set(batch)):
        bh = gemmi.Mtz.Batch()
        bh.number = b
        bh.cell = gemmi.UnitCell(*cell)
        bh.dataset_id = 1
        bh.wavelength = float(np.median(np.array(lam)[np.array(batch) == b]))
        mtz.batches.append(bh)
    mtz.history = ['ref2mtz.py: unmerged TOF-Laue data from d*TREK lists']
    mtz.sort()
    mtz.write_to_file(args.out)
    print(f'{args.out}: {len(inten)} reflections, {len(set(batch))} batches, '
          f'space group {sg.hm}, cell {cell}')


if __name__ == '__main__':
    main()
