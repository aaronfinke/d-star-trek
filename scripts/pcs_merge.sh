#!/bin/sh
#
# pcs_merge.sh: wavelength-normalize, scale and merge the integrated lists
# of the runs of one LANSCE PCS TOF-Laue data set (the LAUENORM + SCALA step
# of the original PCS workflow):
#
#   dtlauenorm -> dtscaleaverage -> tools/pcs/ref2mtz.py
#
# Usage: pcs_merge.sh -head crystal.head [options] run1i.ref run2i.ref ...
#
#   -head file       header with the cell and space group (e.g. a refined
#                    header of one of the runs)
#   -o dir           output directory (default: current directory)
#   -wave "min max"  wavelength range for the normalization (default: all)
#   -reso "lo hi"    resolution range for dtscaleaverage (default "30 1.1")
#   -knots n         knots of the wavelength curve (default 10)
#
# Output (in dir): dtlauenorm.ref (normalized, unmerged; sBatch = run
# number), dtlauenorm_curve.txt, dtscaleaverage.log (merging statistics),
# scaled_unavg.ref and unmerged.mtz (with a LAMBDA column, for AIMLESS).
#
# To use LAUENORM instead, convert the lists with
#   dtref2afout -scale 1000 run1i.ref ...
# and run LAUENORM with a laue.dat as in the original workflow; with the
# sigmas of this d*TREK, use an I/sigma limit for the scale factors of
# about 2 instead of the 7 of the legacy laue.dat.

BIN=${DTREK_BIN:-$(cd "$(dirname "$0")/../DTTREK" && pwd)}
TOOLS=$(cd "$(dirname "$0")/../tools/pcs" && pwd)
OUT=.
HEAD=""
WAVE=""
RESO="30 1.1"
KNOTS=10

while [ $# -gt 0 ]; do
    case "$1" in
        -head)   HEAD=$2; shift 2 ;;
        -o)      OUT=$2; shift 2 ;;
        -wave)   WAVE=$2; shift 2 ;;
        -reso)   RESO=$2; shift 2 ;;
        -knots)  KNOTS=$2; shift 2 ;;
        -*)      echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)       break ;;
    esac
done
if [ -z "$HEAD" ] || [ $# -lt 2 ]; then
    sed -n '3,28p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

mkdir -p "$OUT" || exit 1
HEAD=$(cd "$(dirname "$HEAD")" && pwd)/$(basename "$HEAD")
LISTS=""
for f in "$@"; do
    LISTS="$LISTS $(cd "$(dirname "$f")" && pwd)/$(basename "$f")"
done
cd "$OUT" || exit 1

WOPT=""
[ -n "$WAVE" ] && WOPT="-wave $WAVE"

echo "dtlauenorm ($# runs)"
"$BIN/dtlauenorm" -head "$HEAD" $WOPT -knots $KNOTS $LISTS > dtlauenorm.log 2>&1 \
    || { echo "ERROR: dtlauenorm failed, see $OUT/dtlauenorm.log" >&2; exit 1; }
grep -E "Normalization curve \(|^  All" dtlauenorm.log | sed 's/^ */  /'

echo "dtscaleaverage"
"$BIN/dtscaleaverage" "$HEAD" dtlauenorm.ref -reso $RESO -batchscale \
    -ref scaled_unavg.ref > dtscaleaverage.log 2>&1 \
    || { echo "ERROR: dtscaleaverage failed, see $OUT/dtscaleaverage.log" >&2; exit 1; }
grep -E "Number of unique|% completeness|Rmerge  |Output <I/sigI>|Average redundancy" \
    dtscaleaverage.log | tail -5 | sed 's/^ */  /'

echo "ref2mtz.py"
python3 -I "$TOOLS/ref2mtz.py" -o unmerged.mtz --head "$HEAD" dtlauenorm.ref \
    || echo "WARNING: ref2mtz.py failed (needs python3 with gemmi and numpy)" >&2
echo "Done: $OUT"
