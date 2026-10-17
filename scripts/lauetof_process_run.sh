#!/bin/sh
#
# lauetof_process_run.sh: process one neutron TOF-Laue run of a panel
# instrument (e.g. SNS MANDI, ESS NMX) with d*TREK:
#
#   [dteventbin] -> dtfind -> dtindex -> dtrefine (3 stages) -> dtintegrate
#   -> dtref2afout -> dtexport
#
# Usage: lauetof_process_run.sh [options] -cell "a b c al be ga" -sg N data.h5
#
#   data.h5     a NeXus event file (histogrammed with dteventbin first) or
#               an NXlauetof histogram file (e.g. dteventbin output or the
#               ESS NMX reduced data).
#
# Options:
#   -o dir             output directory (default: current directory)
#   -cell "a b ..."    known unit cell (required)
#   -sg N              space group number (required)
#   -detcal file       dteventbin -detcal (ISAW/Mantid .DetCal calibration)
#   -goniocal file     dteventbin -goniocal (recorded only)
#   -bin "options"     further dteventbin options (e.g. "-logbin 0.005")
#   -reso "lo hi"      resolution range in A (default "15 2.5")
#   -wave "max min"    wavelength range in A (default "4.0 2.0")
#   -fsigma f          dtfind -sigma (default 3)
#   -maxcell f         dtindex -maxcell (default 1.05 x the longest axis)
#   -tofradius n       dtintegrate -tofradius (default 5)
#   -batch n           MTZ batch number (default: the digits of the run name,
#                      e.g. 13378 for MANDI_13378, so that runs do not clash)
#
# The last refinement stage refines an unconstrained (triclinic) cell, as
# Mantid does for its UB matrix: small systematic errors of the panel
# geometry are absorbed by the cell, which places the predictions of the
# outer panels much better.  The output header keeps the space group.
#
# Output (in dir): find.ref, refine.head, <run>i.ref (the integrated list),
# <run>i.afout (LAUENORM input with theta) and <run>i.mtz (unmerged MTZ for
# AIMLESS, one batch with the crystal orientation; written if dtexport was
# built, which needs gemmi), plus the logs.  To put several runs in one MTZ
# file for AIMLESS (scale with SCALES CONSTANT, see TOF_LAUE.md):
#
#   dtexport -head r1/r2.head -out all.mtz r*/refine.head r*/*i.afout

BIN=${DTREK_BIN:-$(cd "$(dirname "$0")/../DTTREK" && pwd)}
OUT=.
CELL=""
SG=""
DETCAL=""
GONIOCAL=""
BINOPT=""
RESO="15 2.5"
WAVE="4.0 2.0"
FSIGMA=3
MAXCELL=""
TOFRAD=5
BATCH=""

while [ $# -gt 0 ]; do
    case "$1" in
        -o)         OUT=$2; shift 2 ;;
        -cell)      CELL=$2; shift 2 ;;
        -sg)        SG=$2; shift 2 ;;
        -detcal)    DETCAL=$2; shift 2 ;;
        -goniocal)  GONIOCAL=$2; shift 2 ;;
        -bin)       BINOPT=$2; shift 2 ;;
        -reso)      RESO=$2; shift 2 ;;
        -wave)      WAVE=$2; shift 2 ;;
        -fsigma)    FSIGMA=$2; shift 2 ;;
        -maxcell)   MAXCELL=$2; shift 2 ;;
        -tofradius) TOFRAD=$2; shift 2 ;;
        -batch)     BATCH=$2; shift 2 ;;
        -*)         echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)          break ;;
    esac
done
if [ $# -ne 1 ] || [ -z "$SG" ] || [ -z "$CELL" ]; then
    sed -n '3,41p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
DATA=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
[ -f "$DATA" ] || { echo "ERROR: no file $1" >&2; exit 1; }
abspath() { [ -n "$1" ] && echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
DETCAL=$(abspath "$DETCAL")
GONIOCAL=$(abspath "$GONIOCAL")

mkdir -p "$OUT" || exit 1
OUT=$(cd "$OUT" && pwd)
cd "$OUT" || exit 1

run() {
    log=$1; shift
    echo "  $*"
    "$@" > "$log" 2>&1 || { echo "ERROR: $(basename "$1") failed, see $OUT/$log" >&2; exit 1; }
}

RUN=$(basename "$DATA"); RUN=${RUN%.h5}; RUN=${RUN%.nxs}; RUN=${RUN%.hdf5}; RUN=${RUN%_lauetof}
echo "Processing $RUN in $OUT"

# Event data (NXevent_data groups <bank>_events) -> NXlauetof histograms
H5=$DATA
if h5ls "$DATA/entry" 2>/dev/null | grep -q '_events  *Group'; then
    H5=$OUT/${RUN}_lauetof.h5
    set --
    [ -n "$DETCAL" ] && set -- "$@" -detcal "$DETCAL"
    [ -n "$GONIOCAL" ] && set -- "$@" -goniocal "$GONIOCAL"
    run "${RUN}_0bin.log" "$BIN/dteventbin" "$DATA" -out "$H5" "$@" $BINOPT
fi

R1=${RESO% *}; R2=${RESO#* }
W1=${WAVE% *}; W2=${WAVE#* }
# Size of the panel mosaic and number of TOF channels, from the header
# d*TREK builds for the file
run "${RUN}_0sum.log" "$BIN/dtpcsproject" "$H5" "${RUN}_sum.img" -sum
DIM=$(head -c 65536 "${RUN}_sum.img" | LC_ALL=C sed -n 's/^DETECTOR_DIMENSIONS=\([^;]*\);.*/\1/p' | head -1)
NX=$(echo $DIM | awk '{print $1}'); NY=$(echo $DIM | awk '{print $2}')
NCH=$(echo $DIM | awk '{print $3}')
[ -n "$NX" ] && [ -n "$NY" ] && [ -n "$NCH" ] \
    || { echo "ERROR: cannot read the dimensions of $H5" >&2; exit 1; }
LAST=$((NCH - 1))
run "${RUN}_1find.log" "$BIN/dtfind" "$H5" -seq 0 $LAST -pad 4 -window 20 20 \
    -sigma $FSIGMA -min 2 -rect 0 0 $NX $NY -out find.head -ref find.ref
[ -n "$MAXCELL" ] || MAXCELL=$(echo $CELL | awk '{m=$1; if ($2>m) m=$2; if ($3>m) m=$3; print int(1.05*m+1)}')
run "${RUN}_2index.log" "$BIN/dtindex" find.head find.ref -cell $CELL -spacegroup $SG \
    -known -fourier -reso $R1 $R2 -maxcell $MAXCELL -nodeice -errors 5 5 5 3 3 3 \
    -out index.head
START=index.head

run "${RUN}_3refine1.log" "$BIN/dtrefine" "$START" find.ref +CrysRot -sigma 3 \
    -wave $W1 $W2 -reso $R1 $R2 -rej 50 50 50 -cycles 50 -go -rej 3 3 2 -go -out r1.head
run "${RUN}_4refine2.log" "$BIN/dtrefine" r1.head find.ref +CrysAll -sigma 3 \
    -wave $W1 $W2 -reso $R1 $R2 -rej 3 3 2 -cycles 50 -go -out r2.head
run "${RUN}_5refine3.log" "$BIN/dtrefine" r2.head find.ref +CrysAll -CrysConstraints -sigma 3 \
    -wave $W1 $W2 -reso $R1 $R2 -rej 3 3 2 -cycles 50 -go -go -out refine.head

run "${RUN}_6integrate.log" "$BIN/dtintegrate" refine.head -seq 0 $LAST -reso $R1 $R2 \
    -wave $W1 $W2 -tofradius $TOFRAD -ref "${RUN}i.ref"
run "${RUN}_7afout.log" "$BIN/dtref2afout" -scale auto -theta -keepneg -out "${RUN}i.afout" "${RUN}i.ref"

# Unmerged MTZ for AIMLESS, with the orientation of refine.head and the
# space-group-constrained cell of r2.head
MTZ=""
if [ -x "$BIN/dtexport" ]; then
    [ -n "$BATCH" ] || BATCH=$(echo "$RUN" | sed -n 's/^\(.*[^0-9]\)\{0,1\}\([0-9][0-9]*\)$/\2/p' | sed 's/^0*\([0-9]\)/\1/')
    [ -n "$BATCH" ] && [ "$BATCH" -gt 0 ] || BATCH=1
    run "${RUN}_8mtz.log" "$BIN/dtexport" -head r2.head -batch $BATCH -out "${RUN}i.mtz" \
        refine.head "${RUN}i.afout"
    MTZ=", $OUT/${RUN}i.mtz (batch $BATCH)"
else
    echo "  dtexport is not built (it needs the gemmi library): no ${RUN}i.mtz"
fi

grep -h "rmsResid (mm)\|rmsResid (TOF" "${RUN}_5refine3.log" | tail -2 | sed 's/^/  /'
tr -d '\r' < refine.head | grep "^CRYSTAL_UNIT_CELL=" | sed 's/^/  /'
echo "Done: $OUT/${RUN}i.ref, $OUT/${RUN}i.afout$MTZ"
