#!/bin/sh
#
# pcs_process_run.sh: process one LANSCE PCS TOF-Laue run with d*TREK, as the
# laue_process.txt recipe of the original PCS workflow did:
#
#   dtfind -> dtindex -> dtrefine (3 stages) -> dtintegrate -> dtref2afout
#
# Usage: pcs_process_run.sh [options] start.head run.(h5|nx.hdf) flat.(h5|nx.hdf)
#
#   start.head  header with the instrument calibration (L1, delay, beam
#               centre, detector goniometer), the cell, space group and
#               crystal goniometer setting of the run (e.g. the
#               dtextractheader output edited as in laue_process.txt, or a
#               refined header of a previous run with CRYSTAL_GONIO_VALUES
#               changed).
#   run         the PCS NeXus file of the run (HDF4 files are converted with
#               pcs_h4toh5.sh first).
#   flat        the V/Nb run used for the detector response correction.
#
# Options:
#   -o dir            output directory (default: current directory)
#   -sg N             space group (default: CRYSTAL_SPACEGROUP of start.head)
#   -maxcell f        dtindex -maxcell (default 50)
#   -ireso "lo hi"    dtindex resolution range (default "4 1")
#   -ithresh f        dtindex -threshold (default 100)
#   -rthresh f        dtrefine -threshold (default 10)
#   -reso "hi lo"     dtintegrate resolution range (default "1.0 20")
#   -wave "max min"   wavelength range (default "7.0 0.6")
#   -noindex          skip dtindex (use the orientation of start.head)
#
# The intensity thresholds are lower than in the legacy recipes (1000-4000
# for dtindex, 100 for dtrefine): the intensities of this d*TREK are not
# scaled up as those of the original program were.
#
# Output (in dir): <run>find.ref, <run>refine.head, <run>i.ref (the
# integrated list) and <run>i.afout (LAUENORM input), plus the logs.

BIN=${DTREK_BIN:-$(cd "$(dirname "$0")/../DTTREK" && pwd)}
OUT=.
SG=""
MAXCELL=50
IRESO="4 1"
ITHRESH=100
RTHRESH=10
RESO="1.0 20"
WAVE="7.0 0.6"
INDEX=1

while [ $# -gt 0 ]; do
    case "$1" in
        -o)        OUT=$2; shift 2 ;;
        -sg)       SG=$2; shift 2 ;;
        -maxcell)  MAXCELL=$2; shift 2 ;;
        -ireso)    IRESO=$2; shift 2 ;;
        -ithresh)  ITHRESH=$2; shift 2 ;;
        -rthresh)  RTHRESH=$2; shift 2 ;;
        -reso)     RESO=$2; shift 2 ;;
        -wave)     WAVE=$2; shift 2 ;;
        -noindex)  INDEX=0; shift ;;
        -*)        echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)         break ;;
    esac
done
if [ $# -ne 3 ]; then
    sed -n '3,40p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
HEAD=$1
DATA=$2
FLAT=$3

mkdir -p "$OUT" || exit 1
OUT=$(cd "$OUT" && pwd)

# HDF4 -> HDF5 (into the output directory)
to_h5() {
    case "$1" in
        *.h5) echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")" ;;
        *)    sh "$(dirname "$0")/pcs_h4toh5.sh" -o "$OUT" "$1" > /dev/null || return 1
              b=$(basename "$1"); b=${b%.gz}; b=${b%.hdf}; b=${b%.nx}
              echo "$OUT/$b.h5" ;;
    esac
}
DATA5=$(to_h5 "$DATA") || { echo "ERROR converting $DATA" >&2; exit 1; }
FLAT5=$(to_h5 "$FLAT") || { echo "ERROR converting $FLAT" >&2; exit 1; }
RUN=$(basename "$DATA5" .h5)

[ -n "$SG" ] || SG=$(tr -d '\r' < "$HEAD" | sed -n 's/^CRYSTAL_SPACEGROUP=\([0-9]*\).*/\1/p' | head -1)
[ -n "$SG" ] || { echo "ERROR: no space group (use -sg)" >&2; exit 1; }

cd "$OUT" || exit 1
tr -d '\r' < "$HEAD" | sed "s#^SCAN_TEMPLATE=.*#SCAN_TEMPLATE=$DATA5;#" > "$RUN.head"
grep -q '^SCAN_TEMPLATE=' "$RUN.head" || echo "SCAN_TEMPLATE=$DATA5;" >> "$RUN.head"

run() {
    log=$1; shift
    echo "  $*"
    "$@" > "$log" 2>&1 || { echo "ERROR: $1 failed, see $OUT/$log" >&2; exit 1; }
}

echo "Processing $RUN in $OUT"
run "${RUN}_1find.log" "$BIN/dtfind" "$RUN.head" -cor "$FLAT5" -type 3 -seq 0 95 \
    -pad 4 -window 20 20 -sigma 4 -min 5 -rect 2 12 253 1910 \
    -out "${RUN}find.head" -ref "${RUN}find.ref"

START="${RUN}find.head"
if [ $INDEX -eq 1 ]; then
    run "${RUN}_2index.log" "$BIN/dtindex" "${RUN}find.head" "${RUN}find.ref" -reso $IRESO \
        -known -threshold $ITHRESH -spacegroup $SG -maxcell $MAXCELL -fourier \
        -errors 5 5 5 3 3 3 -nodeice -out "${RUN}index.head"
    START="${RUN}index.head"
fi

W1=${WAVE% *}; W2=${WAVE#* }
run "${RUN}_3refine1.log" "$BIN/dtrefine" "$START" "${RUN}find.ref" +CrysRot \
    -threshold $RTHRESH -cycles 40 -reso 100 1.0 -sigma 3 -wave $W1 $W2 -rej 50 50 50 \
    -go -out "${RUN}r1.head"
run "${RUN}_4refine2.log" "$BIN/dtrefine" "${RUN}r1.head" "${RUN}find.ref" +CrysRot \
    +DetTrans1 +DetTrans2 -threshold $RTHRESH -cycles 40 -reso 100 1.0 -sigma 3 \
    -wave $W1 $W2 -rej 1.5 1.5 1 -go -out "${RUN}r2.head"
run "${RUN}_5refine3.log" "$BIN/dtrefine" "${RUN}r2.head" "${RUN}find.ref" +CrysRot +DetAll \
    -threshold $RTHRESH -cycles 40 -reso 100 1.0 -sigma 3 -wave $W1 $W2 -rej 1 1 1 \
    -go -out "${RUN}refine.head"

run "${RUN}_6integrate.log" "$BIN/dtintegrate" "${RUN}refine.head" -flat -cor "$FLAT5" \
    -type 204 -seq 0 95 -reso $RESO -minpeakrad 2.0 -h 100 -100 -k 100 -100 -l 100 -100 \
    -rect 0 18 255 1910 -wave $W1 $W2 -pad 4 -window 20 20 -ref "${RUN}i.ref"

run "${RUN}_7afout.log" "$BIN/dtref2afout" -scale 1000 -out "${RUN}i.afout" "${RUN}i.ref"

grep -h "rmsResid" "${RUN}_5refine3.log" | sed 's/^/  /'
grep -h "integrated (" "${RUN}_6integrate.log" | sed 's/^ */  /'
echo "Done: $OUT/${RUN}i.ref, $OUT/${RUN}i.afout"
