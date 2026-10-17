#!/bin/sh
#
# lauetof_lauenorm.sh: wavelength-normalize and merge the integrated lists of
# a TOF-Laue data set with LAUENORM, using a LAUENORM control file such as
# the lnorm.dat that Mantid's mandi_createmtz.py writes.
#
# Usage: lauetof_lauenorm.sh [options] control.dat run1i.ref run2i.ref ...
#
#   control.dat  LAUENORM control cards.  The number of film packs of card 3
#                ("NORMALISE n") is set to the number of lists.
#
# Options:
#   -o dir        output directory (default: current directory)
#   -setup file   shell script that puts lauenorm in the PATH and sets up
#                 CCP4 (e.g. Daresbury_laue/laue-setup.sh)
#   -pgdata file  LAUENORM point group data (default: lib/dat/pglib.dat
#                 beside the bin directory of lauenorm)
#   -out name     output MTZ file name (default lauenorm.mtz)
#
# The lists are converted with "dtref2afout -scale auto -theta" to
# laueNorm001, laueNorm002, ... (one film pack per run).
# Output (in dir): the MTZ file and lauenorm.log.

BIN=${DTREK_BIN:-$(cd "$(dirname "$0")/../DTTREK" && pwd)}
OUT=.
SETUP=""
PGDATA_IN=""
MTZ=lauenorm.mtz

while [ $# -gt 0 ]; do
    case "$1" in
        -o)       OUT=$2; shift 2 ;;
        -setup)   SETUP=$2; shift 2 ;;
        -pgdata)  PGDATA_IN=$2; shift 2 ;;
        -out)     MTZ=$2; shift 2 ;;
        -*)       echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)        break ;;
    esac
done
if [ $# -lt 2 ]; then
    sed -n '3,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
if [ -n "$SETUP" ]; then
    . "$SETUP" > /dev/null 2>&1
fi
command -v lauenorm > /dev/null || { echo "ERROR: lauenorm is not in the PATH (use -setup)" >&2; exit 1; }
[ -n "$CCP4" ] || { echo "ERROR: CCP4 is not set up (use -setup)" >&2; exit 1; }

CTRL=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); shift
mkdir -p "$OUT" || exit 1
LISTS=""
for f in "$@"; do
    LISTS="$LISTS $(cd "$(dirname "$f")" && pwd)/$(basename "$f")"
done
cd "$OUT" || exit 1

n=0
for f in $LISTS; do
    n=$((n + 1))
    pack=$(printf '%03d' $n)
    "$BIN/dtref2afout" -scale auto -theta -out "laueNorm$pack" "$f" \
        || { echo "ERROR: dtref2afout failed for $f" >&2; exit 1; }
    eval "LAUE$pack=\$PWD/laueNorm$pack; export LAUE$pack"
done
sed "s/^NORMALISE *[0-9]*/NORMALISE $n/" "$CTRL" > lauenorm.dat

LBIN=$(dirname "$(command -v lauenorm)")
PGDATA=${PGDATA_IN:-$(cd "$LBIN/../.." 2>/dev/null && pwd)/lib/dat/pglib.dat}
[ -f "$PGDATA" ] || { echo "ERROR: no $PGDATA (use -pgdata)" >&2; exit 1; }
SYMOP=${SYMOP:-$CCP4/lib/data/symop.lib}
SYMINFO=${SYMINFO:-$CCP4/lib/data/syminfo.lib}
HKLOUT=$PWD/$MTZ
HKLMULT=$PWD/hklmult.out
MULTDIAG=$PWD/multdiag.out
export PGDATA SYMOP SYMINFO HKLOUT HKLMULT MULTDIAG

echo "lauenorm ($n runs)"
lauenorm < lauenorm.dat > lauenorm.log 2>&1 \
    || { echo "ERROR: lauenorm failed, see $OUT/lauenorm.log" >&2; exit 1; }
grep -A12 "LAUE MERGING R-FACTORS" lauenorm.log | tail -12 | grep "RFACT\|N =" | sed 's/^ */  /'
grep "NUMBER OF REFLECTIONS WRITTEN" lauenorm.log | sed 's/^ */  /'
echo "Done: $OUT/$MTZ"
