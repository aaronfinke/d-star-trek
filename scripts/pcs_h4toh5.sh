#!/bin/sh
#
# pcs_h4toh5.sh: convert LANSCE PCS HDF4 NeXus files (*.nx.hdf, *.nx.hdf.gz)
# to HDF5 so that d*TREK can read them.  The HDF5 file keeps the HDF4 group
# and dataset names (/run_N/tof_detector_K/data, ...).
#
# Usage: pcs_h4toh5.sh [-o outdir] file.nx.hdf [file.nx.hdf ...]
#
# Output: outdir/<name>.h5 (default outdir: the directory of each input file).
# Existing outputs are not overwritten.  Set H4TOH5 to the h4toh5 executable
# if it is not on the PATH.

H4TOH5=${H4TOH5:-h4toh5}
OUTDIR=""

if [ "$1" = "-o" ]; then
    OUTDIR=$2
    shift 2
fi
if [ $# -eq 0 ]; then
    echo "Usage: $0 [-o outdir] file.nx.hdf [file.nx.hdf ...]" >&2
    exit 1
fi
if ! command -v "$H4TOH5" > /dev/null 2>&1; then
    echo "ERROR: $H4TOH5 not found; set H4TOH5 to the h4toh5 executable." >&2
    exit 1
fi

status=0
for in in "$@"; do
    dir=${OUTDIR:-$(dirname "$in")}
    base=$(basename "$in")
    base=${base%.gz}
    base=${base%.hdf}
    base=${base%.nx}
    out="$dir/$base.h5"
    if [ -e "$out" ]; then
        echo "$out exists, skipping"
        continue
    fi
    src=$in
    tmp=""
    case "$in" in
        *.gz)
            tmp=$(mktemp "${TMPDIR:-/tmp}/pcs_h4toh5.XXXXXX") || exit 1
            gunzip -c "$in" > "$tmp" || { rm -f "$tmp"; status=1; continue; }
            src=$tmp
            ;;
    esac
    if "$H4TOH5" "$src" "$out"; then
        echo "$in -> $out"
    else
        echo "ERROR converting $in" >&2
        rm -f "$out"
        status=1
    fi
    [ -n "$tmp" ] && rm -f "$tmp"
done
exit $status
