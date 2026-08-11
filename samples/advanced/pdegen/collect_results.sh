#!/bin/bash
# ============================================================================
# collect_results.sh — snapshot a result set together with its provenance
#
# A directory of CSVs is not a dataset: in a month nobody remembers which code,
# which library and which environment produced it, and the run gets repeated.
# This records all of that next to the numbers, so a result set can be trusted
# without re-running it.
#
# Score-P experiment directories are folded into the result set too. Their
# default names carry a timestamp and nothing else, so on their own nobody can
# tell which scale point a profile belongs to; here they are moved inside the
# result directory and listed with their times, in run order.
#
# USE (on the cluster, from samples/advanced/pdegen):
#     ./collect_results.sh results_amg_comm_44449535 [label] [scorep-dir ...]
#
# With no scorep-dir given, every scorep-* in the current directory is taken.
#
# Produces  <resdir>/PROVENANCE.md  and a tarball ready to transfer.
# ============================================================================
set -u

RESDIR=${1:?usage: collect_results.sh <results-dir> [label] [scorep-dir ...]}
LABEL=${2:-}
shift 2 2>/dev/null || shift $#
SCOREP_DIRS=("$@")
if [ ${#SCOREP_DIRS[@]} -eq 0 ]; then
    shopt -s nullglob
    SCOREP_DIRS=(scorep-*)
    shopt -u nullglob
fi

if [ ! -d "$RESDIR" ]; then
    echo "No such directory: $RESDIR" >&2
    exit 1
fi

JOBID=$(basename "$RESDIR" | sed 's/.*_//')
OUT=$RESDIR/PROVENANCE.md

# Fold the Score-P experiments into the result set, oldest first: the sbatch
# runs the scale points in order, so run order is the only handle available to
# match a timestamped profile to its rank count.
PROFDIR=$RESDIR/scorep
if [ ${#SCOREP_DIRS[@]} -gt 0 ]; then
    mkdir -p "$PROFDIR"
    echo "Folding ${#SCOREP_DIRS[@]} Score-P experiment(s) into $PROFDIR:"
    while IFS= read -r d; do
        [ -d "$d" ] || continue
        printf '  %s  %s\n' "$(date -r "$d" '+%Y-%m-%d %H:%M:%S')" "$d"
        mv "$d" "$PROFDIR/"
    done < <(ls -1dtr "${SCOREP_DIRS[@]}" 2>/dev/null)
fi

{
    echo "# Result set $JOBID"
    [ -n "$LABEL" ] && echo -e "\n**$LABEL**"
    echo
    echo "Collected: $(date -Iseconds) on $(hostname)"
    echo

    echo "## Instrumentation"
    echo
    # The single most misread property of a result set: whether the binary was
    # instrumented and whether profiling was actually on during the run.
    if [ -d "$PROFDIR" ]; then
        echo "Score-P experiments were collected, so profiling was ACTIVE:"
        echo "absolute timings in the CSVs include instrumentation overhead and are"
        echo "comparable between schemes, not against a clean build."
        echo
        echo "Profiles in run order (the sbatch runs the scale points in order,"
        echo "so the n-th profile corresponds to the n-th entry of RANK_POINTS;"
        echo "check against the timestamps in the per-step run.out to be sure):"
        echo
        echo '```'
        ls -1dtr "$PROFDIR"/*/ 2>/dev/null | while read -r d; do
            printf '%s  %8s  %s\n' "$(date -r "$d" '+%Y-%m-%d %H:%M:%S')" \
                "$(du -sh "$d" | cut -f1)" "$(basename "$d")"
        done
        echo '```'
        echo
        echo "Read them with:"
        echo '```'
        echo "cube_stat -p <profile>/profile.cubex"
        echo "cube_dump -m time -c all <profile>/profile.cubex"
        echo '```'
    else
        echo "No Score-P experiments collected; profiling was probably off."
    fi
    echo
    echo '```'
    echo "SCOREP_ENABLE_PROFILING=${SCOREP_ENABLE_PROFILING:-<unset, defaults to true>}"
    echo "SCOREP_ENABLE_TRACING=${SCOREP_ENABLE_TRACING:-<unset, defaults to false>}"
    echo '```'
    echo
    echo "## MPI runtime settings"
    echo
    echo "These change which implementation the schemes actually run on, so they"
    echo "have to be held constant across a campaign. HCOLL in particular"
    echo "accelerates collectives, which the two neighborhood schemes use and the"
    echo "point-to-point and one-sided ones do not: switching it changes what is"
    echo "being compared, not just how fast it is."
    echo
    echo '```'
    env | grep -E "^(OMPI_MCA_|UCX_|SLURM_CPU_BIND)" | sort || echo "(none set)"
    echo '```'
    echo

    echo "## Code"
    echo
    echo '```'
    for repo in "$HOME/Desktop/scorep/amg4psblas" "$HOME/Desktop/scorep/psblas3"; do
        if [ -d "$repo/.git" ]; then
            printf '%-12s %s  %s\n' "$(basename "$repo")" \
                "$(git -C "$repo" rev-parse --short HEAD 2>/dev/null)" \
                "$(git -C "$repo" rev-parse --abbrev-ref HEAD 2>/dev/null)"
            if [ -n "$(git -C "$repo" status --porcelain 2>/dev/null)" ]; then
                echo "             UNCOMMITTED CHANGES PRESENT"
            fi
        else
            printf '%-12s not a git tree\n' "$(basename "$repo")"
        fi
    done
    echo '```'
    echo

    echo "## Build configuration"
    echo
    echo '```'
    grep -E "^PSBLASDIR|^FC=|^CC=|^FCOPT" "$HOME/Desktop/scorep/amg4psblas/Make.inc" 2>/dev/null
    grep -E "^FCUDEFINES|^CUDA_DIR" \
        "$HOME/opt/psblas-cuda/include/Make.inc.psblas" 2>/dev/null
    echo '```'
    echo

    echo "## Job"
    echo
    echo '```'
    sacct -j "$JOBID" --format=JobID,JobName,Partition,QOS,NNodes,NTasks,Elapsed,State 2>/dev/null \
        || echo "sacct unavailable"
    echo '```'
    echo

    echo "## Parameters"
    echo
    echo '```'
    grep -E "^(DIM|NREP|NLEV|ITMAX|RANK_POINTS)=" amg_comm_strong.sh 2>/dev/null
    echo '```'
    echo

    echo "## Contents"
    echo
    echo '```'
    find "$RESDIR" -type f -printf "%10s  %p\n" | sort -k2
    echo '```'
} > "$OUT"

# Keep the exact script that produced the numbers: parameters drift.
cp -p amg_comm_strong.sh "$RESDIR/amg_comm_strong.sh.used" 2>/dev/null || true

TARBALL=${RESDIR%/}.tar.gz
tar czf "$TARBALL" "$RESDIR"

echo "Provenance written to $OUT"
echo "Tarball: $TARBALL  ($(du -h "$TARBALL" | cut -f1))"
echo
echo "Transfer with:"
echo "  rsync -avz <user>@<host>:$(readlink -f "$TARBALL") ."
