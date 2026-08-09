#!/bin/bash
# ============================================================================
# amg_comm_strong.sbatch   STRONG scaling, communication schemes on AMG
#
# Fixed total problem size; MPI ranks grow across scale points. Every scale
# point runs inside ONE reserved allocation, so the queue is paid once.
#
# The executable sweeps, in a single invocation:
#   * uniform      the 5 schemes applied to the whole hierarchy
#   * sensitivity  baseline everywhere, one level at a time moved onto another
#                  scheme -- the marginal value of a scheme AT a level, which
#                  is what a per-level policy is built on
# so there is no per-scheme or per-level loop here, and no reason to come back
# for a second allocation.
#
# Both CSVs are written per scale point; nothing needs re-running to plot.
# ============================================================================
#SBATCH --job-name=amg_comm_strong
#SBATCH --output=amg_comm_strong_%j.out
#SBATCH --error=amg_comm_strong_%j.err
#SBATCH --nodes=8                   # reserve the maximum up front
#SBATCH --ntasks-per-node=80
#SBATCH --cpus-per-task=1           # CPU-only: one core per rank
#SBATCH --time=01:00:00
#SBATCH --qos=acc_debug
# Account intentionally not hardcoded. Pass it at submit time:
#   sbatch -A <account> amg_comm_strong.sh
# or export SBATCH_ACCOUNT=<account> once in your cluster profile.

# ============================================================================
# USER CONFIGURATION
# ============================================================================
DIM=200                             # FIXED problem size (idim^3 unknowns)
NREP=7                              # repetitions per configuration
NLEV=5                              # max levels requested of the hierarchy
ITMAX=500
RANK_POINTS="80 160 320 640"        # total ranks per scale point (multiples of 80)

EXE=$SLURM_SUBMIT_DIR/runs/amg_d_comm_test

# ============================================================================
# ENVIRONMENT
# ============================================================================
module purge
module load bsc/1.0
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc

export OMPI_MCA_coll_hcoll_enable=0

RESDIR=$SLURM_SUBMIT_DIR/results_amg_comm_${SLURM_JOB_ID}
mkdir -p $RESDIR

echo "=== AMG communication schemes, STRONG scaling (CPU-only) ==="
echo "  fixed_dim=$DIM nrep=$NREP max_levels=$NLEV itmax=$ITMAX"
echo "  reserved_nodes=$SLURM_NNODES rank_points=[$RANK_POINTS]"
echo "  exe=$EXE"
echo "============================================================"

if [ ! -x "$EXE" ]; then
    echo "FATAL: $EXE not found or not executable. Build it before submitting:"
    echo "  cd samples/advanced/pdegen && make amg_d_comm_test"
    exit 1
fi

FAILED=0

for NRANKS in $RANK_POINTS; do
    NNODES=$(( (NRANKS + 79) / 80 ))
    STEP_DIR=$RESDIR/${NRANKS}ranks
    mkdir -p $STEP_DIR

    echo ""
    echo ">>> STRONG point: $NRANKS ranks ($NNODES nodes), fixed dim=$DIM"

    srun -N $NNODES -n $NRANKS --ntasks-per-node=80 --cpus-per-task=1 \
        $EXE $DIM $NREP $NLEV $ITMAX \
        --mode=all --csv=$STEP_DIR/amg_comm.csv \
        > $STEP_DIR/run.out 2>&1
    RC=$?
    echo ">>> exit=$RC output=$STEP_DIR/run.out"

    # Fail loudly rather than silently producing a dataset nobody can trust:
    # if a scheme did not reach every level, the comparison is meaningless.
    if ! grep -q "SCHEME PROPAGATION: OK" $STEP_DIR/run.out; then
        echo "!!! WARNING: scheme propagation not confirmed at $NRANKS ranks"
        FAILED=1
    fi
    # All configurations must converge identically; a differing iteration count
    # means the halo exchange changed the arithmetic, not just its schedule.
    NITERS=$(grep -oE "it +[0-9]+" $STEP_DIR/run.out | awk '{print $2}' | sort -u | wc -l)
    if [ "$NITERS" != "1" ]; then
        echo "!!! WARNING: iteration count is not constant at $NRANKS ranks ($NITERS distinct values)"
        FAILED=1
    fi
done

echo ""
if [ "$FAILED" = "0" ]; then
    echo "=== AMG COMM STRONG DONE, all checks passed. Results: $RESDIR ==="
else
    echo "=== AMG COMM STRONG DONE WITH WARNINGS. Results: $RESDIR ==="
fi
find $RESDIR -name "*.csv" -printf "  %p  (%s bytes)\n"
