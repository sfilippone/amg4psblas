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
#SBATCH --nodes=16                  # reserve the maximum up front
#SBATCH --ntasks-per-node=112       # full GPP node: 2 x Sapphire Rapids 8480+
#SBATCH --cpus-per-task=1           # CPU-only: one core per rank
#SBATCH --time=02:00:00
#SBATCH --qos=gp_ehpc               # production; gp_debug caps at 2h / 3584 proc
#SBATCH --account=ehpc859           # the GPP project; ehpc580 is ACC and over budget
#
# Verify the node width before trusting RANK_POINTS below:
#     sinfo -o "%P %c"
# If a GPP node is not 112 cores, ntasks-per-node, RANK_POINTS and the NNODES
# arithmetic all have to move together.

# ============================================================================
# USER CONFIGURATION
# ============================================================================
# STRONG or WEAK. They answer different questions and both are worth having:
#
#   strong  fixed problem, growing ranks. Matches the SpMV paper, and shows the
#           knee where communication starts to dominate. Caveat: as ranks grow
#           the AMG hierarchy itself changes shape -- fewer rows per rank at
#           every level, coarse levels degenerating -- so a scale point differs
#           from the next by more than just the rank count.
#
#   weak    fixed rows per rank, problem grown with the ranks. The hierarchy
#           keeps its shape, so the rank count is isolated from the coarsening.
#           This is the cleaner design for "how do the schemes behave as the
#           machine grows", which is the question a per-level policy needs.
SCALING=${SCALING:-strong}

NREP=7                              # repetitions per configuration
NLEV=5                              # max levels requested of the hierarchy
ITMAX=500

RANK_POINTS="112 224 448 896 1792"  # 1, 2, 4, 8, 16 full nodes

# strong: one size for every scale point.
# 400^3 = 64M unknowns; at 1792 ranks that is ~36k rows/rank on the fine level
# and still ~550 on level 3, so the hierarchy stays meaningful to the end.
DIM_STRONG=400

# weak: idim grows as ranks^(1/3) to hold rows/rank constant at ~36.5k,
# matching the strong run at its largest point so the two are comparable there.
DIM_WEAK="160 202 254 320 403"

# Two passes are needed and they answer different questions:
#   PROFILE=false  clean timings -- the numbers that compare schemes
#   PROFILE=true   Score-P attribution -- where the time goes, per level
# The binary is built with scorep-mpifort, so profiling is ON unless disabled:
# a "clean" run is only clean if these are set explicitly.
PROFILE=${PROFILE:-false}

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

if [ "$PROFILE" = "true" ]; then
    export SCOREP_ENABLE_PROFILING=true
    export SCOREP_ENABLE_TRACING=false
    export SCOREP_TOTAL_MEMORY=128M
else
    export SCOREP_ENABLE_PROFILING=false
    export SCOREP_ENABLE_TRACING=false
fi

RESDIR=$SLURM_SUBMIT_DIR/results_amg_${SCALING}_${SLURM_JOB_ID}
mkdir -p $RESDIR

# Region filter: exclude the high-frequency tiny helpers so the profile reflects
# real MPI time rather than instrumentation of the descriptor bookkeeping.
FILTER=$RESDIR/scorep.filt
cat > $FILTER <<'FILTEOF'
SCOREP_REGION_NAMES_BEGIN
  EXCLUDE
    psb_indx_map_mod::*
    psb_desc_mod::*
    psb_error_mod::*
    psb_gen_block_map_mod::*
    psi_penv_mod::*
    psb_hash_mod::*
SCOREP_REGION_NAMES_END
FILTEOF

case "$SCALING" in
    strong) echo "=== AMG communication schemes, STRONG scaling (CPU-only) ===" ;;
    weak)   echo "=== AMG communication schemes, WEAK scaling (CPU-only) ===" ;;
    *)      echo "FATAL: SCALING must be 'strong' or 'weak', got '$SCALING'"; exit 1 ;;
esac
echo "  nrep=$NREP max_levels=$NLEV itmax=$ITMAX"
echo "  PROFILE=$PROFILE"
echo "  reserved_nodes=$SLURM_NNODES rank_points=[$RANK_POINTS]"
echo "  exe=$EXE"
echo "============================================================"

if [ ! -x "$EXE" ]; then
    echo "FATAL: $EXE not found or not executable. Build it before submitting:"
    echo "  cd samples/advanced/pdegen && make amg_d_comm_test"
    exit 1
fi

FAILED=0

IDX=0
for NRANKS in $RANK_POINTS; do
    IDX=$((IDX+1))
    NNODES=$(( (NRANKS + 111) / 112 ))

    if [ "$SCALING" = "weak" ]; then
        DIM=$(echo $DIM_WEAK | cut -d' ' -f$IDX)
    else
        DIM=$DIM_STRONG
    fi

    STEP_DIR=$RESDIR/${NRANKS}ranks
    mkdir -p $STEP_DIR

    echo ""
    echo ">>> $SCALING point: $NRANKS ranks ($NNODES nodes), dim=$DIM"

    # Name the experiment after the scale point. Left to itself Score-P writes
    # scorep-<timestamp>, which carries no indication of the rank count and has
    # to be matched back by hand afterwards.
    if [ "$PROFILE" = "true" ]; then
        export SCOREP_EXPERIMENT_DIRECTORY=$STEP_DIR/scorep_${NRANKS}ranks
        export SCOREP_FILTERING_FILE=$FILTER
    fi

    srun -N $NNODES -n $NRANKS --ntasks-per-node=112 --cpus-per-task=1 \
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
    echo "=== AMG COMM ${SCALING} DONE, all checks passed. Results: $RESDIR ==="
else
    echo "=== AMG COMM ${SCALING} DONE WITH WARNINGS. Results: $RESDIR ==="
fi
find $RESDIR -name "*.csv" -printf "  %p  (%s bytes)\n"
