#!/usr/bin/env python3
"""
Aggregate the CSVs written by amg_d_comm_test.

Reads one or more result directories (the per-scale-point directories the sbatch
script creates) and prints, per rank count:

  * the uniform ranking, i.e. which scheme wins over the whole hierarchy;
  * the noise floor, measured rather than assumed: levels whose halo width is
    zero are replicated and exchange nothing, so any apparent effect of changing
    their scheme is pure measurement scatter. That number is the threshold below
    which nothing else in the table can be believed;
  * the per-level sensitivity, with effects inside the noise floor marked so
    they are not read as signal.

Usage:
    ./analyse_amg_comm.py results_amg_comm_123456/*/amg_comm.csv
    ./analyse_amg_comm.py runs/amg_comm.csv
"""

import csv
import sys
from collections import defaultdict


def mean(v):
    return sum(v) / len(v)


def stdev(v):
    if len(v) < 2:
        return 0.0
    m = mean(v)
    return (sum((x - m) ** 2 for x in v) / (len(v) - 1)) ** 0.5


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def levels_path(path):
    return path + ".levels.csv"


def analyse(path):
    rows = load(path)
    if not rows:
        print(f"{path}: empty")
        return

    # Halo width per level: identical across runs, take it from the first one.
    halo = {}
    try:
        for r in load(levels_path(path)):
            halo[int(r["level"])] = int(r["halo_width"])
    except FileNotFoundError:
        pass

    nranks = rows[0]["nranks"]
    idim = rows[0]["idim"]
    nlev = rows[0]["nlev_built"]
    print("=" * 78)
    print(f"{path}")
    print(f"  ranks={nranks}  idim={idim}  levels_built={nlev}")
    if halo:
        desc = "  ".join(f"L{k}:{v}" for k, v in sorted(halo.items()))
        print(f"  halo width per level: {desc}")
    print("=" * 78)

    solve = defaultdict(list)
    iters = set()
    for r in rows:
        key = (r["mode"], r["scheme"], int(r["target_level"]))
        solve[key].append(float(r["krylov_s"]))
        iters.add(int(r["iters"]))

    if len(iters) != 1:
        print(f"  !! iteration count is not constant: {sorted(iters)}")
        print("     the configurations are not solving the same problem;")
        print("     the timings below are not comparable.")
    base_key = ("uniform", "isend_irecv", 0)
    if base_key not in solve:
        print("  !! no uniform baseline run found, cannot compute deltas")
        return
    base = mean(solve[base_key])

    print("\nUNIFORM  (same scheme on every level)")
    print(f"  {'scheme':<28}{'solve [s]':>14}{'sd':>12}{'vs base':>10}")
    uni = [(k, v) for k, v in solve.items() if k[0] == "uniform"]
    for k, v in sorted(uni, key=lambda kv: mean(kv[1])):
        print(f"  {k[1]:<28}{mean(v):>14.4e}{stdev(v):>12.2e}"
              f"{100*(mean(v)-base)/base:>9.1f}%")

    # Noise floor from the levels that cannot possibly be affected.
    dead = [lv for lv, h in halo.items() if h == 0]
    floor = None
    if dead:
        dev = [abs(100 * (mean(v) - base) / base)
               for k, v in solve.items()
               if k[0] == "sensitivity" and k[2] in dead]
        if dev:
            floor = max(dev)
            print(f"\nNOISE FLOOR  {floor:.1f}%")
            print(f"  measured on level(s) {dead}, whose halo width is zero:")
            print("  they are replicated and exchange nothing, so their apparent")
            print("  effect is scatter. Nothing smaller than this is a result.")

    sens = [(k, v) for k, v in solve.items() if k[0] == "sensitivity"]
    if sens:
        print("\nSENSITIVITY  (baseline everywhere except one level)")
        print(f"  {'level':>6}  {'scheme':<28}{'solve [s]':>14}{'vs base':>10}   verdict")
        for k, v in sorted(sens, key=lambda kv: (kv[0][2], mean(kv[1]))):
            delta = 100 * (mean(v) - base) / base
            if floor is not None and abs(delta) <= floor:
                verdict = "noise"
            elif k[2] in dead:
                verdict = "(no exchange at this level)"
            elif delta < 0:
                verdict = "GAIN"
            else:
                verdict = "cost"
            print(f"  {k[2]:>6}  {k[1]:<28}{mean(v):>14.4e}{delta:>9.1f}%   {verdict}")

        if floor is not None:
            gains = [(k[2], k[1], 100 * (mean(v) - base) / base)
                     for k, v in sens
                     if k[2] not in dead and 100 * (mean(v) - base) / base < -floor]
            print()
            if gains:
                print("  Levels where a non-baseline scheme beats the baseline by more")
                print("  than the noise floor -- the case for a per-level policy:")
                for lv, sc, d in sorted(gains, key=lambda g: g[2]):
                    print(f"    level {lv}  {sc}  {d:.1f}%")
            else:
                print("  No per-level gain above the noise floor in this run.")
    print()


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    for path in argv[1:]:
        if path.endswith(".levels.csv"):
            continue
        analyse(path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
