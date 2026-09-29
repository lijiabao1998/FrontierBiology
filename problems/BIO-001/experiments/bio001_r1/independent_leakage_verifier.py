#!/usr/bin/env python3
"""Independent reconstruction of corrected BIO-001 r1 evidence.

No import of the production evaluator. E2 and the post-hoc donor comparison
have separate cohorts. Row indices, not nonunique cell attributes, identify
observations. The expected E2 verdict is derived from the reconstructed empty
CLEAN train set, never from the result fields being verified. The historical
v4 mismatch remains in results/r1/history; it used the wrong donor split with
3,929 overlapping train/test rows.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results" / "r1"


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def gen(seed=42, n_genes=60, n_groups=8, n_donors=4, n_batches=3, n_cells=6000):
    rng = random.Random(seed)
    ggrp = [g * n_groups // n_genes for g in range(n_genes)]
    geff = [rng.uniform(-1.0, 1.0) for _ in range(n_groups)]
    doff = [rng.uniform(-0.4, 0.4) for _ in range(n_donors)]
    boff = [rng.uniform(-0.3, 0.3) for _ in range(n_batches)]
    out = []
    for _ in range(n_cells):
        g = rng.randrange(n_genes)
        d = rng.randrange(n_donors)
        b = rng.randrange(n_batches)
        out.append((g, ggrp[g], d, b,
                    geff[ggrp[g]] + rng.gauss(0, 0.25) + doff[d] + boff[b]))
    return out, ggrp


def reconstruct():
    cells, ggrp = gen()
    rng = random.Random(12345)
    held = set()
    for group in sorted(set(ggrp)):
        genes = [gene for gene, grp in enumerate(ggrp) if grp == group]
        rng.shuffle(genes)
        held.update(genes[:max(1, int(len(genes) * 0.2))])
    held_groups = {ggrp[gene] for gene in held}
    e2_test = [i for i, cell in enumerate(cells) if cell[0] in held]
    e2_train = [i for i, cell in enumerate(cells) if cell[0] not in held]
    clean_train = [i for i in e2_train if cells[i][1] not in held_groups]

    # Independently aggregate only E2 training rows by group.
    means = {}
    for group in sorted(set(ggrp)):
        values = [cells[i][4] for i in e2_train if cells[i][1] == group]
        means[group] = sum(values) / len(values)
    predicted = [means[cells[i][1]] for i in e2_test]
    observed = [cells[i][4] for i in e2_test]
    leaky_pearson = pearson(predicted, observed)
    leaky_mse = sum((x-y)**2 for x,y in zip(predicted, observed)) / len(observed)

    # Donor diagnostic: held GROUPS 6/7, not the E2 held-gene split.
    # All test groups are absent from train, so production's
    # donor_mean[d] + (group_mean.get(grp, global_mean) - global_mean)
    # reduces exactly to donor_mean[d]. CLEAN sees only donor 3 and returns
    # a constant. Both exposure and predictor differ; this is post-hoc only.
    donor_test = [i for i,c in enumerate(cells) if c[1] in (6,7) and c[2] in (0,1,2)]
    donor_train = [i for i,c in enumerate(cells) if c[1] not in (6,7)]
    donor_clean = [i for i in donor_train if cells[i][2] == 3]
    donor_means = {}
    for donor in range(4):
        values = [cells[i][4] for i in donor_train if cells[i][2] == donor]
        donor_means[donor] = sum(values) / len(values)
    constant = sum(cells[i][4] for i in donor_clean) / len(donor_clean)
    target = [cells[i][4] for i in donor_test]
    donor_leaky_p = pearson([donor_means[cells[i][2]] for i in donor_test], target)
    donor_clean_p = pearson([constant] * len(donor_test), target)
    empty = len(clean_train) == 0
    return {
        "n_cells": len(cells), "n_test_E2": len(e2_test),
        "clean_train_size": len(clean_train),
        "e2_no_row_overlap": not bool(set(e2_train) & set(e2_test)),
        "leaky_pearson": round(leaky_pearson, 4), "leaky_mse": round(leaky_mse, 3),
        "donor_leak": round(donor_leaky_p, 4),
        "donor_clean": round(donor_clean_p, 4),
        "donor_gap": round(donor_leaky_p-donor_clean_p, 4),
        "donor_test_size": len(donor_test),
        "donor_train_leaky_size": len(donor_train),
        "donor_train_clean_size": len(donor_clean),
        "donor_no_row_overlap": not bool(set(donor_train) & set(donor_test)),
        "expected_verdict": "E2_INCONCLUSIVE_SIBLING_FREE_CLEAN_NOT_CONSTRUCTIBLE"
                            if empty else "UNEXPECTED_CONSTRUCTIBLE_DESIGN",
        "expected_status": "FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR"
                           if empty else "UNEXPECTED_CONSTRUCTIBLE_DESIGN",
    }


def verify(committed):
    actual = reconstruct()
    demo = committed.get("synthetic_demo", {})
    prereg = committed.get("pre_registered_checks", {})
    design = demo.get("donor_design", {})
    checks = {
        "V1_clean_train_empty": actual["clean_train_size"] == 0,
        "V2_empty_estimator_recorded": (demo.get("frozen_train_size") == 0
            and "frozen_pearson" in demo and demo["frozen_pearson"] is None
            and "group_leak_gap" in demo and demo["group_leak_gap"] is None),
        "V3_e2_status_matches": (prereg.get("E2_group_leak_gap_ge_0.2") is False
            and prereg.get("E2_status") == actual["expected_status"]
            and committed.get("e2_status") == actual["expected_status"]),
        "V4_e2_counts_and_rows": (demo.get("n_cells") == actual["n_cells"]
            and demo.get("n_test_E2") == actual["n_test_E2"]
            and demo.get("same_test_rows_invariant") is True
            and demo.get("n_test_cohorts_identical") is True
            and actual["e2_no_row_overlap"]),
        "V5_e2_leaky_metrics": (demo.get("leaky_pearson") == actual["leaky_pearson"]
            and demo.get("leaky_mse") == actual["leaky_mse"]),
        "V6_donor_metrics": (demo.get("donorleak_pearson") == actual["donor_leak"]
            and demo.get("donorclean_pearson") == actual["donor_clean"]
            and demo.get("donor_leak_gap_vs_donorclean") == actual["donor_gap"]),
        "V7_donor_cohort": (design.get("held_groups") == [6,7]
            and design.get("test_donors") == [0,1,2]
            and design.get("n_test") == actual["donor_test_size"]
            and design.get("n_test_clean") == actual["donor_test_size"]
            and design.get("n_train_leaky") == actual["donor_train_leaky_size"]
            and design.get("n_train_clean") == actual["donor_train_clean_size"]
            and demo.get("same_test_rows_invariant_donor") is True
            and actual["donor_no_row_overlap"]),
        "V8_committed_verdict_matches_independent_result":
            committed.get("verdict") == actual["expected_verdict"],
    }
    return {"independent_counts": actual, "checks": checks,
            "verdict": "INDEPENDENT_VERIFICATION_MATCH" if all(checks.values())
                       else "INDEPENDENT_VERIFICATION_MISMATCH",
            "scope": "r1 E2 negative result and fixed post-hoc donor estimand only; not r2",
            "historical_mismatch": "history/independent_leakage_verification.v4.json"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULTS / "bio001_r1_results.json")
    parser.add_argument("--check-only", action="store_true", help="Do not write verification evidence")
    args = parser.parse_args()
    out = verify(json.loads(args.result.read_text(encoding="utf-8")))
    if not args.check_only:
        (RESULTS / "independent_leakage_verification.json").write_bytes(
            (json.dumps(out, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps(out, ensure_ascii=False, indent=2))
    return 0 if all(out["checks"].values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
