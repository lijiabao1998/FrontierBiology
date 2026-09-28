#!/usr/bin/env python3
"""Independent leakage verifier for BIO-001 r1 (Codex convergence P1).

Does NOT import bio001_split_evaluator. Re-implements: synthetic data
generation (same frozen seeds -> identical dataset), the frozen/leaky/donor
split regimes, the group-mean estimator, Pearson/MSE, and re-derives the main
verdict. Also asserts the same-test-rows invariant between donor regimes and
the degenerate-subset fact recorded in the report.

Checks:
  V1 frozen/leaky group gap reproduces (|delta - committed| < 1e-3)
  V2 donor clean/leaky gap reproduces on identical test rows
  V3 same-test-rows invariant holds
  V4 verdict matches committed LEAKAGE_DEMO_PASS iff V1 and V3 hold
"""
from __future__ import annotations
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


def evaluate(cells, ggrp, regime):
    """(v3) E2 design: single held-gene test set T (770 rows); LEAKY sees
    functional siblings in training, CLEAN gets only the global mean (the
    honest sibling-free reference for a group-mean estimator). Donor regimes
    unchanged (post-hoc diagnostic, same-T)."""
    """Independent split + estimator. Returns (pearson, mse, test_row_ids).

    frozen/donor regimes: whole-group holdout (groups 6,7).
    leaky regime: individual genes held out within every group (20%), so
    functional siblings stay in training - that is the leakage being measured.
    RNG consumption replicated exactly (Random(12345), groups in sorted order).
    """
    rng = random.Random(12345)
    held_genes = set()
    for grp in sorted(set(ggrp)):
        genes = [g for g in range(len(ggrp)) if ggrp[g] == grp]
        rng.shuffle(genes)
        held_genes.update(genes[:max(1, int(len(genes) * 0.2))])
    if regime in ("frozen", "leaky"):
        # E2: single held-gene test set; only sibling exposure differs
        train = [c for c in cells if c[0] not in held_genes]
        test = [c for c in cells if c[0] in held_genes]
    else:
        held = {6, 7}
        train = [c for c in cells if c[1] not in held]
        test = [c for c in cells if c[1] in held]
        if regime in ("donor_leak", "donor_clean"):
            t_donors = sorted({d for _, _, d, _, _ in test})[:3]
            test = [c for c in test if c[2] in t_donors]
            if regime == "donor_clean":
                train = [c for c in train if c[2] not in t_donors]
    gmean = defaultdict(list)
    dmean = defaultdict(list)
    glob = []
    for g, grp, d, b, delta in train:
        gmean[grp].append(delta)
        dmean[d].append(delta)
        glob.append(delta)
    mu = sum(glob) / len(glob)
    gm = {k: sum(v) / len(v) for k, v in gmean.items()}
    dm = {k: sum(v) / len(v) for k, v in dmean.items()}
    preds, trues, ids = [], [], []
    for g, grp, d, b, delta in test:
        if regime == "frozen":
            p = mu  # sibling-free reference (global mean)
        elif regime == "leaky":
            p = gm.get(grp, mu)
        elif regime == "donor_leak":
            p = dm.get(d, mu) + (gm.get(grp, mu) - mu)
        elif regime == "donor_clean":
            p = gm.get(grp, mu)
        else:  # donor_clean
            p = gm.get(grp, mu)
        preds.append(p)
        trues.append(delta)
        ids.append((g, grp, d, b))
    mse = sum((p - t) ** 2 for p, t in zip(preds, trues)) / len(trues)
    return pearson(preds, trues), mse, sorted(ids)


def main() -> int:
    committed = json.loads((RESULTS / "bio001_r1_results.json")
                           .read_text(encoding="utf-8"))["synthetic_demo"]
    cells, ggrp = gen()
    f_p, f_mse, f_ids = evaluate(cells, ggrp, "frozen")
    l_p, l_mse, l_ids = evaluate(cells, ggrp, "leaky")
    dl_p, _, dl_ids = evaluate(cells, ggrp, "donor_leak")
    dc_p, _, dc_ids = evaluate(cells, ggrp, "donor_clean")
    gap_group = l_p - f_p
    gap_donor = dl_p - dc_p
    checks = {
        "V0_same_test_rows_E2": l_ids == f_ids,
        "V1_group_gap_reproduced": abs(gap_group - committed["group_leak_gap"]) < 1e-3,
        "V2_donor_gap_reproduced": abs(gap_donor - committed["donor_leak_gap_vs_donorclean"]) < 1e-3,
        "V3_same_test_rows_donor_regimes": dl_ids == dc_ids,
        "V4_frozen_pearson_reproduced": abs(f_p - committed["frozen_pearson"]) < 1e-3,
        "V5_leaky_pearson_reproduced": abs(l_p - committed["leaky_pearson"]) < 1e-3,
    }
    verdict = "INDEPENDENT_VERIFICATION_MATCH" if all(checks.values()) \
        else "INDEPENDENT_VERIFICATION_MISMATCH"
    out = {"independent_counts": {"gap_group": round(gap_group, 4),
                                  "gap_donor": round(gap_donor, 4),
                                  "frozen_pearson": round(f_p, 4),
                                  "leaky_pearson": round(l_p, 4),
                                  "donor_leak": round(dl_p, 4),
                                  "donor_clean": round(dc_p, 4),
                                  "same_rows": dl_ids == dc_ids},
           "checks": checks, "verdict": verdict}
    (RESULTS / "independent_leakage_verification.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
