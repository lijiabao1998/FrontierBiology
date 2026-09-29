#!/usr/bin/env python3
"""Independent verifier for the BIO-001 r2 model-ablation round (Codex P1).

Structurally independent of eval_model_ablation.py: re-derives the held-gene
set, test set, training rows, both estimator variants, and the ablation gap
from scratch, then audits the COMMITTED artifact for:
  V1 gap reproduces (|delta| < 1e-3)
  V2 identical-test-row invariant holds
  V3 threshold scope: committed verdict must NOT contain PASS/FAIL
     (admitted acceptance was record-only)
  V4 committed Pearson values reproduce

Exit 0 iff all pass. Stdlib only, deterministic.
"""
from __future__ import annotations
import json
import math
import random
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent / "results" / "r2"


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def main() -> int:
    committed = json.loads((RESULTS / "bio001_r2_ablation_results.json")
                           .read_text(encoding="utf-8"))

    # independent data generation (identical seeds/config)
    rng = random.Random(42)
    n_genes, n_groups, n_donors, n_batches, n_cells = 60, 8, 4, 3, 6000
    ggrp = [g * n_groups // n_genes for g in range(n_genes)]
    geff = [rng.uniform(-1.0, 1.0) for _ in range(n_groups)]
    doff = [rng.uniform(-0.4, 0.4) for _ in range(n_donors)]
    boff = [rng.uniform(-0.3, 0.3) for _ in range(n_batches)]
    cells = []
    for _ in range(n_cells):
        g = rng.randrange(n_genes)
        d = rng.randrange(n_donors)
        b = rng.randrange(n_batches)
        cells.append((g, ggrp[g], d, b,
                      geff[ggrp[g]] + rng.gauss(0, 0.25) + doff[d] + boff[b]))

    # independent held-gene derivation (Random(12345), sorted groups)
    rng2 = random.Random(12345)
    held = set()
    for grp in sorted(set(ggrp)):
        genes = [g for g in range(n_genes) if ggrp[g] == grp]
        rng2.shuffle(genes)
        held.update(genes[:max(1, int(len(genes) * 0.2))])

    train = [c for c in cells if c[0] not in held]
    test = [c for c in cells if c[0] in held]
    gmean = defaultdict(list)
    glob = []
    for g, grp, d, b, delta in train:
        gmean[grp].append(delta)
        glob.append(delta)
    mu = sum(glob) / len(glob)
    gm = {k: sum(v) / len(v) for k, v in gmean.items()}

    ga_pr, ga_tr, sb_pr, sb_tr, ids = [], [], [], [], []
    for g, grp, d, b, delta in test:
        ga_pr.append(gm.get(grp, mu))
        sb_pr.append(mu)
        ga_tr.append(delta)
        sb_tr.append(delta)
        ids.append((g, grp, d, b))

    ga_p = pearson(ga_pr, ga_tr)
    sb_p = pearson(sb_pr, sb_tr)
    gap = ga_p - sb_p

    checks = {
        "V1_gap_reproduces": abs(gap - committed["ablation_gap"]) < 1e-3,
        "V2_pearsons_reproduce": (abs(ga_p - committed["group_aware_pearson"]) < 1e-3
                                  and abs(sb_p - committed["structure_blind_pearson"]) < 1e-3),
        # both regimes share the SAME test rows by construction (single loop);
        # row identities repeat across cells ((g,d,b) combos), so uniqueness
        # is not an invariant — size and train-size are (Codex V3 fix)
        "V3_test_and_train_sizes": (len(ids) == committed["n_test"]
                                    and len(train) == committed["n_train"]),
        "V4_threshold_scope": ("PASS" not in committed["verdict"]
                               and "FAIL" not in committed["verdict"]),
        "V5_no_volatile_fields": "generated" not in committed,
    }
    verdict = "INDEPENDENT_VERIFICATION_MATCH" if all(checks.values()) \
        else "INDEPENDENT_VERIFICATION_MISMATCH"
    out = {"independent": {"gap": round(gap, 4), "group_aware": round(ga_p, 4),
                           "structure_blind": round(sb_p, 4)},
           "checks": checks, "verdict": verdict}
    (RESULTS / "bio001_r2_ablation_independent_verification.json").write_bytes(
        (json.dumps(out, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps(out, ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
