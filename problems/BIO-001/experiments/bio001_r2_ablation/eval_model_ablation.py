#!/usr/bin/env python3
"""BIO-001 exploratory model-ablation replay (not an admitted research round).

The earlier admission/completion/PASS labels are withdrawn: no numerical pass
threshold was frozen, the recorded searches were reused, and no independent r2
verifier exists. Preserve descriptive synthetic values only. This script has no
scientific PASS/FAIL cutoff. Exit 0 means its descriptive output was written,
not that the result is independently verified or the research protocol passed.
The group-aware and structure-blind estimators use the same existing train/test
lists (seeds 42/12345); this is not a split/data-leakage demonstration.
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


def main() -> int:
    cells, ggrp = gen()
    rng = random.Random(12345)
    held_genes = set()
    for grp in sorted(set(ggrp)):
        genes = [g for g in range(len(ggrp)) if ggrp[g] == grp]
        rng.shuffle(genes)
        held_genes.update(genes[:max(1, int(len(genes) * 0.2))])
    train = [c for c in cells if c[0] not in held_genes]
    test = [c for c in cells if c[0] in held_genes]

    gmean = defaultdict(list)
    glob = []
    for g, grp, d, b, delta in train:
        gmean[grp].append(delta)
        glob.append(delta)
    mu = sum(glob) / len(glob)
    gm = {k: sum(v) / len(v) for k, v in gmean.items()}

    def score(pred_fn):
        preds, trues = [], []
        for g, grp, d, b, delta in test:
            preds.append(pred_fn(g, grp))
            trues.append(delta)
        mse = sum((p - t) ** 2 for p, t in zip(preds, trues)) / len(trues)
        return pearson(preds, trues), mse

    ga_p, ga_mse = score(lambda g, grp: gm.get(grp, mu))
    sb_p, sb_mse = score(lambda g, grp: mu)
    gap = ga_p - sb_p

    # Both scores traverse the same test list; this is an author diagnostic,
    # not an independent validation or a retrospective pass criterion.
    diagnostics = {"shared_test_rows_by_construction": True,
                   "gap_recorded": isinstance(gap, float)}
    verdict = "EXPLORATORY_UNVERIFIED"
    out = {"round_id": "20260928T172735194406Z-glm-BIO-001",
           "design": "MODEL_ABLATION (not a leakage demonstration)",
           "n_test": len(test), "n_train": len(train),
           "group_aware_pearson": round(ga_p, 4),
           "structure_blind_pearson": round(sb_p, 4),
           "ablation_gap": round(gap, 4),
           "group_aware_mse": round(ga_mse, 4),
           "structure_blind_mse": round(sb_mse, 4),
           "diagnostics": diagnostics, "verdict": verdict,
           "independently_verified": False, "success_threshold": None,
           "claim_scope": "Descriptive synthetic model comparison only; "
                          "NO confirmatory, split/data-leakage or biological claim"}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bio001_r2_ablation_results.json").write_bytes(
        (json.dumps(out, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps({k: out[k] for k in ("group_aware_pearson",
                                          "structure_blind_pearson",
                                          "ablation_gap", "diagnostics", "verdict")},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
