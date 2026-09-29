#!/usr/bin/env python3
"""Independent leakage verifier for BIO-001 r1 (v4, convergence-3).

Does NOT import bio001_split_evaluator. Independently re-derives the r1
outcome under its FINAL (post-remediation) semantics:

  - sibling-free CLEAN training set is EMPTY (E2 INCONCLUSIVE root cause)
  - committed group_leak_gap is null and E2_status is
    FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR
  - leaky Pearson reproduces; donor post-hoc diagnostic reproduces on
    identical test rows
  - the COMMITTED top-level verdict matches the verdict independently
    recomputed from the production verdict logic (Codex P1)

Writes canonical LF bytes (Codex P2). Exit 0 iff all checks pass.
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


def main() -> int:
    committed = json.loads((RESULTS / "bio001_r1_results.json")
                           .read_text(encoding="utf-8"))
    demo = committed["synthetic_demo"]
    cells, ggrp = gen()

    # held genes, replicated exactly (Random(12345), groups in sorted order)
    rng = random.Random(12345)
    held = set()
    for grp in sorted(set(ggrp)):
        genes = [g for g in range(len(ggrp)) if ggrp[g] == grp]
        rng.shuffle(genes)
        held.update(genes[:max(1, int(len(genes) * 0.2))])
    held_groups = {ggrp[g] for g in held}
    sibling_genes = {g for g in range(len(ggrp))
                     if ggrp[g] in held_groups and g not in held}

    # 1) sibling-free CLEAN training set must be EMPTY (E2 root cause)
    clean_train = [c for c in cells
                   if c[0] not in held and c[0] not in sibling_genes]

    # 2) leaky regime reproduction (siblings present)
    train = [c for c in cells if c[0] not in held]
    gmean = defaultdict(list)
    glob = []
    for g, grp, d, b, delta in train:
        gmean[grp].append(delta)
        glob.append(delta)
    mu = sum(glob) / len(glob)
    gm = {k: sum(v) / len(v) for k, v in gmean.items()}
    preds, trues = [], []
    for g, grp, d, b, delta in cells:
        if g in held:
            preds.append(gm.get(grp, mu))
            trues.append(delta)
    l_p = pearson(preds, trues)

    # 3) donor post-hoc diagnostic on identical rows
    def donor_regime(leaky: bool):
        tr = [c for c in cells if c[0] not in held]
        te = [c for c in cells if c[1] in held_groups]
        t_donors = sorted({d for _, _, d, _, _ in te})[:3]
        te = [c for c in te if c[2] in t_donors]
        if not leaky:
            tr = [c for c in tr if c[2] not in t_donors]
        dm = defaultdict(list)
        g2 = defaultdict(list)
        gl = []
        for g, grp, d, b, delta in tr:
            dm[d].append(delta)
            g2[grp].append(delta)
            gl.append(delta)
        mu2 = sum(gl) / len(gl)
        dmm = {k: sum(v) / len(v) for k, v in dm.items()}
        gmm = {k: sum(v) / len(v) for k, v in g2.items()}
        pr, tu, ids = [], [], []
        for g, grp, d, b, delta in te:
            p = (dmm.get(d, mu2)) + \
                ((gmm.get(grp, mu2) - mu2) if leaky else 0.0)
            pr.append(p)
            tu.append(delta)
            ids.append((g, grp, d, b))
        return pearson(pr, tu), sorted(ids)

    dl_p, dl_ids = donor_regime(leaky=True)
    dc_p, dc_ids = donor_regime(leaky=False)

    # 4) independently recompute the expected top-level verdict from the
    #    production verdict logic and compare with the COMMITTED verdict
    expected = ("LEAKAGE_DEMO_PASS"
                if committed["pre_registered_checks"].get("E2_group_leak_gap_ge_0.2") is True
                and demo.get("same_test_rows_invariant")
                else "E2_INCONCLUSIVE_SIBLING_FREE_CLEAN_NOT_CONSTRUCTIBLE"
                if committed["pre_registered_checks"].get("E2_status")
                == "FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR"
                else "LEAKAGE_DEMO_FAIL")

    checks = {
        "V1_clean_train_empty": len(clean_train) == 0,
        "V2_inconclusive_reproduced": (demo.get("group_leak_gap") is None
                                       and committed.get("e2_status")
                                       == "FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR"),
        # V3 documents a REAL discrepancy: the independent implementation
        # yields donor gap 0.905 vs production 0.1616 (the two donor_leak
        # estimators differ in group-adjustment exposure). Recorded honestly
        # as a documented discrepancy; reconciliation is queued for the next
        # session. NOT forced to match.
        "V3_donor_gap_reproduced": abs(
            (dl_p - dc_p) - demo["donor_leak_gap_vs_donorclean"]) < 1e-3,
        "V3_documented_discrepancy": {
            "independent_gap": round(dl_p - dc_p, 4),
            "production_gap": demo["donor_leak_gap_vs_donorclean"],
            "interpretation": "donor diagnostic is estimator-sensitive; the "
                              "POST_HOC donor gap must not be quoted without "
                              "stating the estimator variant"},
        "V4_same_test_rows_donor_regimes": dl_ids == dc_ids,
        "V5_leaky_pearson_reproduced": abs(l_p - demo["leaky_pearson"]) < 1e-3,
        "V6_committed_top_level_verdict_matches": committed["verdict"] == expected,
    }
    verdict = "INDEPENDENT_VERIFICATION_MATCH" if all(checks.values()) \
        else "INDEPENDENT_VERIFICATION_MISMATCH"
    out = {"independent_counts": {"clean_train_size": len(clean_train),
                                  "leaky_pearson": round(l_p, 4),
                                  "donor_leak": round(dl_p, 4),
                                  "donor_clean": round(dc_p, 4),
                                  "donor_gap": round(dl_p - dc_p, 4),
                                  "same_rows": dl_ids == dc_ids,
                                  "expected_verdict": expected,
                                  "committed_verdict": committed["verdict"]},
           "checks": checks, "verdict": verdict}
    (RESULTS / "independent_leakage_verification.json").write_bytes(
        (json.dumps(out, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps(out, ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
