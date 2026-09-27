#!/usr/bin/env python3
"""BIO-001 round 1: real Norman-2019 metadata audit + split/leakage evaluator.

Part A (REAL DATA): audit GSE133344 filtered_cell_identities.csv.gz (downloaded
from GEO, 1.9 MB). Extracts target gene(s) from guide_identity (convention
"<GENE>_<guide>__<GENE>_<guide>"; guides named NegCtrl<n> are controls),
classifies control / single / dual perturbations, reports per-perturbation cell
counts, and assigns perturbations to 5 deterministic folds via BLAKE2b hashing
(deterministic across machines/runs). The expression matrix (1.1 GB mtx.gz) is
NOT downloaded this round: DATA_TRACTABILITY_BLOCKED (stdlib-only constraint,
see dataset_map.md).

Part B (SYNTHETIC GROUND TRUTH, seeded): a cell table with known functional
groups, donor and batch offsets, and per-gene response deltas. The SAME split
validator + evaluator that will later consume real expression data is run
under three regimes:

  frozen  : whole functional groups held out (function-grouped holdout;
            the group-mean estimator has no same-group training cells and must
            fall back to the global mean)
  leaky   : individual genes held out (their functional siblings remain in
            training, so the group-mean estimator "knows" the held-out effect)
  donorleak : groups frozen but donors overlap between train and test, with a
            per-donor offset estimator memorizing donor offsets

Verdict (pre-registered in round.json): demo succeeds if
  leaky_pearson - frozen_pearson >= 0.2  and  donorleak_pearson - donorclean_pearson >= 0.1.
A high leaky score with a low frozen score QUANTIFIES the leakage inflation that
PerturbVAE and Systema report for published benchmarks. Stdlib only.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIT = HERE.parent / "lit_data"
RESULTS = HERE.parent.parent / "results" / "r1"


def pearson(xs, ys) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def fold_of(key: str, n_folds: int = 5) -> int:
    h = hashlib.blake2b(key.encode("utf-8"), digest_size=8).hexdigest()
    return int(h, 16) % n_folds


# ---------------- Part A: real metadata audit ----------------

def target_genes(guide_identity: str) -> list[str]:
    """Target genes of a guide_identity (Codex-review corrected parser).

    Norman convention: each '__'-half lists the cell's two guide constructs as
    '<gene>_<guide>' tokens; gene names may contain underscores (SET, C11orf53
    style) and guides may be NegCtrl<n> or numeric sgRNA ids. So the target set
    is ALL non-control, non-numeric tokens across both halves (duplicates
    collapse). Examples (real):
      SET_KLF1__SET_KLF1          -> {SET, KLF1}        (dual)
      NegCtrl0_BAK1__NegCtrl0_BAK1-> {BAK1}             (single)
      NegCtrl10_NegCtrl0          -> {}                 (control)
      TGFBR2_IGDCC3_2             -> {TGFBR2, IGDCC3}   (dual, sgRNA id '2')
    """
    genes = []
    for half in guide_identity.split("__"):
        for tok in half.split("_"):
            if tok.startswith("NegCtrl") or tok.isdigit():
                continue
            if tok and tok not in genes:
                genes.append(tok)
    return sorted(genes)


def audit_norman() -> dict:
    path = LIT / "GSE133344_filtered_cell_identities.csv.gz"
    cells = 0
    good = 0
    per_pert = Counter()
    cls = Counter()
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.DictReader(f):
            cells += 1
            if row["good_coverage"] == "True":
                good += 1
            gi = row["guide_identity"]
            tg = target_genes(gi)
            if not tg:
                cls["control"] += 1
                pert = "control"
            elif len(tg) == 1:
                cls["single"] += 1
                pert = tg[0]
            else:
                cls["dual"] += 1
                pert = "+".join(tg)
            per_pert[pert] += 1
    singles = {p: c for p, c in per_pert.items() if "+" not in p and p != "control"}
    folds = Counter(fold_of(p) for p in singles)
    counts = sorted(singles.values())
    n = len(counts)
    return {
        "file": path.name, "cells_total": cells, "cells_good_coverage": good,
        "class_counts": dict(cls),
        "single_perturbations": len(singles),
        "single_perturbation_cells": sum(singles.values()),
        "dual_perturbations": len(per_pert) - len(singles) - (1 if "control" in per_pert else 0),
        "cells_per_single_perturbation": {
            "min": counts[0], "median": counts[n // 2], "max": counts[-1],
            "p10": counts[n // 10], "p90": counts[9 * n // 10]},
        "smallest_5": dict(sorted(singles.items(), key=lambda kv: kv[1])[:5]),
        "largest_5": dict(sorted(singles.items(), key=lambda kv: -kv[1])[:5]),
        "blake2b_fold_sizes_over_singles": {str(k): folds[k] for k in sorted(folds)},
    }


# ---------------- Part B: synthetic ground truth + leakage demo ----------------

def make_synthetic(seed: int = 42, n_genes: int = 60, n_groups: int = 8,
                   n_donors: int = 4, n_batches: int = 3, n_cells: int = 6000):
    rng = random.Random(seed)
    gene_group = [g * n_groups // n_genes for g in range(n_genes)]
    group_effect = [rng.uniform(-1.0, 1.0) for _ in range(n_groups)]
    donor_offset = [rng.uniform(-0.4, 0.4) for _ in range(n_donors)]
    batch_offset = [rng.uniform(-0.3, 0.3) for _ in range(n_batches)]
    cells = []
    for _ in range(n_cells):
        g = rng.randrange(n_genes)
        d = rng.randrange(n_donors)
        b = rng.randrange(n_batches)
        delta = (group_effect[gene_group[g]] + rng.gauss(0, 0.25)
                 + donor_offset[d] + batch_offset[b])
        cells.append((g, gene_group[g], d, b, delta))
    return cells, gene_group, group_effect


def run_split(cells, gene_group, mode: str, held_groups=(6, 7), hold_frac: float = 0.2,
              donor_holdout: bool = False):
    """Returns (pearson, mse, n_test) of the group-mean / global-mean baseline.

    donor_holdout=True (with mode='frozen' or 'donorleak') reserves donors 3+
    for test only: the estimator never sees their offsets, so any donor-mean
    usage on test cells is genuine leakage rather than a mere feature benefit
    (Codex-review corrected comparison).
    """
    rng = random.Random(12345)  # split RNG independent of data RNG
    if mode in ("frozen", "donorleak"):
        test_groups = set(held_groups)
        train = [c for c in cells if c[1] not in test_groups]
        test = [c for c in cells if c[1] in test_groups]
    elif mode == "leaky":
        all_groups = sorted({c[1] for c in cells})
        held_genes = set()
        for grp in all_groups:
            genes = [g for g in range(len(gene_group)) if gene_group[g] == grp]
            rng.shuffle(genes)
            held_genes.update(genes[:max(1, int(len(genes) * hold_frac))])
        train = [c for c in cells if c[0] not in held_genes]
        test = [c for c in cells if c[0] in held_genes]
    if donor_holdout:
        # reserve ONE donor for test only (with 4 donors, holding out all
        # test-observed donors would empty the training pool)
        held_donor = max({d for _, _, d, _, _ in cells})
        train = [c for c in train if c[2] != held_donor]
        test = [c for c in test if c[2] == held_donor]
    # estimators
    grp_mean = defaultdict(list)
    glob, don = [], defaultdict(list)
    for g, grp, d, b, delta in train:
        grp_mean[grp].append(delta)
        glob.append(delta)
        don[d].append(delta)
    global_mean = sum(glob) / len(glob)
    grp_mean = {k: sum(v) / len(v) for k, v in grp_mean.items()}
    donor_mean = {k: sum(v) / len(v) for k, v in don.items()}
    use_donor = mode == "donorleak"
    preds, trues = [], []
    for g, grp, d, b, delta in test:
        if use_donor and d in donor_mean:
            pred = donor_mean[d] + (grp_mean.get(grp, global_mean) - global_mean)
        else:
            pred = grp_mean.get(grp, global_mean)
        preds.append(pred)
        trues.append(delta)
    mse = sum((p - t) ** 2 for p, t in zip(preds, trues)) / len(trues)
    return pearson(preds, trues), mse, len(test)


def main() -> int:
    out = {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    out["norman_audit"] = audit_norman()
    cells, gene_group, _ = make_synthetic()
    frozen_p, frozen_mse, n_frozen = run_split(cells, gene_group, "frozen")
    leaky_p, leaky_mse, n_leaky = run_split(cells, gene_group, "leaky")
    dl_p, dl_mse, _ = run_split(cells, gene_group, "donorleak")
    dc_p, dc_mse, _ = run_split(cells, gene_group, "frozen", donor_holdout=True)
    out["synthetic_demo"] = {
        "n_cells": len(cells), "n_test_frozen": n_frozen, "n_test_leaky": n_leaky,
        "frozen_pearson": round(frozen_p, 4), "frozen_mse": round(frozen_mse, 4),
        "leaky_pearson": round(leaky_p, 4), "leaky_mse": round(leaky_mse, 4),
        "donorleak_pearson": round(dl_p, 4),
        "donorclean_pearson": round(dc_p, 4),
        "group_leak_gap": round(leaky_p - frozen_p, 4),
        "donor_leak_gap_vs_donorclean": round(dl_p - dc_p, 4),
        "donor_leak_gap_note": "Codex-review corrected: donor leakage is now "
                               "measured against a donor-held-out split (test "
                               "donors unseen in training), not merely against "
                               "the estimator feature being switched off.",
    }
    gap_group = leaky_p - frozen_p
    gap_donor = dl_p - dc_p
    out["pre_registered_checks"] = {
        "E1_frozen_split_low_signal": frozen_p < 0.5,
        "E2_group_leak_gap_ge_0.2": gap_group >= 0.2,
        "E3_donor_leak_gap_ge_0.1": gap_donor >= 0.1,
        "E4_leaky_mse_better_than_frozen": leaky_mse < frozen_mse,
    }
    out["verdict"] = ("LEAKAGE_DEMO_PASS" if all(out["pre_registered_checks"].values())
                      else "LEAKAGE_DEMO_FAIL")
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bio001_r1_results.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"norman": {k: out["norman_audit"][k] for k in
                                 ("cells_total", "class_counts", "single_perturbations",
                                  "cells_per_single_perturbation")},
                      "demo": out["synthetic_demo"],
                      "checks": out["pre_registered_checks"],
                      "verdict": out["verdict"]},
                     ensure_ascii=False, indent=2))
    return 0 if all(out["pre_registered_checks"].values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
