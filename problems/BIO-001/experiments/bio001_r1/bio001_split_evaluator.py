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

Part B (SYNTHETIC, seeded): preserve the corrected E2 negative result.
For one held-gene test cohort, excluding all functional siblings empties CLEAN
training. E2 is INCONCLUSIVE; it does not establish leakage. The donor comparison
is a POST_HOC diagnostic on held groups 6/7 and donors 0/1/2. LEAKY learns donor
means from non-held groups; CLEAN sees only donor 3 and predicts its training
mean. This diagnostic changes both exposure and estimator and is not a causal
isolation of donor leakage. No real expression matrix is evaluated.

Outputs omit the wall clock and use canonical LF bytes for exact replay.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
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
    """Score a historical synthetic regime.

    With donor_holdout=True, test groups are 6/7, restricted to donors 0/1/2.
    LEAKY training uses groups 0..5 with all donors and predicts donor means
    (test groups are unseen). CLEAN keeps only donor 3 training rows and uses
    a constant global mean. The resulting same-cohort difference is POST_HOC,
    not a preregistered or causally isolated leakage result.
    """
    rng = random.Random(12345)  # split RNG independent of data RNG
    if mode in ("frozen", "donorleak", "donorclean"):
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
        # Codex convergence P1: CLEAN and LEAKY must differ in DONOR EXPOSURE
        # while scoring IDENTICAL test rows T. T = held-group cells of donors
        # {0,1,2} (multi-donor). LEAKY trains on all non-held rows (T donors'
        # rows included, donor means learnable). CLEAN trains ONLY on rows of
        # the non-T donor (donor 3) with the donor-blind estimator, so no
        # T-donor context enters training. Predictions scored on the same T.
        t_donors = sorted({d for _, _, d, _, _ in test})[:3]
        test = [c for c in test if c[2] in t_donors]
        if mode == "donorclean":
            train = [c for c in train if c[2] not in t_donors]
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
    test_ids = sorted((g, grp, d, b) for g, grp, d, b, _ in test)
    return pearson(preds, trues), mse, len(test), test_ids


def main() -> int:
    out = {}
    out["norman_audit"] = audit_norman()
    cells, gene_group, _ = make_synthetic()
    # Codex P1 (convergence-2): the preregistered E2 comparison must use ONE
    # test set. T = cells of the held-out genes (individual-gene holdout,
    # 770 rows across ALL groups). The two regimes differ ONLY in whether the
    # held genes' functional siblings enter training:
    #   leaky : train = all cells except held genes  (siblings present)
    #   clean : train also removes the held genes' whole groups (no sibling)
    rng = random.Random(12345)
    held_genes = set()
    for grp in sorted(set(ggrp_val for ggrp_val in
                          [g * 8 // 60 for g in range(60)])):
        genes = [g for g in range(60) if g * 8 // 60 == grp]
        rng.shuffle(genes)
        held_genes.update(genes[:max(1, int(len(genes) * 0.2))])
    held_groups = {g * 8 // 60 for g in held_genes}

    def run_same_T(siblings: bool):
        # Codex P1 (convergence-3): CLEAN must ACTUALLY remove functional-
        # sibling rows from training — not merely ignore the group means.
        # Siblings = non-held genes in groups containing held genes. Under
        # the preregistered per-group 20% holdout, EVERY group contains held
        # genes, so sibling removal empties the CLEAN training set entirely:
        # no valid sibling-free estimator exists. This is recorded as the
        # honest outcome (E2 INCONCLUSIVE), not papered over.
        sibling_genes = {g for g in range(60)
                         if g * 8 // 60 in held_groups and g not in held_genes}
        train = [c for c in cells if c[0] not in held_genes]
        if not siblings:
            train = [c for c in train if c[0] not in sibling_genes]
        test = [c for c in cells if c[0] in held_genes]
        gmean = defaultdict(list)
        glob = []
        for g, grp, d, b, delta in train:
            gmean[grp].append(delta)
            glob.append(delta)
        if not glob:
            # sibling-free CLEAN has NO training data under the preregistered
            # design: no valid estimator can be constructed
            return None, None, len(test), sorted(
                (g, grp, d, b) for g, grp, d, b, _ in test), 0
        mu = sum(glob) / len(glob)
        gm = {k: sum(v) / len(v) for k, v in gmean.items()}
        preds, trues = [], []
        for g, grp, d, b, delta in test:
            preds.append(gm.get(grp, mu) if siblings else mu)
            trues.append(delta)
        mse = sum((p - t) ** 2 for p, t in zip(preds, trues)) / len(trues)
        ids = sorted((g, grp, d, b) for g, grp, d, b, _ in test)
        return pearson(preds, trues), mse, len(test), ids, len(train)

    leaky_p, leaky_mse, n_leaky, leaky_rows, leaky_train = run_same_T(siblings=True)
    frozen_result = run_same_T(siblings=False)
    frozen_p = frozen_result[0]
    frozen_mse = frozen_result[1] if frozen_p is not None else None
    frozen_train_n = frozen_result[4]
    same_T = frozen_result is not None and frozen_result[3] == leaky_rows
    clean_reference_note = (
        f"CLEAN training set size after sibling removal: {frozen_train_n}. "
        "Under the preregistered per-group 20% holdout every group contains "
        "held genes, so sibling removal EMPTIES the CLEAN training set: no "
        "valid sibling-free estimator exists. Per owner instruction the "
        "preregistered E2 is recorded as FAILED/INCONCLUSIVE; the revised "
        "design (model ablation) is retained as exploratory/unverified; "
        "its earlier admission and PASS labels are withdrawn.")
    dl_p, dl_mse, dl_n, dl_rows = run_split(cells, gene_group, "donorleak",
                                         donor_holdout=True)
    dc_p, dc_mse, dc_n, dc_rows = run_split(cells, gene_group, "donorclean",
                                         donor_holdout=True)
    same_T_donor = dl_rows == dc_rows
    out["synthetic_demo"] = {
        "n_cells": len(cells),
        "n_test_E2": n_leaky,
        "E2_design": "single held-gene test set T; regimes differ only in "
                     "whether functional siblings enter training",
        "same_test_rows_invariant": same_T,
        "clean_reference_note": clean_reference_note,
        "n_test_cohorts_identical": True,  # single T by construction (770 rows both regimes)
        "frozen_pearson": (round(frozen_p, 4) if frozen_p is not None else None),
        "frozen_train_size": frozen_train_n,
        "leaky_pearson": round(leaky_p, 4), "leaky_mse": round(leaky_mse, 4),
        "donorleak_pearson": round(dl_p, 4),
        "donorclean_pearson": round(dc_p, 4),
        "same_test_rows_invariant_donor": same_T_donor,
        "donor_design": {
            "held_groups": [6, 7], "test_donors": [0, 1, 2],
            "n_test": dl_n, "n_test_clean": dc_n,
            "n_train_leaky": sum(c[1] not in (6, 7) for c in cells),
            "n_train_clean": sum(c[1] not in (6, 7) and c[2] == 3 for c in cells),
            "claim_status": "POST_HOC_EXPOSURE_AND_ESTIMATOR_COMPARISON"
        },
        "group_leak_gap": (round(leaky_p - frozen_p, 4) if frozen_p is not None
                           else None),  # INCONCLUSIVE: no valid CLEAN estimator
        "donor_leak_gap_vs_donorclean": round(dl_p - dc_p, 4),
        "donor_exposure_note": "POST_HOC: test groups 6/7 and donors 0/1/2; "
                               "LEAKY predicts donor means from non-held groups/all donors; "
                               "CLEAN predicts the constant mean from non-held groups/donor 3. "
                               "Exposure and estimator both differ; no isolated leakage claim",
    }
    gap_group = (leaky_p - frozen_p) if frozen_p is not None else None
    gap_donor = dl_p - dc_p
    gap_group = (leaky_p - frozen_p) if frozen_p is not None else None
    gap_donor = dl_p - dc_p
    # Codex P1: the admitted acceptance froze ONLY the group gap >= 0.2.
    # Everything else is a post-hoc diagnostic and must not be called
    # preregistered (history is not rewritten).
    out["pre_registered_checks"] = {
        "E2_group_leak_gap_ge_0.2": (gap_group >= 0.2 if gap_group is not None
                                     else False),  # INCONCLUSIVE -> not satisfied
    }
    out["pre_registered_checks"]["E2_status"] = ("DEMONSTRATED" if gap_group is not None
                                                 else "FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR")
    out["post_hoc_diagnostics"] = {
        "D1_frozen_split_low_signal": (frozen_p < 0.5 if frozen_p is not None
                                       else None),  # INCONCLUSIVE
        "D3_donor_leak_gap_ge_0.1": gap_donor >= 0.1,
        "D4_leaky_mse_better_than_frozen": (leaky_mse < frozen_mse
                                            if frozen_mse is not None else None),
    }
    e2_constructible = frozen_p is not None
    out["verdict"] = ("LEAKAGE_DEMO_PASS" if (e2_constructible and
                      out["pre_registered_checks"]["E2_group_leak_gap_ge_0.2"] and same_T)
                      else "E2_INCONCLUSIVE_SIBLING_FREE_CLEAN_NOT_CONSTRUCTIBLE"
                      if not e2_constructible else "LEAKAGE_DEMO_FAIL")
    out["e2_status"] = ("CONSTRUCTED" if e2_constructible else
                        "FAILED_INCONCLUSIVE_NO_SIBLING_FREE_ESTIMATOR")
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bio001_r1_results.json").write_bytes(
        (json.dumps(out, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps({"norman": {k: out["norman_audit"][k] for k in
                                 ("cells_total", "class_counts", "single_perturbations",
                                  "cells_per_single_perturbation")},
                      "demo": out["synthetic_demo"],
                      "pre_registered": out["pre_registered_checks"],
                      "post_hoc": out["post_hoc_diagnostics"],
                      "verdict": out["verdict"]},
                     ensure_ascii=False, indent=2))
    # E2 INCONCLUSIVE is a legitimate recorded outcome (negative result):
    # exit 0 so the evidence is preserved and pushed, not hidden
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
