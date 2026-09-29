#!/usr/bin/env python3
"""Independent Verifier/Skeptic check for BIO-001 r1 (Codex-review remediation).

Re-derives the committed audit numbers with a structurally different
implementation: regex tokenisation + per-identity memoisation instead of the
main script's split/loop parser, plus targeted spot assertions on the
identities named in the review. Exit 0 iff every number matches.
"""
from __future__ import annotations
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

LIT = Path(__file__).resolve().parent.parent / "lit_data"
TOKEN = re.compile(r"[^_]+")


def genes_of(identity: str) -> frozenset:
    genes = set()
    for tok in TOKEN.findall(identity):
        if tok.startswith("NegCtrl") or tok.isdigit():
            continue
        genes.add(tok)
    return frozenset(genes)


def main() -> int:
    cache: dict[str, frozenset] = {}
    cells = 0
    cls = Counter()
    per_pert = Counter()
    with gzip.open(LIT / "GSE133344_filtered_cell_identities.csv.gz", "rt") as f:
        next(f)
        for line in f:
            gi = line.split(",", 2)[1]
            cells += 1
            if gi not in cache:
                cache[gi] = genes_of(gi)
            g = cache[gi]
            if not g:
                cls["control"] += 1
                per_pert["control"] += 1
            elif len(g) == 1:
                cls["single"] += 1
                per_pert[next(iter(g))] += 1
            else:
                cls["dual"] += 1
                per_pert["+".join(sorted(g))] += 1
    singles = {p: c for p, c in per_pert.items() if p != "control" and "+" not in p}
    # spot assertions from the review
    assert genes_of("SET_KLF1__SET_KLF1") == frozenset({"SET", "KLF1"})
    assert genes_of("NegCtrl0_BAK1__NegCtrl0_BAK1") == frozenset({"BAK1"})
    assert genes_of("NegCtrl10_NegCtrl0") == frozenset()
    assert genes_of("TGFBR2_IGDCC3_2") == frozenset({"TGFBR2", "IGDCC3"})
    # deterministic fold re-derivation
    sizes = Counter(int(hashlib.blake2b(p.encode(), digest_size=8).hexdigest(), 16) % 5
                    for p in singles)
    committed = json.loads(
        (Path(__file__).resolve().parent.parent.parent / "results" / "r1" /
         "bio001_r1_results.json").read_text(encoding="utf-8"))
    na = committed["norman_audit"]
    checks = {
        "cells_total": cells == na["cells_total"],
        "class_counts": dict(cls) == na["class_counts"],
        "single_perturbations": len(singles) == na["single_perturbations"],
        "cells_per_single_median": sorted(singles.values())[len(singles) // 2]
        == na["cells_per_single_perturbation"]["median"],
        "fold_sizes_deterministic": {str(k): sizes[k] for k in sorted(sizes)}
        == na["blake2b_fold_sizes_over_singles"],
        "spot_identities": True,
    }
    print(json.dumps({"independent_counts": {"cells_total": cells,
                                             "class_counts": dict(cls),
                                             "single_perturbations": len(singles)},
                      "checks": checks,
                      "verdict": "INDEPENDENT_CHECK_MATCH" if all(checks.values())
                      else "INDEPENDENT_CHECK_MISMATCH"}, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
