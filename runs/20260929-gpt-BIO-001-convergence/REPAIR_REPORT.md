# BIO-001 evidence reconciliation｜2026-09-29

Scope and repair acceptance were recorded in REPAIR_SCOPE.md before edits.
Source PR #1 / base `9590a99939b7e54f48f6aa3f8e31c540e6fe1c4f`; governance pin
`07d2b13051b83215182e411e1612f92f1912d8fb`. This correction supersedes the old
active summaries; it does not add a research round or new literature evidence.

## Corrected conclusions

- r1 E2 remains INCONCLUSIVE: CLEAN has zero training rows and the group gap is
  null. The original acceptance is preserved; no threshold or scientific split
  was adjusted to obtain a PASS.
- r2 admission, finished state and PASS are withdrawn. Its current record is
  DRAFT / EXPLORATORY_UNVERIFIED_ADMISSION_WITHDRAWN, explicitly describing the
  already executed exploratory work. Original acceptance/search/result records
  remain in historical_record_at_9590a99.json. No backdated search or new cutoff
  was supplied; the replay emits descriptive values with no success threshold.
- The donor discrepancy was a wrong-cohort bug in the independent verifier.
  The old verifier trained on 5,230 non-held-gene rows and tested on 4,509 rows
  spanning all groups; 3,929 rows overlapped. Production held out groups 6/7:
  4,504 LEAKY training rows, 1,100 CLEAN training rows, and 1,105 test rows, with
  zero overlap. The independent reconstruction now uses that documented
  estimand. Since test groups are unseen, the group adjustment cancels and the
  LEAKY predictor is donor_mean[d]; CLEAN is a constant. The reproduced gap is
  0.1616. Both exposure and estimator differ, so this is post-hoc description,
  not an isolated donor-leakage claim. The old 0.905/MISMATCH JSON is preserved.
- README, STATUS, runs index, active round records and r1 report now agree.
  Historical reports are clearly separated. BIO-001 remains OPEN; expression
  matrix analysis, biological mechanism and r2 independent validation are absent.

## Validation actually performed

Environment: Python 3.13.5 on Windows, standard library only; existing seeds
42 and 12345, existing public metadata inputs, no package install/download.
Commands and outcomes are also recorded in
`problems/BIO-001/results/r1/convergence_validation.json`.

1. Existing r1 evaluator: exit 0 while preserving E2 INCONCLUSIVE. Exit 0 means
   successful evidence serialization, not scientific success.
2. Independent metadata implementation: 6/6 checks match; output retained in
   independent_metadata_verification.json.
3. Independent r1 reconstruction: 8/8 checks match. The expected verdict is
   derived from its own CLEAN cohort, not trusted result status fields.
4. r2 descriptive replay: exit 0, EXPLORATORY_UNVERIFIED, no threshold and no
   volatile timestamp.
5. Convergence regression: 7/7 tests pass, including seven deliberate metric,
   cohort and status corruptions; invalid evidence makes the verifier CLI exit 1;
   changed/empty manifest negative cases fail; four producer/verifier replays
   preserve exact LF result bytes. The first test run exposed a test-harness
   dotted-path bug at the literal key suffix `0.2`; the failure is retained in
   convergence_validation.json and was fixed without changing scientific data.
   The follow-up metadata negative control alters an audit total, seeds stale MATCH
   output, then requires exit 1 and replacement by the newly computed MISMATCH.
6. Manifest: 15/15 entries match working bytes. `verify_manifest.py --git-ref :`
   verifies the staged snapshot and `--git-ref HEAD` verifies committed blobs,
   including the manifest itself. Text artifacts changed by this repair use LF;
   original CRLF source normalization accounts for part of the diff.
7. The pinned registry validator passes all ten problem records. This validates
   metadata structure only; r2 is not admitted. Original r1 search fields remain
   historical; this maintenance correction does not claim a new search.

The proposed extra evidence CI job could not be pushed: GitHub rejected workflow
modification because the OAuth App lacks workflow scope. The workflow change was
removed; credentials and repository permissions were not changed. Existing records
CI remains metadata-only. Scientific replay and exact Git-blob checks above were
performed locally; Linux evidence CI is NOT_RUN. CI and parent review must be
assessed on the pushed HEAD. The author has not merged or approved it.

## Remaining limits

No independent r2 scientific verification is claimed. No current-round search or
preregistration was recovered for r2. Any renewed research requires a new honest
preflight and frozen design. Existing data-license/source limitations in r1 are
unchanged; no clinical, causal or real-expression result is asserted. Independent
review of this correction and the merge decision belong to the parent reviewer.

## Automated review follow-up

Review 5353487166 on c00ed96 identified one P2: metadata replay only printed its
new payload and could leave the saved verification artifact stale. The checker
now writes the same payload to independent_metadata_verification.json using LF
bytes, whether comparison succeeds or fails. The successful replay preserves exact
bytes, and an isolated corrupted-audit case replaces a seeded stale MATCH with
MISMATCH and exits 1. The seven-test suite passes. This changes evidence persistence,
not the metadata calculation, scientific inputs, E2 verdict or r2 status.
