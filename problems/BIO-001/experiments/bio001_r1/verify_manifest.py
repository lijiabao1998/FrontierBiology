#!/usr/bin/env python3
"""Verify the BIO-001 manifest against working files or exact Git blobs.

Use --git-ref HEAD for committed bytes or --git-ref : for staged bytes.
Exit nonzero for missing files, malformed entries or any mismatch.
"""
import argparse
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
REPO = ROOT.parent.parent
MANIFEST = "results/r1/hashes.txt"


def check_entries(lines, read_bytes):
    statuses = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        try:
            digest, name = line.split(maxsplit=1)
            name = name.removeprefix("*")
            valid = len(digest) == 64 and all(ch in "0123456789abcdef" for ch in digest)
            if not valid or Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("invalid manifest entry")
            actual = hashlib.sha256(read_bytes(name)).hexdigest()
            statuses.append(("OK" if actual == digest else "MISMATCH", name))
        except (OSError, ValueError, subprocess.CalledProcessError):
            statuses.append(("MISSING_OR_INVALID", line))
    if not statuses:
        statuses.append(("EMPTY_MANIFEST", MANIFEST))
    return statuses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--git-ref", help="Verify Git blobs at this ref; ':' means index")
    args = parser.parse_args()
    def read_bytes(name):
        if args.git_ref is None:
            return (ROOT / name).read_bytes()
        prefix = "" if args.git_ref == ":" else args.git_ref
        path = (Path("problems/BIO-001") / name).as_posix()
        return subprocess.check_output(["git", "-C", str(REPO), "show", prefix + ":" + path],
                                       stderr=subprocess.PIPE)
    try:
        lines = read_bytes(MANIFEST).decode("utf-8").splitlines()
    except (OSError, UnicodeError, subprocess.CalledProcessError) as exc:
        print("MANIFEST_UNAVAILABLE", type(exc).__name__)
        return 1
    statuses = check_entries(lines, read_bytes)
    for status, name in statuses:
        print(status, name)
    return 0 if all(status == "OK" for status, _ in statuses) else 1


if __name__ == "__main__":
    raise SystemExit(main())
