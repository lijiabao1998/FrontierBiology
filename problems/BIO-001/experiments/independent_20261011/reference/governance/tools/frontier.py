#!/usr/bin/env python3
"""Dependency-free registry and per-round gate. It checks records, not scientific truth."""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import urlparse

STATUSES = {'OPEN', 'PARTIAL', 'CLAIMED_RESOLVED', 'COMPLETED_EXTERNAL', 'COMPLETED_INTERNAL', 'PAUSED', 'RETRACTED'}
TERMINAL = {'COMPLETED_EXTERNAL', 'COMPLETED_INTERNAL'}
VERDICTS = {'NO_RESOLUTION_FOUND', 'PARTIAL_PROGRESS', 'CLAIMED_RESOLVED', 'RESOLVED_EXTERNAL', 'BLOCKED'}
CATEGORIES = {'general', 'discipline', 'solution', 'criticism'}
UTC = dt.timezone.utc


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load(path: Path) -> dict:
    with path.open(encoding='utf-8') as f:
        obj = json.load(f)
    require(isinstance(obj, dict), f'{path}: expected JSON object')
    return obj


def save(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def stamp(value: str) -> dt.datetime:
    t = dt.datetime.fromisoformat(value.replace('Z', '+00:00'))
    require(t.tzinfo is not None, 'timestamp must include UTC offset')
    return t.astimezone(UTC)


def nonempty(obj: dict, keys: tuple[str, ...], context: str) -> None:
    for k in keys:
        require(isinstance(obj.get(k), str) and bool(obj[k].strip()), f'{context}: missing {k}')


def source(ref: dict) -> None:
    nonempty(ref, ('title', 'url', 'kind', 'supports'), 'source')
    p = urlparse(ref['url'])
    require(p.scheme == 'https' and bool(p.netloc), 'source needs HTTPS URL')
    require(p.netloc not in {'example.com', 'example.org'}, 'placeholder source')


def validate_resolution(obj: dict, problem_id: str) -> None:
    nonempty(obj, ('problem_id', 'scope_statement', 'basis', 'confirmed_by', 'independent_check'), 'resolution')
    require(obj['problem_id'] == problem_id, 'resolution belongs to another problem')
    require(obj.get('scope_match') is True, 'a subproblem does not close its parent')
    require(obj.get('independently_verified') is True, 'resolution needs independent verification')
    refs = obj.get('sources', [])
    require(len({r.get('url') for r in refs}) >= 2, 'resolution needs primary source plus corroboration/certificate')
    for r in refs:
        source(r)
    require(any(r['kind'] == 'primary' for r in refs), 'resolution missing primary source')
    stamp(obj['confirmed_at'])


def validate_preflight(r: dict, fresh: bool = False) -> str:
    nonempty(r, ('round_id', 'problem_id', 'agent', 'branch', 'base_sha', 'acceptance', 'not_done'), 'round')
    require(re.fullmatch(r'[a-z0-9-]+/.+', r['branch']) is not None, 'work on an agent branch, not main')
    require(re.fullmatch(r'[0-9a-f]{40}', r['base_sha']) is not None, 'record exact base commit SHA')
    started = stamp(r['started_at'])
    p = r.get('preflight', {})
    checked = stamp(p['checked_at'])
    require(started <= checked <= started + dt.timedelta(hours=24), 'each round needs its own start-time search')
    if fresh:
        now = dt.datetime.now(UTC)
        require(now - dt.timedelta(hours=24) <= checked <= now + dt.timedelta(minutes=5), 'stale or future search; repeat preflight')
    for name in ('sources_read', 'scope_compared', 'retractions_checked', 'licenses_checked', 'safety_checked', 'baseline_identified', 'acceptance_frozen'):
        require(p.get(name) is True, f'preflight incomplete: {name}')
    nonempty(p, ('limitations', 'scope_comparison', 'baseline'), 'preflight')
    queries = p.get('queries', [])
    require(CATEGORIES <= {q.get('category') for q in queries}, 'need general / discipline / solution / criticism searches')
    require(len({q.get('query') for q in queries}) >= 4, 'need four distinct actual queries')
    for q in queries:
        nonempty(q, ('query', 'engine', 'category', 'outcome'), 'query')
    refs = p.get('sources', [])
    require(len({s.get('url') for s in refs}) >= 2, 'read at least two distinct sources')
    for ref in refs:
        source(ref)
    require(any(s['kind'] == 'primary' for s in refs), 'at least one source must be primary')
    verdict = p.get('verdict')
    require(verdict in VERDICTS, 'unknown novelty verdict')
    budget = r.get('budget', {})
    require(0 < budget.get('wall_minutes', 0) <= 120, 'bounded run: 1–120 minutes')
    require(budget.get('usd', -1) >= 0, 'explicit USD budget required')
    require(0 < budget.get('max_trials', 0), 'explicit trial limit required')
    if verdict == 'RESOLVED_EXTERNAL':
        validate_resolution(p.get('resolution', {}), r['problem_id'])
    return verdict


def problem_path(root: Path, problem_id: str) -> Path:
    require(re.fullmatch(r'(MATH|PHYS|BIO|CHEM)-\d{3}', problem_id) is not None, 'invalid problem id')
    return root / 'problems' / problem_id / 'problem.json'


def validate_problem(p: dict) -> None:
    nonempty(p, ('id', 'title', 'domain', 'statement', 'known_result', 'open_gap', 'first_task', 'evaluator', 'validation_limits', 'completion_criterion', 'screening_note'), 'problem')
    require(p['status'] in STATUSES, 'invalid problem status')
    require(p['priority'] in {'A', 'B', 'C'}, 'priority must be A/B/C')
    date = dt.date.fromisoformat(p['checked_on'])
    require(date <= dt.datetime.now(UTC).date(), 'future literature check')
    require(len(p.get('screening_queries', [])) >= 1, 'record bootstrap search query')
    require(len(p.get('sources', [])) >= 1, 'problem needs a source')
    for ref in p['sources']:
        source(ref)
    if p['status'] in TERMINAL:
        validate_resolution(p.get('resolution', {}), p['id'])


def validate_repo(root: Path) -> None:
    lab = load(root / 'lab.json')
    for filename in ('README.md', 'AGENTS.md', 'STATUS.md'):
        require((root / filename).is_file(), f'missing {filename}')
    if lab['domain'] == 'governance':
        print('Governance documents present; run unit tests separately.')
        return
    files = sorted((root / 'problems').glob('*/problem.json'))
    require(len(files) == lab['expected_problem_count'], 'problem count does not match lab.json')
    seen = set()
    for path in files:
        p = load(path)
        validate_problem(p)
        require(p['id'] not in seen and path.parent.name == p['id'], 'duplicate or mismatched problem id')
        require(p['domain'] == lab['domain'], 'cross-domain problem contamination')
        seen.add(p['id'])
    for path in sorted((root / 'runs').glob('*/round.json')):
        r = load(path)
        require(r['problem_id'] in seen, 'round points to absent problem')
        if r.get('state') != 'DRAFT':
            validate_preflight(r)
    print(f'PASS: {len(files)} problem records; metadata checks only, not proof of openness or discovery.')


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()


def check_diff(root: Path, base: str, head: str) -> None:
    changed = git(root, 'diff', '--name-only', f'{base}...{head}').splitlines()
    affected = set()
    for path in changed:
        m = re.match(r'problems/((?:MATH|PHYS|BIO|CHEM)-\d{3})/(?:experiments|proofs|results)/', path)
        if m:
            affected.add(m.group(1))
    rounds = [load(root / p) for p in changed if re.fullmatch(r'runs/[^/]+/round.json', p) and (root / p).is_file()]
    for pid in affected:
        eligible = []
        for r in rounds:
            if r.get('problem_id') != pid or r.get('state') == 'DRAFT':
                continue
            verdict = validate_preflight(r)
            if verdict in {'NO_RESOLUTION_FOUND', 'PARTIAL_PROGRESS'}:
                eligible.append(r)
        require(bool(eligible), f'{pid}: research changes lack a completed per-round search record')
    print(f'PASS: preflight coverage for {len(affected)} changed research scopes.')


def start(root: Path, pid: str, agent: str) -> None:
    p = load(problem_path(root, pid))
    require(p['status'] not in TERMINAL, 'already completed: select another problem or open a replication task')
    require(re.fullmatch(r'[a-z0-9-]+', agent) is not None, 'use an agent slug')
    now = dt.datetime.now(UTC)
    rid = now.strftime('%Y%m%dT%H%M%S%fZ') + '-' + agent + '-' + pid
    record = {'round_id': rid, 'problem_id': pid, 'agent': agent, 'branch': agent + '/' + pid.lower(), 'base_sha': git(root, 'rev-parse', 'HEAD'), 'started_at': now.isoformat(), 'state': 'DRAFT', 'acceptance': '', 'not_done': '尚未檢索、尚未實驗；本記錄不是研究成果。', 'budget': {'wall_minutes': 30, 'usd': 0, 'max_trials': 100}, 'preflight': {'checked_at': None, 'verdict': 'BLOCKED', 'queries': [], 'sources': [], 'limitations': '', 'scope_comparison': '', 'baseline': ''}, 'artifacts': [], 'result': None}
    path = root / 'runs' / rid / 'round.json'
    save(path, record)
    print(path)
    print('BLOCKED: browse now, read primary sources, complete preflight, then run admit. No autonomous search has occurred.')


def admit(root: Path, path: Path) -> None:
    r = load(path)
    verdict = validate_preflight(r, fresh=True)
    pp = problem_path(root, r['problem_id'])
    p = load(pp)
    require(p['status'] not in TERMINAL, 'target already completed')
    if verdict == 'RESOLVED_EXTERNAL':
        p['status'] = 'COMPLETED_EXTERNAL'
        p['resolution'] = r['preflight']['resolution']
        p['checked_on'] = stamp(r['preflight']['checked_at']).date().isoformat()
        save(pp, p)
        r['state'] = 'CLOSED_EXTERNAL'
        save(path, r)
        print('COMPLETED_EXTERNAL: record updated on this working branch; stop searching for a new solution. PR still requires review.')
    elif verdict in {'CLAIMED_RESOLVED', 'BLOCKED'}:
        r['state'] = 'PAUSED'
        save(path, r)
        raise ValueError('STOP: investigate the resolution claim or restore access; no discovery run admitted')
    else:
        r['state'] = 'ADMITTED'
        r['admitted_at'] = dt.datetime.now(UTC).isoformat()
        save(path, r)
        print('ADMITTED: one budgeted round. Preserve code, seeds, checksums, stdout, failures and limitations.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='cmd', required=True)
    v = sub.add_parser('validate'); v.add_argument('root', type=Path)
    d = sub.add_parser('check-diff'); d.add_argument('root', type=Path); d.add_argument('base'); d.add_argument('head')
    s = sub.add_parser('start'); s.add_argument('root', type=Path); s.add_argument('problem'); s.add_argument('--agent', required=True)
    a = sub.add_parser('admit'); a.add_argument('root', type=Path); a.add_argument('record', type=Path)
    args = parser.parse_args()
    try:
        if args.cmd == 'validate': validate_repo(args.root)
        elif args.cmd == 'check-diff': check_diff(args.root, args.base, args.head)
        elif args.cmd == 'start': start(args.root, args.problem, args.agent)
        elif args.cmd == 'admit': admit(args.root, args.record)
        return 0
    except (ValueError, KeyError, TypeError, OSError, subprocess.CalledProcessError) as e:
        print(f'BLOCKED: {e}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
