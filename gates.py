"""Run the SELECTION.md §2 hard gates locally and write rows into scores.csv.

abguard is authoritative and runs on submission. This pre-checks what can be
checked from sequence alone, so no submission slot is lost to a rejection we
could have predicted ourselves.

Honesty rule, per SELECTION.md §8 ("verified, not assumed"): a gate this script
cannot actually evaluate is reported as **unchecked**, never as a pass. Three of
the nine are unchecked by design -- see UNCHECKED below.

Input (either form):

  CSV   columns: design_id, sequence [, light] [, track] [, job_name] [, contig]
                 a two-chain design may also put both chains in `sequence`
                 separated by ':' or '/'
  FASTA headers: >design_id [track=de_novo] [job=...] [contig=...]
                 chains separated by ':' or '/' in the sequence body

Run:
  cd ~/antibody-hack && source .venv/bin/activate
  python gates.py --self-test                 # sanity-check on trastuzumab
  python gates.py candidates.csv              # print the table
  python gates.py candidates.csv --write      # also append rows to scores.csv
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import sys
from pathlib import Path

try:
    from anarcii import Anarcii
except ImportError:  # import path has differed between versions
    from anarcii.pipeline import Anarcii

import anarcii

HERE = Path(__file__).parent
SCORES_CSV = HERE / "scores.csv"

# Trastuzumab reference, from smoke.py. CDR-H3 = SRWGGDGFYAMDY.
TRAST_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRAST_CDRH3 = "SRWGGDGFYAMDY"

# IMGT boundaries, inclusive (same table as smoke.py).
IMGT_REGIONS = {
    "FR1": (1, 26), "CDR1": (27, 38), "FR2": (39, 55), "CDR2": (56, 65),
    "FR3": (66, 104), "CDR3": (105, 117), "FR4": (118, 128),
}
# The intradomain disulfide: IMGT 23 and 104 are the conserved cysteines.
CONSERVED_CYS = (23, 104)

# Thresholds. G1/G2/G4 are structural facts; G3's is ours, see note.
MAX_TOTAL_AA = 250          # G4
MIN_ANARCII_SCORE = 10.0    # G3 proxy, see the comment on gate_g3()
CDRH3_IDENTITY_MAX = 0.70   # O1
MIN_CDRH3_MUTATIONS = 2     # D2

# Gates this script cannot evaluate from sequence alone.
UNCHECKED = {
    "D1": "needs a known-antibody database (abguard/BLAST) for identity <70%",
    "O2": "needs the provided data points; pass --data-points FILE to check",
    "G1": "format is inferred from chain typing, not asserted by a schema",
}

# scores.csv column order -- must match the committed header exactly.
COLUMNS = [
    "design_id", "job_name", "contig", "generated_at", "gate_pass",
    "gate_failed", "length_aa", "cdr3_seq", "ipTM", "pLDDT", "motif_rmsd",
    "tnp_flags", "epitope_ok", "tier", "portfolio_group", "submitted",
    "score_source_tool", "score_source_version", "notes",
]


# --- small alignment, so identity/mutation counts are not length-naive -----

def _align(a: str, b: str, match=1, mismatch=-1, gap=-2):
    """Needleman-Wunsch. Returns (aligned_a, aligned_b).

    Implemented here rather than pulled from Biopython because the
    PairwiseAligner API has moved between releases and this is 20 lines.
    """
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = gap * i
    for j in range(1, m + 1):
        score[0][j] = gap * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score[i][j] = max(
                score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch),
                score[i - 1][j] + gap,
                score[i][j - 1] + gap,
            )
    out_a, out_b, i, j = [], [], n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + (
                match if a[i - 1] == b[j - 1] else mismatch):
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] + gap:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def identity(a: str, b: str) -> float:
    """Matches / length of the longer sequence. Stated explicitly because
    'identity' is ambiguous and the gate is a hard threshold."""
    if not a or not b:
        return 0.0
    aa, bb = _align(a, b)
    matches = sum(1 for x, y in zip(aa, bb) if x == y and x != "-")
    return matches / max(len(a), len(b))


def mutation_count(a: str, b: str) -> int:
    """Alignment columns that are not an identical match: substitutions +
    indels. An indel counts as one mutation per residue."""
    if not a or not b:
        return max(len(a), len(b))
    aa, bb = _align(a, b)
    return sum(1 for x, y in zip(aa, bb) if x != y)


# --- numbering -------------------------------------------------------------

def regions(numbering):
    out = {name: "" for name in IMGT_REGIONS}
    for (pos, _ins), aa in numbering:
        if aa == "-":
            continue
        for name, (lo, hi) in IMGT_REGIONS.items():
            if lo <= pos <= hi:
                out[name] += aa
                break
    return out


def number_domains(model, seq: str, max_domains: int = 2):
    """Number `seq`, then re-number whatever ANARCII left over.

    An scFv is one string holding two domains; ANARCII returns the first it
    finds. Anything substantial after query_end is numbered as a second domain,
    so H+L in a single chain is detected rather than silently half-read.
    """
    domains, offset, rest = [], 0, seq
    while rest and len(domains) < max_domains:
        rec = list(model.number([rest]).values())[0]
        if rec.get("numbering") is None or rec.get("chain_type") == "F":
            break
        start, end = rec.get("query_start"), rec.get("query_end")
        if start is None or end is None:
            break
        domains.append({
            "chain_type": rec["chain_type"],
            "score": float(rec["score"]),
            "numbering": rec["numbering"],
            "regions": regions(rec["numbering"]),
            "span": (offset + start, offset + end),
            "seq": rest[start:end + 1],
        })
        tail_start = end + 1
        if len(rest) - tail_start < 80:  # too short to be another V domain
            break
        offset += tail_start
        rest = rest[tail_start:]
    return domains


# --- the gates -------------------------------------------------------------

def gate_g1_format(domains, n_chains):
    """VHH, scFv or Fab. Inferred: one H domain = VHH; H+L in one chain =
    scFv; H+L across two chains = Fab-like."""
    types = [d["chain_type"] for d in domains]
    heavy = sum(1 for t in types if t == "H")
    light = sum(1 for t in types if t in ("K", "L"))
    if not domains:
        return False, "no antibody domain found", ""
    if n_chains == 1 and heavy == 1 and light == 0:
        return True, "", "VHH"
    if n_chains == 1 and heavy == 1 and light == 1:
        return True, "", "scFv"
    if n_chains == 2 and heavy == 1 and light == 1:
        return True, "", "Fab"
    return False, f"not VHH/scFv/Fab (H={heavy} L={light} chains={n_chains})", ""


def gate_g2_cysteines(domains, chains):
    """No unpaired cysteines.

    Sequence alone cannot prove a disulfide, so this is the checkable half:
    the conserved IMGT 23/104 pair must be present in every domain, and any
    further cysteine must have a partner (even count). Reported as such.
    """
    problems = []
    numbered_total = 0
    for d in domains:
        pos_of = {}
        for (pos, _ins), aa in d["numbering"]:
            if aa != "-":
                pos_of.setdefault(pos, aa)
        missing = [p for p in CONSERVED_CYS if pos_of.get(p) != "C"]
        if missing:
            problems.append(
                f"{d['chain_type']}: conserved Cys missing at IMGT "
                f"{','.join(str(p) for p in missing)}")
        cys = sum(1 for aa in pos_of.values() if aa == "C")
        numbered_total += len(d["seq"])
        extra = cys - (len(CONSERVED_CYS) - len(missing))
        if extra % 2:
            problems.append(f"{d['chain_type']}: {extra} extra Cys, odd -> unpaired")
    # cysteines sitting outside any numbered domain (linkers, tags)
    outside = sum(c.count("C") for c in chains) - sum(
        d["seq"].count("C") for d in domains)
    if outside:
        problems.append(f"{outside} Cys outside the numbered domains")
    return (not problems), "; ".join(problems)


def gate_g3_confidence(domains):
    """SELECTION.md G3 is 'immunum numbering confidence >= 0.1'.

    We do not have immunum; ANARCII's score is a different scale entirely
    (good chains land ~20-35, not 0-1), so this is a *proxy*, not the gate.
    It fails only what ANARCII cannot recognise as an antibody at all, which
    is the part of G3 we can honestly check. abguard remains authoritative.
    """
    if not domains:
        return False, "ANARCII found no antibody domain"
    weak = [f"{d['chain_type']}={d['score']:.1f}" for d in domains
            if d["score"] < MIN_ANARCII_SCORE]
    if weak:
        return False, f"ANARCII score below {MIN_ANARCII_SCORE} ({', '.join(weak)})"
    return True, ""


def gate_g4_length(total_aa):
    ok = total_aa <= MAX_TOTAL_AA
    return ok, "" if ok else f"{total_aa} aa > {MAX_TOTAL_AA}"


def gate_o1_cdrh3_identity(cdrh3):
    if not cdrh3:
        return False, "no CDR-H3 found"
    ident = identity(cdrh3, TRAST_CDRH3)
    ok = ident < CDRH3_IDENTITY_MAX
    return ok, "" if ok else f"CDR-H3 identity to trastuzumab {ident:.0%} >= 70%"


def gate_d2_cdrh3_mutations(cdrh3):
    if not cdrh3:
        return False, "no CDR-H3 found"
    n = mutation_count(cdrh3, TRAST_CDRH3)
    ok = n >= MIN_CDRH3_MUTATIONS
    return ok, "" if ok else f"only {n} CDR-H3 mutations from trastuzumab, need >= 2"


def gate_o2_point_mutations(chains, data_points):
    """O2: strictly more than 2 point mutations from ANY provided data point.

    Without the data points this is unchecked, not passed.
    """
    if not data_points:
        return None, UNCHECKED["O2"]
    query = max(chains, key=len)
    worst = min(mutation_count(query, dp) for dp in data_points)
    ok = worst > 2
    return ok, "" if ok else f"only {worst} point mutations from the nearest data point"


# --- per-candidate driver --------------------------------------------------

def normalise_track(raw):
    t = (raw or "").strip().lower().replace("-", "_").replace(" ", "_")
    if t in ("de_novo", "denovo", "novo"):
        return "de_novo"
    if t in ("optimisation", "optimization", "opt"):
        return "optimisation"
    return ""


def evaluate(model, cand, data_points):
    chains = cand["chains"]
    track = cand["track"]
    total_aa = sum(len(c) for c in chains)

    domains = []
    for c in chains:
        # One chain may hold more than one V domain (scFv, tandem constructs);
        # number up to 3 so G2's cysteine accounting is not thrown off by a
        # domain left unread.
        domains.extend(number_domains(model, c, max_domains=3 if len(chains) == 1 else 1))

    heavy = next((d for d in domains if d["chain_type"] == "H"), None)
    cdrh3 = heavy["regions"]["CDR3"] if heavy else ""

    results = {}   # gate -> (True | False | None, detail)
    ok, detail, fmt = gate_g1_format(domains, len(chains))
    results["G1"] = (ok, detail)
    results["G2"] = gate_g2_cysteines(domains, chains)
    results["G3"] = gate_g3_confidence(domains)
    results["G4"] = gate_g4_length(total_aa)

    if track == "de_novo":
        results["D1"] = (None, UNCHECKED["D1"])
        results["D2"] = gate_d2_cdrh3_mutations(cdrh3)
    elif track == "optimisation":
        results["O1"] = gate_o1_cdrh3_identity(cdrh3)
        results["O2"] = gate_o2_point_mutations(chains, data_points)
    else:
        results["TRACK"] = (False, "no track assigned, and §2 requires it "
                                   "before generation")

    failed = [g for g, (ok, _) in results.items() if ok is False]
    unchecked = [g for g, (ok, _) in results.items() if ok is None]
    notes = "; ".join(f"{g}: {d}" for g, (ok, d) in results.items()
                      if ok is not True and d)

    return {
        "design_id": cand["design_id"],
        "job_name": cand["job_name"],
        "contig": cand["contig"],
        "generated_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "gate_pass": "" if failed else ("PASS" if not unchecked else "PASS*"),
        "gate_failed": ",".join(failed),
        "length_aa": total_aa,
        "cdr3_seq": cdrh3,
        "ipTM": "", "pLDDT": "", "motif_rmsd": "", "tnp_flags": "",
        "epitope_ok": "", "tier": "", "portfolio_group": "", "submitted": "",
        "score_source_tool": f"gates.py (anarcii {anarcii.__version__})",
        "score_source_version": "gates.py/1",
        "notes": notes,
        "_format": fmt,
        "_track": track or "UNASSIGNED",
        "_unchecked": unchecked,
        "_domains": [f"{d['chain_type']}:{d['score']:.1f}" for d in domains],
    }


# --- input / output --------------------------------------------------------

def split_chains(raw):
    parts = [p.strip().upper() for p in raw.replace("/", ":").split(":")]
    return [p for p in parts if p]


def load_fasta(path):
    cands, hdr, buf = [], None, []

    def flush():
        if hdr is None:
            return
        fields = hdr.split()
        design_id = fields[0]
        kv = dict(f.split("=", 1) for f in fields[1:] if "=" in f)
        cands.append({
            "design_id": design_id,
            "chains": split_chains("".join(buf)),
            "track": normalise_track(kv.get("track")),
            "job_name": kv.get("job", ""),
            "contig": kv.get("contig", ""),
        })

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            flush()
            hdr, buf = line[1:], []
        elif line:
            buf.append(line)
    flush()
    return cands


def load_csv(path):
    cands = []
    with open(path, newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh), 1):
            row = {(k or "").strip().lower(): (v or "").strip()
                   for k, v in row.items()}
            seq = row.get("sequence") or row.get("heavy") or ""
            chains = split_chains(seq)
            light = row.get("light") or ""
            if light:
                chains += split_chains(light)
            if not chains:
                continue
            cands.append({
                "design_id": row.get("design_id") or f"row{i}",
                "chains": chains,
                "track": normalise_track(row.get("track")),
                "job_name": row.get("job_name", ""),
                "contig": row.get("contig", ""),
            })
    return cands


def append_scores(rows):
    existing_header = None
    if SCORES_CSV.exists():
        with open(SCORES_CSV, newline="") as fh:
            existing_header = next(csv.reader(fh), None)
    if existing_header and existing_header != COLUMNS:
        raise SystemExit(
            f"{SCORES_CSV.name} header does not match the expected columns.\n"
            f"  file:     {existing_header}\n  expected: {COLUMNS}\n"
            "Refusing to append and corrupt the provenance table.")
    new = not SCORES_CSV.exists() or existing_header is None
    with open(SCORES_CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nappended {len(rows)} row(s) to {SCORES_CSV}")


def report(rows):
    print(f"\n{'design_id':<22} {'track':<13} {'fmt':<5} {'aa':>4} "
          f"{'verdict':<8} detail")
    print("-" * 100)
    for r in rows:
        verdict = "FAIL" if r["gate_failed"] else r["gate_pass"]
        detail = r["gate_failed"] or (
            f"unchecked: {','.join(r['_unchecked'])}" if r["_unchecked"] else "")
        print(f"{r['design_id']:<22} {r['_track']:<13} {r['_format'] or '-':<5} "
              f"{r['length_aa']:>4} {verdict:<8} {detail}")
        if r["notes"]:
            for part in r["notes"].split("; "):
                print(f"{'':<22} - {part}")

    passed = [r for r in rows if not r["gate_failed"]]
    print(f"\n{len(passed)}/{len(rows)} passed the checkable gates.")
    print("PASS* = passed everything checkable, with gates left unchecked.")
    print("abguard is authoritative and runs on submission (SELECTION.md §2).")


def self_test(model):
    """Trastuzumab through the gates both ways.

    As de novo it should clear G1-G4 and fail D2 (zero CDR-H3 mutations from
    itself); as optimisation it should fail O1 (100% CDR-H3 identity). If
    those two come out the other way the gate logic is inverted.
    """
    print("--- self-test: trastuzumab VH against both tracks ---")
    rows = []
    for track in ("de_novo", "optimisation"):
        rows.append(evaluate(model, {
            "design_id": f"trastuzumab-VH[{track}]",
            "chains": [TRAST_VH], "track": track,
            "job_name": "self-test", "contig": "",
        }, data_points=[TRAST_VH]))
    report(rows)

    checks = [
        ("de_novo fails D2", "D2" in rows[0]["gate_failed"]),
        ("optimisation fails O1", "O1" in rows[1]["gate_failed"]),
        ("optimisation fails O2", "O2" in rows[1]["gate_failed"]),
        ("G1 typed it VHH", rows[0]["_format"] == "VHH"),
        ("G2/G3/G4 all pass",
         not {"G2", "G3", "G4"} & set(rows[0]["gate_failed"].split(","))),
        ("CDR-H3 == SRWGGDGFYAMDY", rows[0]["cdr3_seq"] == TRAST_CDRH3),
    ]
    print()
    bad = 0
    for name, ok in checks:
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}")
        bad += not ok
    print("\nself-test:", "all checks passed" if not bad else f"{bad} FAILED")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", nargs="?", help="CSV or FASTA of candidates")
    ap.add_argument("--write", action="store_true",
                    help=f"append the rows to {SCORES_CSV.name}")
    ap.add_argument("--data-points", metavar="FILE",
                    help="FASTA of provided data points, for gate O2")
    ap.add_argument("--self-test", action="store_true",
                    help="run the gates on trastuzumab and verify the logic")
    args = ap.parse_args()

    if not args.candidates and not args.self_test:
        ap.error("give a candidates file, or --self-test")

    model = Anarcii(seq_type="antibody", mode="accuracy", cpu=True, verbose=False)

    if args.self_test:
        rc = self_test(model)
        if not args.candidates:
            return rc

    data_points = []
    if args.data_points:
        data_points = ["".join(c["chains"]) for c in load_fasta(args.data_points)]
        print(f"loaded {len(data_points)} data point(s) for gate O2")

    path = Path(args.candidates)
    if not path.exists():
        raise SystemExit(f"no such file: {path}")
    cands = (load_fasta(path) if path.suffix.lower() in (".fa", ".fasta", ".faa")
             else load_csv(path))
    if not cands:
        raise SystemExit(f"no candidates parsed from {path}")
    print(f"{len(cands)} candidate(s) from {path}")

    rows = [evaluate(model, c, data_points) for c in cands]
    report(rows)
    if args.write:
        append_scores(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
