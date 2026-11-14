#!/usr/bin/env python3
"""Pair every submitted design with the trastuzumab light chain as an scFv.

Why this exists: submission.csv carries a bare 120-aa VH and no light chain at
all. The brief requires VHH, scFv or Fab, and a lone VH domain from a paired
antibody is none of those -- gates.py:gate_g1_format classifies one H domain
with no L as "VHH", so our own pre-check passes it, and abguard's type-aware
comparison ("a VHH design is compared with known VHHs or heavy chains") may
well pass it too. The gate is not the problem. Epitope preservation is: it is
one of the three properties Adaptyv measures against the measured trastuzumab
baseline, trastuzumab's HER2 epitope is contacted by both VH and VL CDRs, and a
VH-only construct cannot present it. SELECTION.md §1 commits every design to a
credible attempt at all three properties simultaneously; twenty VH-only designs
forfeit one before anything is measured.

The designs themselves are unchanged. Only CDRH3 was ever varied, the VH
framework is trastuzumab's own, and the light chain is trastuzumab's own. So
pairing restores the baseline architecture without touching the ranking: rank
order in submission_scfv.csv is exactly the rank order of submission.csv.

VL provenance (SELECTION.md §8, "verified, not assumed"): PDB 1N8Z chain A,
"Herceptin Fab (antibody) - light chain", fetched from RCSB. The variable
domain boundary is set by ANARCII numbering, not by hand -- numbering the full
chain returns chain_type K over IMGT 1-127 for the first 107 residues, and the
kappa constant domain that follows is unnumbered and dropped. Re-derive with
--verify-vl, which refuses to write anything unless ANARCII agrees.

Linker: (GGGGS)x3, 15 aa. Total 120 + 15 + 107 = 242 aa, inside the brief's
250-aa cap ("whole length of a single scFv or VHH, or combined for VH and VL")
with 8 aa of margin. Orientation is VH-linker-VL.

Cysteines stay paired: 2 in VH, 2 in VL, 4 per scFv, each domain's own pair.

Run:  python3 make_scfv.py [--verify-vl]
Out:  submission_scfv.csv   full submission, one scFv per row
      scfv_gates_input.csv  design_id,sequence,track for gates.py
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# PDB 1N8Z chain A, variable domain only (IMGT 1-127 per ANARCII; the kappa
# constant domain beginning RTVAAPSVFIFPP... is excluded).
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIY"
    "SASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

LINKER = "GGGGS" * 3
MAX_TOTAL = 250


def verify_vl():
    """Re-derive the VL domain boundary with ANARCII and check our constant.

    Keeps smoke.py's import fallback: the anarcii import path has differed
    between versions.
    """
    try:
        from anarcii import Anarcii
    except ImportError:
        from anarcii.pipeline import Anarcii

    # Full 1N8Z chain A, variable + kappa constant, as RCSB returns it.
    full_light = TRASTUZUMAB_VL + (
        "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDS"
        "KDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    )
    model = Anarcii(seq_type="antibody", mode="accuracy", cpu=True, verbose=False)
    result = model.number([full_light])
    record = next(iter(result.values()))

    numbered = [(n, r) for n, r in record["numbering"] if r != "-"]
    domain = "".join(r for _, r in numbered)

    problems = []
    if record["chain_type"] != "K":
        problems.append(f"chain_type is {record['chain_type']!r}, expected 'K'")
    if record["score"] < 0.1:
        problems.append(f"numbering confidence {record['score']:.4f} below 0.1")
    if domain != TRASTUZUMAB_VL:
        problems.append(
            f"ANARCII numbers a {len(domain)}-aa variable domain that differs "
            f"from TRASTUZUMAB_VL ({len(TRASTUZUMAB_VL)} aa)"
        )
    if problems:
        sys.exit("VL verification failed:\n  " + "\n  ".join(problems))

    print(
        f"VL verified: chain_type {record['chain_type']}, "
        f"confidence {record['score']:.4f}, "
        f"IMGT {numbered[0][0][0]}-{numbered[-1][0][0]}, "
        f"{len(domain)} aa, {domain.count('C')} cysteines"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=os.path.join(HERE, "submission_o2.csv"),
                    help="ranked list to pair; defaults to the O2-filtered one")
    ap.add_argument("--out", default=os.path.join(HERE, "submission_scfv.csv"))
    ap.add_argument("--verify-vl", action="store_true",
                    help="re-derive the VL boundary with ANARCII before writing")
    args = ap.parse_args()

    if args.verify_vl:
        verify_vl()

    src = args.source
    if not os.path.exists(src):
        sys.exit(f"{src} not found")

    rows = list(csv.DictReader(open(src, newline="")))
    if not rows:
        sys.exit(f"{src} is empty")

    # design_id does not name a molecule: scores.csv holds 234 rows under only
    # 168 distinct ids, because the `novel` and `exploratory` runs reuse ids for
    # different CDRH3s (her2_cdrh3_98 is WGGGRQYTYD in one and TTAAVGYAHD in the
    # other). 7 of the 20 submitted ids are ambiguous that way. submission.csv
    # carries no `run` column, so it cannot be traced back to a job -- against
    # §8's "weights, seeds and settings are recorded per job". The pair
    # (design_id, cdrh3) does resolve all 7 uniquely, so join on that and carry
    # run/weights/seed through to the submission.
    scores_path = os.path.join(HERE, "scores.csv")
    provenance = {}
    if os.path.exists(scores_path):
        for r in csv.DictReader(open(scores_path)):
            provenance.setdefault((r["design_id"], r["cdrh3"]), r)

    out_rows = []
    for row in rows:
        vh = row["vh"].strip()
        scfv = vh + LINKER + TRASTUZUMAB_VL

        # Fail loudly rather than writing a design that cannot be submitted.
        if len(scfv) > MAX_TOTAL:
            sys.exit(
                f"{row['design_id']}: scFv is {len(scfv)} aa, over the "
                f"{MAX_TOTAL}-aa cap"
            )
        if scfv.count("C") % 2 != 0:
            sys.exit(f"{row['design_id']}: odd cysteine count in scFv")
        if row["cdrh3"] not in vh:
            sys.exit(f"{row['design_id']}: cdrh3 not found in its own vh")

        src_row = provenance.get((row["design_id"], row["cdrh3"]))
        if provenance and src_row is None:
            sys.exit(
                f"{row['design_id']} / {row['cdrh3']} is not in scores.csv -- "
                f"cannot establish which run produced it"
            )

        out_rows.append(
            {
                "rank": row["rank"],
                "design_id": row["design_id"],
                "run": (src_row or {}).get("run", ""),
                "weights": (src_row or {}).get("weights", ""),
                "seed": (src_row or {}).get("seed", ""),
                "track": row["track"],
                "format": "scFv",
                "cdrh3": row["cdrh3"],
                "vh": vh,
                "linker": LINKER,
                "vl": TRASTUZUMAB_VL,
                "scfv": scfv,
                "scfv_len": len(scfv),
                "n_cys": scfv.count("C"),
                "vl_source": "PDB 1N8Z chain A, variable domain (ANARCII IMGT 1-127)",
            }
        )

    dst = args.out
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)

    gates_dst = os.path.join(HERE, "scfv_gates_input.csv")
    with open(gates_dst, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["design_id", "sequence", "track"])
        for row in out_rows:
            writer.writerow([row["design_id"], row["scfv"], row["track"]])

    lengths = {r["scfv_len"] for r in out_rows}
    print(f"wrote {len(out_rows)} scFvs to {os.path.basename(dst)}")
    print(f"  lengths {sorted(lengths)} (cap {MAX_TOTAL})")
    print(f"  rank order preserved from {os.path.basename(src)}")
    print(f"wrote scfv_gates_input.csv for gates.py")


if __name__ == "__main__":
    main()
