#!/usr/bin/env python3
"""Gate O2, actually checked: point mutations from the nearest provided data point.

SELECTION.md §2 gate O2 requires **strictly more than 2** point mutations from
any provided data point. gates.py reports it `unchecked` because the provided
data was never pulled off the VM, and §7 logs it as the one remaining gate that
can still cost submission slots. This checks it.

Data, from the brief's own citation (p. 4, "Published in Chinery et al. ·
Zenodo"), DOI 10.5281/zenodo.22310365, BSD-3:

  her2_aff_large train/test/val.csv   524,346 CDRH3 variants, FACS, three
                                      affinity classes. 367,042 + 78,652 +
                                      78,652 = 524,346, matching the count the
                                      brief states exactly -- which is how we
                                      know this is the set O2 refers to.
  spr_submitted.csv                   768 designs from the same study, with
                                      paired H/L and measured SPR labels.

Both give CDRH3 as a 10-mer in trastuzumab's own VH framework, and our designs
vary nothing but CDRH3 in that same framework. So Hamming distance on the
10-mer **is** the full-VH point-mutation count: no alignment needed, and no
indels are possible. The comparison is exact, not a proxy.

Distances are computed over every data point, not a sample, with numpy.

Run:  python3 check_o2.py --variants DIR --spr FILE [--designs FILE]
Out:  o2_report.csv          per design: nearest data point and verdict
      o2_data_points.fasta   the nearest data points, as full VH, so
                             `gates.py --data-points` reproduces this verdict
                             without re-reading 524k rows
"""
import argparse
import csv
import os
import sys

import numpy as np

# The same counter gates.py applies, so O2 is not re-implemented here with
# subtly different semantics (substitutions + one per indel residue).
from gates import mutation_count

HERE = os.path.dirname(os.path.abspath(__file__))

# PDB 1N8Z chain B / spr_submitted.csv IgG_1 "Positive Control".
TRAST_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
assert len(TRAST_VH) == 120, len(TRAST_VH)
TRAST_CDRH3 = "WGGDGFYAMD"
O2_MIN_MUTATIONS = 2  # strictly more than this


def encode(seqs, width):
    """Pack equal-length sequences into a uint8 matrix for exact Hamming."""
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)
    return arr.reshape(len(seqs), width)


def load_variants(paths):
    """CDRH3 10-mers from the her2_aff_large splits."""
    seqs = []
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                s = (row.get("seq") or "").strip().upper()
                if s:
                    seqs.append(s)
    return seqs


def load_spr(path):
    """CDRH3 10-mers and full VH from spr_submitted.csv."""
    h3s, vhs = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            h3 = (row.get("H3") or "").strip().upper()
            vh = (row.get("H") or "").strip().upper()
            if h3:
                h3s.append(h3)
            if vh:
                vhs.append(vh)
    return h3s, vhs


def graft(cdrh3):
    """Put a CDRH3 10-mer back into trastuzumab's VH framework."""
    i = TRAST_VH.find(TRAST_CDRH3)
    if i < 0:
        sys.exit("trastuzumab CDRH3 not found in TRAST_VH -- framework constant is wrong")
    return TRAST_VH[:i] + cdrh3 + TRAST_VH[i + len(TRAST_CDRH3):]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True,
                    help="directory holding her2_aff_large train/test/val.csv")
    ap.add_argument("--spr", required=True, help="spr_submitted.csv")
    ap.add_argument("--designs", default=os.path.join(HERE, "submission_scfv.csv"))
    args = ap.parse_args()

    vpaths = [os.path.join(args.variants, n)
              for n in ("train.csv", "val.csv", "test.csv")]
    missing = [p for p in vpaths if not os.path.exists(p)]
    if missing:
        sys.exit("missing variant split(s): " + ", ".join(missing))

    variants = load_variants(vpaths)
    spr_h3s, spr_vhs = load_spr(args.spr)
    print(f"variants: {len(variants)} rows, {len(set(variants))} distinct")
    print(f"spr_submitted: {len(spr_h3s)} CDRH3, {len(set(spr_h3s))} distinct")

    # 10-mers get exact Hamming via numpy. Anything else (spr_submitted holds a
    # few 9-mers, i.e. a deletion) is an indel case, so it goes through
    # gates.py's own alignment-based mutation_count rather than a second
    # implementation of the same rule.
    pool = sorted(set(variants) | set(spr_h3s))
    pool10 = [s for s in pool if len(s) == len(TRAST_CDRH3)]
    pool_indel = [s for s in pool if len(s) != len(TRAST_CDRH3)]
    print(f"data points: {len(pool)} distinct "
          f"({len(pool10)} same-length, {len(pool_indel)} indel)")

    # Any spr VH that differs from trastuzumab OUTSIDE CDRH3 would break the
    # "CDRH3 distance == VH distance" equivalence this check rests on. Verify.
    framework_mismatch = []
    for vh in set(spr_vhs):
        if len(vh) != len(TRAST_VH):
            framework_mismatch.append(("length", len(vh)))
            continue
        i = TRAST_VH.find(TRAST_CDRH3)
        if vh[:i] != TRAST_VH[:i] or vh[i + 10:] != TRAST_VH[i + 10:]:
            framework_mismatch.append(("framework", vh))
    if framework_mismatch:
        # These cannot lower any design's minimum. Our designs carry
        # trastuzumab's framework exactly, so for such a data point the full-VH
        # distance is (framework differences) + (CDRH3 differences), which is
        # never below the CDRH3 difference alone -- and that data point's own
        # CDRH3 is already in the pool as a 10-mer. O2 is a minimum, so the
        # minimum stands. Reported, not silently dropped.
        print(f"NOTE: {len(framework_mismatch)} spr VH differ from trastuzumab "
              f"outside CDRH3. Their full-VH distance is framework + CDRH3 "
              f"differences, never less than CDRH3 alone, and their CDRH3s are "
              f"already in the pool -- so they cannot lower any minimum below.")

    pool_arr = encode(pool10, len(TRAST_CDRH3))

    designs = list(csv.DictReader(open(args.designs, newline="")))
    if not designs:
        sys.exit(f"no designs in {args.designs}")

    rows, failures = [], []
    for d in designs:
        cdrh3 = d["cdrh3"].strip().upper()
        if len(cdrh3) != len(TRAST_CDRH3):
            sys.exit(f"{d['design_id']}: CDRH3 is {len(cdrh3)} aa, expected 10")
        q = np.frombuffer(cdrh3.encode("ascii"), dtype=np.uint8)
        dist = (pool_arr != q).sum(axis=1)
        best = int(dist.min())
        nearest = pool10[int(dist.argmin())]
        n_at_best = int((dist == best).sum())

        # The indel data points, through gates.py's own counter.
        for dp in pool_indel:
            n = mutation_count(cdrh3, dp)
            if n < best:
                best, nearest, n_at_best = n, dp, 1
            elif n == best:
                n_at_best += 1

        passes = best > O2_MIN_MUTATIONS

        rows.append({
            "rank": d["rank"],
            "design_id": d["design_id"],
            "cdrh3": cdrh3,
            "mut_from_trastuzumab": int(sum(
                a != b for a, b in zip(cdrh3, TRAST_CDRH3))),
            "min_mut_from_data_point": best,
            "nearest_data_point": nearest,
            "n_data_points_at_min": n_at_best,
            "o2_verdict": "PASS" if passes else "FAIL",
        })
        if not passes:
            failures.append((d["rank"], d["design_id"], best, nearest))

    dst = os.path.join(HERE, "o2_report.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # Nearest data points only, grafted to full VH: O2 is a minimum, so the
    # nearest point decides it. Feeding gates.py these reproduces the verdict.
    fasta = os.path.join(HERE, "o2_data_points.fasta")
    with open(fasta, "w") as fh:
        for cdrh3 in sorted({r["nearest_data_point"] for r in rows}):
            fh.write(f">data_point_{cdrh3}\n{graft(cdrh3)}\n")

    print()
    print(f"{'rank':>4}  {'design_id':<18} {'CDRH3':<11} {'vs Tz':>5} "
          f"{'min':>4} {'nearest':<11} verdict")
    for r in rows:
        print(f"{r['rank']:>4}  {r['design_id']:<18} {r['cdrh3']:<11} "
              f"{r['mut_from_trastuzumab']:>5} {r['min_mut_from_data_point']:>4} "
              f"{r['nearest_data_point']:<11} {r['o2_verdict']}")

    print()
    n_fail = len(failures)
    print(f"O2 (> {O2_MIN_MUTATIONS} mutations from any provided data point): "
          f"{len(rows) - n_fail}/{len(rows)} pass, {n_fail} fail")
    print(f"compared against all {len(pool)} distinct provided CDRH3 data points")
    print(f"wrote o2_report.csv and o2_data_points.fasta")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
