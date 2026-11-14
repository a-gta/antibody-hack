#!/usr/bin/env python3
"""Write the co-folding input for the top N designs of the current ranking.

Why this exists: results/vm/refold-shortlist.csv was written at 15:30, before
scores.csv existed at 16:24, and it overlaps the submitted top 20 by exactly
one design (her2_cdrh3_110, rank 9). Its other 23 entries rank 25-233 and are
not being submitted. So the fold budget was spent on designs the ranking does
not select, and 19 of the 20 submitted designs have never been folded.

This regenerates the fold input FROM the ranking, so the designs we fold are
the designs we submit. Folding these gives min_ipae / iptm / binder_plddt for
the shortlist that actually decides the submission -- i.e. SELECTION.md's
pre-registered R1-R3, applied where they change the answer.

N defaults to 24: the top 20 plus 4 spares, matching the size of the shortlist
this replaces.

Column order is copied from results/vm/refold-shortlist.csv so the existing
submit path takes this file unchanged.

Run:  python3 make_fold_input.py [N]
Out:  results/vm/refold-top<N>.csv
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ["design", "cdrh3", "n_mut_vs_trastuzumab", "hotspot_obedience", "vh"]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 24

    scores = os.path.join(HERE, "scores.csv")
    if not os.path.exists(scores):
        sys.exit("scores.csv not found -- run score_designs.py first")

    rows = list(csv.DictReader(open(scores)))
    # scores.csv is written in rank order, but sort on the rank column rather
    # than trusting file order.
    rows.sort(key=lambda r: int(r["rank"]))
    top = rows[:n]
    if len(top) < n:
        print("warning: only %d ranked designs available" % len(top))

    out = os.path.join(HERE, "results", "vm", "refold-top%d.csv" % n)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in top:
            w.writerow({
                "design": r["design_id"],
                "cdrh3": r["cdrh3"],
                "n_mut_vs_trastuzumab": r["n_mut_vs_trastuzumab"],
                "hotspot_obedience": r["hotspot_obedience"],
                "vh": r["vh"],
            })

    print("wrote %s (%d designs)" % (os.path.relpath(out, HERE), len(top)))
    print("%4s %18s %12s %4s" % ("rank", "design", "cdrh3", "mut"))
    print("-" * 44)
    for r in top:
        print("%4s %18s %12s %4s" % (
            r["rank"], r["design_id"], r["cdrh3"], r["n_mut_vs_trastuzumab"]))

    # The overlap with the file this replaces, stated rather than assumed.
    old = os.path.join(HERE, "results", "vm", "refold-shortlist.csv")
    if os.path.exists(old):
        prev = {r["design"] for r in csv.DictReader(open(old))}
        new = {r["design_id"] for r in top}
        print("\noverlap with refold-shortlist.csv: %d of %d"
              % (len(prev & new), len(new)))
        print("already folded, no recompute needed: %s"
              % (sorted(prev & new) or "none"))


if __name__ == "__main__":
    main()
