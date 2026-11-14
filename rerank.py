#!/usr/bin/env python3
"""Rebuild the top 20 from the designs that actually pass gate O2.

What this is, precisely: a **filter**, not a new ranking. `check_o2.py` showed
13 of the 20 submitted designs fail O2 (>2 point mutations from any provided
data point), rank 1 among them, so 13 of 20 slots were going to be rejected on
submission. 118 of the 234 ranked designs pass.

What this does NOT do, deliberately. §4 forbids moving a threshold once ranking
has begun, and §1's whole point is that selection is a method and not a story
told afterwards. So:

  * no score is recomputed and no weight is touched;
  * the order is `scores.csv`'s existing rank order, which already has the §7
    tie-break applied (rounded composite, then generator prior, then
    |mutations-6|, then design id);
  * the only operation is dropping the rows that fail a §2 hard gate, then
    taking the first 20 that remain and renumbering 1..20.

A hard gate is not a ranking criterion -- it decides eligibility, and §2 always
said so ("ranking -- gate survivors only", §3). Designs that fail O2 were never
eligible; they were ranked because O2 was unchecked. Removing them applies the
pre-registered rule rather than departing from it. Relative order among the
survivors is untouched: if A outranked B before, it still does.

`design_id` does not identify a molecule (234 rows, 168 distinct ids; the
`novel` and `exploratory` runs reuse ids), so every join is on
(design_id, cdrh3) and `run` is carried through to the output.

Run:  python3 rerank.py [--o2-report o2_report_all234.csv] [-n 20]
Out:  submission_o2.csv   the new ranked list, same columns as submission.csv
                          plus `run`, ready for make_scfv.py --source
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(HERE, "scores.csv"))
    ap.add_argument("--o2-report",
                    default=os.path.join(HERE, "o2_report_all234.csv"),
                    help="check_o2.py output covering all ranked designs")
    ap.add_argument("-n", type=int, default=20,
                    help="list length (§5: 16 assuming a win, plus 4 spares)")
    args = ap.parse_args()

    for p in (args.scores, args.o2_report):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    scores = list(csv.DictReader(open(args.scores, newline="")))
    o2 = list(csv.DictReader(open(args.o2_report, newline="")))

    # scores.csv must already be in rank order for "take the first N" to mean
    # "keep the existing order". Verify rather than assume.
    ranks = [int(r["rank"]) for r in scores]
    if ranks != sorted(ranks):
        sys.exit("scores.csv is not in rank order -- refusing to 'preserve' an "
                 "order that is not there")

    verdict = {(r["design_id"], r["cdrh3"]): r for r in o2}
    unjudged = [r for r in scores if (r["design_id"], r["cdrh3"]) not in verdict]
    if unjudged:
        sys.exit(
            f"{len(unjudged)} ranked design(s) have no O2 verdict "
            f"(e.g. {unjudged[0]['design_id']}/{unjudged[0]['cdrh3']}) -- "
            f"re-run check_o2.py --designs scores.csv first"
        )

    survivors = [r for r in scores
                 if verdict[(r["design_id"], r["cdrh3"])]["o2_verdict"] == "PASS"]
    print(f"{len(scores)} ranked designs, {len(survivors)} pass O2, "
          f"{len(scores) - len(survivors)} dropped")
    if len(survivors) < args.n:
        sys.exit(f"only {len(survivors)} survivors, need {args.n}")

    chosen = survivors[:args.n]

    out = []
    for new_rank, r in enumerate(chosen, start=1):
        v = verdict[(r["design_id"], r["cdrh3"])]
        out.append({
            "rank": new_rank,
            "prev_rank": r["rank"],
            "design_id": r["design_id"],
            "run": r["run"],
            "track": r["track"],
            "cdrh3": r["cdrh3"],
            "vh": r["vh"],
            "composite_score": r["composite_score"],
            "min_mut_from_data_point": v["min_mut_from_data_point"],
        })

    dst = os.path.join(HERE, "submission_o2.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    # Relative order must be unchanged: previous ranks strictly increasing.
    prev = [int(r["prev_rank"]) for r in out]
    if prev != sorted(prev):
        sys.exit("relative order changed -- this should be a filter only")

    kept = {(r["design_id"], r["cdrh3"]) for r in chosen}
    old_sub = os.path.join(HERE, "submission.csv")
    if os.path.exists(old_sub):
        old = list(csv.DictReader(open(old_sub, newline="")))
        held = [r for r in old if (r["design_id"], r["cdrh3"]) in kept]
        print(f"of the previously submitted 20, {len(held)} are retained, "
              f"{len(old) - len(held)} dropped")

    print()
    print(f"{'new':>3} {'was':>4}  {'design_id':<18} {'run':<12} {'CDRH3':<11} "
          f"{'score':>7} {'O2min':>5}")
    for r in out:
        print(f"{r['rank']:>3} {r['prev_rank']:>4}  {r['design_id']:<18} "
              f"{r['run']:<12} {r['cdrh3']:<11} {r['composite_score']:>7} "
              f"{r['min_mut_from_data_point']:>5}")
    print()
    print(f"wrote submission_o2.csv (top {args.n}, order preserved from scores.csv)")


if __name__ == "__main__":
    main()
