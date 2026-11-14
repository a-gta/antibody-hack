#!/usr/bin/env python3
"""Add a diversity constraint to the R1 ranking: D = 3, fixed before running.

Ranking by §3's R1 alone converged on one structural family -- 16 of 20
submitted CDRH3s sharing `ST..V[DGS]YA[FMH]D`, 39 pairs within 2 mutations of
each other, the closest pair differing at a single position. Adaptyv measures
20 designs; if that family does not bind, most of the slots fail together.
Spreading that risk is a decision, not a correction, so it is pre-registered
here and logged in SELECTION.md §7 as a departure from §3-§5.

THE RULE, stated before it is run:

    Walk the R1 order from the top. Admit a design only if its CDRH3 differs
    from every already-admitted design at D or more positions. Stop at 20.

    D = 3.

WHY D = 3, on grounds that do not depend on the outcome. This is the one number
that could be tuned to flatter the result, so it is taken from the challenge's
own definition of a distinct molecule rather than chosen: gate O2 requires
strictly more than 2 point mutations from any provided data point, i.e. at
least 3, and O1 sets the same bar against trastuzumab via 70% identity over a
10-mer. If 3 mutations is the distance at which the organisers consider a loop
a different sequence from a known one, it is the distance at which two of our
loops are different from each other. No other value was tried, and the script
is written to fail rather than relax D if 20 designs cannot be found -- relaxing
it to fill the list is exactly the tuning this reasoning is meant to prevent.

WHAT THIS DOES NOT DO:

  * It does not reorder. Admission walks down R1, so the submitted designs stay
    in R1 order and the assertion below proves it: `min_ipae` is monotonic and
    previous ranks strictly increase. A design is dropped, never promoted past a
    better-scoring one.
  * It does not change eligibility or any gate. The pool is the same 62 gate
    survivors `rerank_fold.py` ranks.
  * It does not claim the diverse 20 are better binders. By R1 they are worse:
    each admitted design displaces a higher-scoring near-duplicate. The claim is
    about the variance of the outcome, not its mean, and §7 says so.

Run:  python3 rerank_diverse.py [-n 20] [--out submission_diverse.csv]
"""
import argparse
import csv
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Fixed before the script was first run. See the rationale above.
DIVERSITY_D = 3


def hamming(a, b):
    if len(a) != len(b):
        # Every design here is a 10-mer in trastuzumab's framework; a length
        # difference would mean the pool is not what this rule assumes.
        sys.exit(f"CDRH3 length mismatch: {a!r} vs {b!r}")
    return sum(x != y for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked", default=os.path.join(HERE, "pool_fold_ranked.csv"),
                    help="full R1-ordered pool from rerank_fold.py -n <all>")
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(HERE, "submission_diverse.csv"))
    args = ap.parse_args()

    if not os.path.exists(args.ranked):
        sys.exit(f"missing {args.ranked} -- generate it with:\n"
                 f"  python3 rerank_fold.py -n 62 --out pool_fold_ranked.csv")
    pool = list(csv.DictReader(open(args.ranked, newline="")))

    ipae = [float(r["min_ipae"]) for r in pool]
    if ipae != sorted(ipae):
        sys.exit("pool is not in R1 order -- refusing to 'walk' an order that "
                 "is not there")
    print(f"pool: {len(pool)} gate survivors in R1 order, "
          f"min_ipAE {ipae[0]:.4f}-{ipae[-1]:.4f}")
    print(f"rule: admit if >= {DIVERSITY_D} mutations from every design already "
          f"admitted; walk R1 from the top; stop at {args.n}")

    chosen, skipped = [], []
    for r in pool:
        if len(chosen) == args.n:
            break
        clash = next((c for c in chosen
                      if hamming(r["cdrh3"], c["cdrh3"]) < DIVERSITY_D), None)
        if clash:
            skipped.append((r, clash, hamming(r["cdrh3"], clash["cdrh3"])))
        else:
            chosen.append(r)

    print(f"\nadmitted {len(chosen)}, skipped {len(skipped)} as near-duplicates")
    if len(chosen) < args.n:
        sys.exit(
            f"\nFAILED: only {len(chosen)} of {args.n} designs are mutually "
            f">= {DIVERSITY_D} mutations apart. D is not relaxed to fill the "
            f"list -- that would be tuning the rule to the data. Either submit "
            f"{len(chosen)}, or fold more of the pool to widen it."
        )

    out = []
    for rank, r in enumerate(chosen, start=1):
        row = {"rank": rank, "r1_rank": r["rank"]}
        row.update({k: v for k, v in r.items() if k != "rank"})
        out.append(row)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    # Proof this filtered rather than reordered.
    r1 = [int(r["r1_rank"]) for r in out]
    if r1 != sorted(r1):
        sys.exit("R1 order not preserved -- this must be a filter only")
    ip = [float(r["min_ipae"]) for r in out]
    if ip != sorted(ip):
        sys.exit("output not monotonic in min_ipae -- R1 was not preserved")

    # The rule's own guarantee, checked on the output rather than trusted.
    worst = min(hamming(a["cdrh3"], b["cdrh3"])
                for a, b in itertools.combinations(out, 2))
    assert worst >= DIVERSITY_D, worst

    print()
    print(f"{'new':>3} {'R1':>3}  {'design_id':<16} {'CDRH3':<11} "
          f"{'min_ipAE':>8} {'ipTM':>6} {'apo_holo':>8}")
    for r in out:
        print(f"{r['rank']:>3} {r['r1_rank']:>3}  {r['design_id']:<16} "
              f"{r['cdrh3']:<11} {r['min_ipae']:>8} {r['iptm']:>6} "
              f"{r['apo_holo_rmsd']:>8}")

    dists = [hamming(a["cdrh3"], b["cdrh3"])
             for a, b in itertools.combinations(out, 2)]
    print(f"\ndiversity: closest pair {worst} mutations, mean pairwise "
          f"{sum(dists)/len(dists):.2f}, pairs within 2: "
          f"{sum(1 for d in dists if d <= 2)}")
    print(f"cost in R1: admitting these 20 reaches R1 rank {out[-1]['r1_rank']} "
          f"instead of {args.n}; min_ipAE range {ip[0]:.4f}-{ip[-1]:.4f}")
    if skipped:
        print(f"\nskipped as near-duplicates (the price of the rule):")
        for r, c, d in skipped[:10]:
            print(f"  R1 {r['rank']:>3} {r['design_id']:<16} {r['cdrh3']}  "
                  f"{d} mutation(s) from {c['design_id']} ({c['cdrh3']})")
        if len(skipped) > 10:
            print(f"  ... and {len(skipped) - 10} more")
    print(f"\nwrote {os.path.basename(args.out)}")


if __name__ == "__main__":
    main()
