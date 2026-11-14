#!/usr/bin/env python3
"""Rank the submission on the OpenFold3 refolding scores, as SELECTION.md §3 says.

This replaces `score_designs.py`'s sequence-and-generator composite with the
ranking §3 specified in the first place:

    R1  min_ipae        primary ordering          (lower is better)
    R2  iptm            tiers, per §4
    R3  binder_plddt    tiers, per §4

§7's first entry abandoned R1-R3 because every ranked row had those columns
empty -- only 2 of 400 designs had been folded before the deadline. That is no
longer true. `results/vm/novel-200-of3ob-refold-metrics.csv`, recovered off the
VM before it was shut down, holds OpenFold3 scores for 142 of the 200 novel
seed-50 designs: one row each, one model (`of3ob`), one folding seed (42), one
diffusion sample, every row `predicted`.

So this is **not a new ranking rule**. It is the pre-registered one, applied now
that its inputs exist, and it supersedes the composite rather than adjusting it.
§4 forbids moving a threshold once ranking has begun; restoring the criterion
the file named before any candidate was ranked is the opposite of that.

What it deliberately does NOT do:

  * It does not rank on `apo_holo_rmsd`, though that column is populated and
    informative (3-10 Å across the pool: these binders move when the target is
    removed). §3 names no such criterion, and adding one now, after seeing the
    numbers, is exactly the story-told-afterwards §1 refuses. It is reported in
    the output so a reader can weigh it, and flagged in §7.
  * It does not blend R1 with R2/R3. §3 makes R1 the ordering and R2/R3 tier
    inputs, so ties on R1 break on iptm (higher better), then binder_plddt
    (higher better), then design id -- a tie-break only, never a reweighting.
  * It does not touch eligibility. A design must already be a gate survivor:
    novel seed 50, O2 PASS. Designs with no fold row cannot be ordered on R1
    and are excluded, which is reported rather than hidden.

Run:  python3 rerank_fold.py [-n 20] [--out submission_fold.csv]
Out:  submission_fold.csv   ranked list, ready for make_scfv.py --source
"""
import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FOLD = os.path.join(HERE, "results", "vm", "novel-200-of3ob-refold-metrics.csv")


def design_id(row):
    """The metrics table stores a path; pull the design id out of it."""
    m = re.search(r"(her2_cdrh3_\d+)", row.get("design") or "")
    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(HERE, "scores.csv"))
    ap.add_argument("--o2-report",
                    default=os.path.join(HERE, "o2_report_all234.csv"))
    ap.add_argument("--fold", default=FOLD)
    ap.add_argument("--run", default="novel")
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(HERE, "submission_fold.csv"))
    args = ap.parse_args()

    for p in (args.scores, args.o2_report, args.fold):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    scores = list(csv.DictReader(open(args.scores, newline="")))
    o2 = {(r["design_id"], r["cdrh3"]): r
          for r in csv.DictReader(open(args.o2_report, newline=""))}
    fold_rows = list(csv.DictReader(open(args.fold, newline="")))

    # One row per design, one model, one seed, one sample -- assert it rather
        # than assume it, because a second sample would make "the" score ambiguous.
    for field, label in (("model", "folding model"), ("seed", "folding seed"),
                         ("sample", "diffusion sample")):
        vals = {r[field] for r in fold_rows}
        if len(vals) != 1:
            sys.exit(f"{label} is not single-valued: {sorted(vals)} -- scores "
                     f"from different settings are not comparable (§8)")
    fold = {}
    for r in fold_rows:
        did = design_id(r)
        if did in fold:
            sys.exit(f"{did} has two fold rows -- cannot pick one score")
        fold[did] = r
    model = fold_rows[0]["model"]
    print(f"{len(fold)} designs folded by {model}, "
          f"seed {fold_rows[0]['seed']}, sample {fold_rows[0]['sample']}")

    # Eligibility: the same gate survivors as before. Nothing is relaxed here.
    elig = [r for r in scores
            if r["run"] == args.run
            and o2.get((r["design_id"], r["cdrh3"]), {}).get("o2_verdict") == "PASS"]
    folded = [r for r in elig if r["design_id"] in fold]
    missing = [r for r in elig if r["design_id"] not in fold]
    print(f"{len(elig)} eligible ({args.run}, O2 PASS); "
          f"{len(folded)} have fold scores, {len(missing)} do not")
    if missing:
        print("  excluded for want of a fold score, so R1 cannot order them: "
              + ", ".join(r["design_id"] for r in missing))
    if len(folded) < args.n:
        sys.exit(f"only {len(folded)} rankable, need {args.n}")

    def key(r):
        f = fold[r["design_id"]]
        return (float(f["min_ipae"]),          # R1: lower is better
                -float(f["iptm"]),             # tie-break, higher is better
                -float(f["binder_plddt"]),     # tie-break, higher is better
                r["design_id"])

    ordered = sorted(folded, key=key)
    chosen = ordered[:args.n]

    out = []
    for rank, r in enumerate(chosen, start=1):
        f = fold[r["design_id"]]
        out.append({
            "rank": rank,
            "prev_rank": r["rank"],          # its place under the composite
            "design_id": r["design_id"],
            "run": r["run"],
            "weights": r["weights"],
            "seed": r["seed"],
            "track": r["track"],
            "cdrh3": r["cdrh3"],
            "vh": r["vh"],
            "min_ipae": f["min_ipae"],
            "iptm": f["iptm"],
            "binder_plddt": f["binder_plddt"],
            "apo_holo_rmsd": f["apo_holo_rmsd"],   # reported, not ranked on
            "folding_model": f["model"],
            "folding_seed": f["seed"],
            "composite_score": r["composite_score"],  # the superseded score
            "min_mut_from_data_point":
                o2[(r["design_id"], r["cdrh3"])]["min_mut_from_data_point"],
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    # R1 must be monotonic in the output: proof the order is min_ipae's, not ours.
    ipae = [float(r["min_ipae"]) for r in out]
    if ipae != sorted(ipae):
        sys.exit("output is not ordered by min_ipae -- R1 was not applied")

    print()
    print(f"{'new':>3} {'was':>4}  {'design_id':<16} {'CDRH3':<11} "
          f"{'min_ipAE':>8} {'ipTM':>6} {'pLDDT':>7} {'apo_holo':>8}")
    for r in out:
        print(f"{r['rank']:>3} {r['prev_rank']:>4}  {r['design_id']:<16} "
              f"{r['cdrh3']:<11} {r['min_ipae']:>8} {r['iptm']:>6} "
              f"{r['binder_plddt']:>7} {r['apo_holo_rmsd']:>8}")

    # How much did the rule actually change? State it, do not bury it.
    prev_list = os.path.join(HERE, "submission_scfv_novel.csv")
    if os.path.exists(prev_list):
        prev = {r["design_id"] for r in csv.DictReader(open(prev_list, newline=""))}
        kept = [r for r in out if r["design_id"] in prev]
        print(f"\nof the previous 20 (composite-ranked), {len(kept)} are "
              f"retained, {args.n - len(kept)} replaced")
    worst = max(out, key=lambda r: float(r["apo_holo_rmsd"]))
    print(f"apo/holo RMSD across the new 20: "
          f"{min(float(r['apo_holo_rmsd']) for r in out):.1f}-"
          f"{max(float(r['apo_holo_rmsd']) for r in out):.1f} Å "
          f"(worst {worst['design_id']}); not a ranking criterion, see §7")
    print(f"\nwrote {os.path.basename(args.out)}")


if __name__ == "__main__":
    main()
