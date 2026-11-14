#!/usr/bin/env python3
"""Rank HER2 / trastuzumab CDRH3 designs on the signals we actually measured.

We could fold only 2 of the 400 designs before the deadline, so this ranks on
sequence-derived and generator-reported signals. Every number below is computed
by code in this file. Nothing is a predicted affinity and nothing is a KD.

Run:  python3 score_designs.py
Out:  scores.csv        every eligible design, every component, with sources
      submission.csv    the ranked top 20
"""
import csv
import os
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
VM = os.path.join(HERE, "results", "vm")

RUNS = [
    ("novel", "novel-200.csv", "novel-200-novelty.txt", "her2-novel-2026-11-09", 50),
    ("exploratory", "exploratory-200.csv", "exploratory-200-novelty.txt", "her2-2026-11-07", 42),
]

# The organisers published a benchmark of these two weight sets on this exact
# task (hackathon task description, p5-6). It is a measured prior on the
# GENERATOR, not on any individual design, so it applies per run:
#
#   weights                   known binders  pass novelty  interface error
#   her2-novel-2026-11-09              4%           86%           6.8 A
#   her2-2026-11-07 exploratory      0.3%           74%           9.1 A
#
# Lower interface error means the generator places the binder more accurately.
# We fold that 2.3 A gap in as a per-run term rather than pretending the two
# runs are interchangeable.
RUN_PRIOR = {"novel": 0.10, "exploratory": 0.0}
RUN_PRIOR_WHY = "organisers' benchmark: 6.8 A vs 9.1 A interface error"

TRASTUZUMAB_CDRH3 = "WGGDGFYAMD"

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

LIABILITIES = [
    ("NG",       r"NG",        "Asn deamidation, the fastest motif"),
    ("NS/NT",    r"N[ST]",     "Asn deamidation, slower"),
    ("N-x-S/T",  r"N[^P][ST]", "N-linked glycosylation sequon"),
    ("DG",       r"DG",        "Asp isomerisation"),
    ("DP",       r"DP",        "acid-labile, fragmentation"),
    ("DS/DT",    r"D[ST]",     "Asp isomerisation, slower"),
    ("M",        r"M",         "Met oxidation"),
    ("unpaired C", r"C",       "free thiol, disulfide scrambling"),
]

MIN_MUT_MARGIN = 5


def read_novelty(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        p = line.split()
        if len(p) >= 4 and p[0].startswith("her2_cdrh3_"):
            out[p[0]] = (p[3] == "yes")
    return out


def num(row, col):
    v = (row.get(col) or "").strip()
    if not v or v.lower() in ("none", "nan", "na"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def liabilities(cdrh3):
    hits = []
    for name, pattern, _why in LIABILITIES:
        n = len(re.findall("(?=" + pattern + ")", cdrh3))
        if n:
            hits.extend([name] * n)
    return hits


def hydropathy(seq):
    return sum(KD.get(a, 0.0) for a in seq) / max(len(seq), 1)


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def aromatic_fraction(seq):
    return sum(seq.count(a) for a in "WFY") / max(len(seq), 1)


def score(row):
    cdr = row["cdrh3"]
    parts = {}
    hits = liabilities(cdr)
    parts["liability"] = -0.12 * len(hits)
    h = hydropathy(cdr)
    parts["hydrophobicity"] = -0.20 * max(0.0, h - hydropathy(TRASTUZUMAB_CDRH3))
    arom = aromatic_fraction(cdr)
    parts["aromatic"] = -0.60 * max(0.0, arom - 0.40)
    parts["charge"] = -0.08 * abs(net_charge(cdr))
    parts["mutation_distance"] = -0.05 * max(0.0, row["_mut"] - 7.0)
    parts["generator_prior"] = RUN_PRIOR.get(row["_run"], 0.0)
    # Round to the precision we report. Unrounded, two designs that both print
    # as 0.8280 differ at ~1e-16 (0.82800000000000007 vs 0.82799999999999996),
    # so the first sort key never ties and the stated tie-break below never
    # fires -- order inside a tie group would be float noise, not a rule.
    return round(1.0 + sum(parts.values()), 4), parts, hits


def main():
    rows = []
    for label, csv_name, nov_name, weights, seed in RUNS:
        path = os.path.join(VM, csv_name)
        if not os.path.exists(path):
            print("  skip " + csv_name + ", not found")
            continue
        nov = read_novelty(os.path.join(VM, nov_name))
        for r in csv.DictReader(open(path)):
            r["_run"] = label
            r["_weights"] = weights
            r["_seed"] = seed
            r["_mut"] = num(r, "n_mut_vs_trastuzumab") or 0.0
            r["_hot"] = num(r, "hotspot_obedience")
            r["_nov"] = nov.get(r["design"], None)
            rows.append(r)

    print("loaded %d designs from %d runs\n" % (len(rows), len(RUNS)))

    gated, reasons = [], Counter()
    for r in rows:
        why = []
        vh = r.get("vh", "")
        if not vh or len(vh) > 250:
            why.append("length>250")
        if vh.count("C") % 2 != 0:
            why.append("unpaired_cys")
        if r["_nov"] is not True:
            why.append("novelty_fail")
        if r["_hot"] != 1.0:
            why.append("hotspot<1.0")
        if r["_mut"] < MIN_MUT_MARGIN:
            why.append("mut<%d" % MIN_MUT_MARGIN)
        r["_gate_failed"] = ";".join(why)
        if why:
            reasons[why[0]] += 1
        else:
            gated.append(r)

    print("gate results")
    print("  %4d eligible" % len(gated))
    for k, v in reasons.most_common():
        print("  %4d rejected: %s" % (v, k))

    seen, uniq = set(), []
    for r in sorted(gated, key=lambda x: x["design"]):
        if r["cdrh3"] in seen:
            continue
        seen.add(r["cdrh3"])
        uniq.append(r)
    print("  %4d distinct CDRH3 loops\n" % len(uniq))

    for r in uniq:
        r["_score"], r["_parts"], r["_hits"] = score(r)
    # Ties are common because the penalty model has a ceiling at 1.0 + prior.
    # Break them on stated grounds, never on file order: prefer the better
    # generator, then 6 mutations (mid-range: novelty margin without drifting
    # from a validated binder), then design id for reproducibility.
    uniq.sort(key=lambda r: (-r["_score"],
                             -RUN_PRIOR.get(r["_run"], 0.0),
                             abs(r["_mut"] - 6.0),
                             r["design"]))

    cols = ["rank", "design_id", "run", "weights", "seed", "track", "cdrh3",
            "cdrh3_len", "n_mut_vs_trastuzumab", "hotspot_obedience",
            "liability_hits", "n_liabilities", "hydropathy", "net_charge",
            "aromatic_fraction", "iptm", "min_ipae", "binder_plddt",
            "folding_model", "composite_score", "score_source_tool",
            "score_source_version", "vh"]
    with open(os.path.join(HERE, "scores.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for i, r in enumerate(uniq, 1):
            w.writerow({
                "rank": i, "design_id": r["design"], "run": r["_run"],
                "weights": r["_weights"], "seed": r["_seed"], "track": "optimisation",
                "cdrh3": r["cdrh3"], "cdrh3_len": r.get("cdrh3_len", ""),
                "n_mut_vs_trastuzumab": int(r["_mut"]),
                "hotspot_obedience": r["_hot"],
                "liability_hits": "|".join(r["_hits"]) or "none",
                "n_liabilities": len(r["_hits"]),
                "hydropathy": "%.3f" % hydropathy(r["cdrh3"]),
                "net_charge": net_charge(r["cdrh3"]),
                "aromatic_fraction": "%.3f" % aromatic_fraction(r["cdrh3"]),
                "iptm": r.get("iptm", ""), "min_ipae": r.get("min_ipae", ""),
                "binder_plddt": r.get("binder_plddt", ""),
                "folding_model": r.get("folding_model", ""),
                "composite_score": "%.4f" % r["_score"],
                "score_source_tool": "score_designs.py",
                "score_source_version": "git HEAD",
                "vh": r.get("vh", ""),
            })

    top = uniq[:20]
    with open(os.path.join(HERE, "submission.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "design_id", "track", "cdrh3", "vh"])
        for i, r in enumerate(top, 1):
            w.writerow([i, r["design"], "optimisation", r["cdrh3"], r.get("vh", "")])

    print("TOP 20")
    print("%3s %18s %12s %12s %3s %4s %6s" % ("#", "design", "run", "cdrh3", "mut", "liab", "score"))
    print("-" * 68)
    for i, r in enumerate(top, 1):
        print("%3d %18s %12s %12s %3d %4d %6.3f" % (
            i, r["design"], r["_run"], r["cdrh3"], int(r["_mut"]), len(r["_hits"]), r["_score"]))

    print("\ntop-20 split by run: %s" % dict(Counter(r["_run"] for r in top)))
    print("top-4 (tested regardless): %s" % [r["design"] for r in top[:4]])
    print("\nwritten: scores.csv, submission.csv")


if __name__ == "__main__":
    main()
