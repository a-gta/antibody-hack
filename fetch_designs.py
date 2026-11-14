"""Pull the designed sequences from the two 1ZVY RFD4-Proteina jobs.

Reads the key from the environment; it is never written to a file or passed
on the command line (argv is visible to other processes via ps).

    export JAPANFOLD_API_KEY=...        # in your own terminal
    python fetch_designs.py             # or pass job ids as arguments
"""

import json
import os
import sys
import urllib.error
import urllib.request

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY = os.environ.get("JAPANFOLD_API_KEY")

JOBS = sys.argv[1:] or [
    "design-72c173",  # rfd4-binder    B1-129, /0, 110-130
    "design-17e757",  # rfd4-scaffold  B1-129, /0, A2-96, 20, A117-125
]

# What check_1zvy.py established, repeated here so the output is checkable.
KEPT = {"design-72c173": 129, "design-17e757": 129 + 95 + 9}


def resolve(job):
    """The workbench shows a job's `name` (design-72c173); the API wants its
    32-hex `id`. Accept either."""
    if len(job) == 32 and all(c in "0123456789abcdef" for c in job):
        return job, job
    listing = get("/v1/jobs?limit=50")
    items = listing.get("jobs") or listing.get("items") or listing.get("data") or []
    for j in items:
        if j.get("name") == job or j.get("id") == job:
            return j["id"], f"{j.get('name')} [{j.get('protocol')}]"
    raise SystemExit(f"no job named {job!r} on this account")


def get(path):
    req = urllib.request.Request(
        f"{BASE}{path}", headers={"Authorization": f"Bearer {KEY}"}
    )
    try:
        with urllib.request.urlopen(req) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")[:400]
        raise SystemExit(f"{path} -> HTTP {e.code}\n{body}")


def main():
    if not KEY:
        raise SystemExit(
            "JAPANFOLD_API_KEY is not set. Export it in your terminal first "
            "(see japanfold.aiand.com/account/token), then re-run."
        )

    for job in JOBS:
        job_id, label = resolve(job)
        res = get(f"/v1/jobs/{job_id}/results")
        designs = res.get("designs") or res.get("results") or []
        print(f"\n=== {label} === {len(designs)} designs")

        for d in sorted(designs, key=lambda x: x.get("final_rank", 99)):
            seq = d.get("designed_sequence") or ""
            rank = d.get("final_rank")
            rmsd = d.get("motif_rmsd")
            chains = d.get("num_chains")
            expected = d.get("chains_expected")

            flags = []
            if expected is not None and chains != expected:
                flags.append(f"CHAINS {chains} != expected {expected}")
            kept = KEPT.get(job)
            if kept is not None and seq:
                # the designed stretch is whatever is not kept framework/target
                flags.append(f"designed ~{len(seq) - kept} aa" if len(seq) > kept
                             else f"length {len(seq)} <= kept {kept}?")

            print(f"\n  #{rank}  motif_rmsd={rmsd}  chains={chains}"
                  + (f"  [{'; '.join(flags)}]" if flags else ""))
            print(f"  {d.get('id')}  {d.get('structure')}")
            if seq:
                for i in range(0, len(seq), 60):
                    print(f"    {i + 1:>4} {seq[i:i + 60]}")
            else:
                print("    (no designed_sequence in the result)")

    print("\nNeither model emits a confidence score. To get one, refold a "
          "sequence as a prediction job and read its pLDDT / iptm.")


if __name__ == "__main__":
    main()
