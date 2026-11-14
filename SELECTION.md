# Selection rule

**Team:** `[NAME]`
**Target:** HER2 · **Baseline antibody:** trastuzumab
**Tracks entered:** `[optimisation | de novo | both]`
**Written at:** `[TIME]` — before any candidate was ranked.

This file fixes how we rank our candidates and which 20 we submit. It is
committed before we have ranked anything, so our selection is a method and not
a story told afterwards. Every departure is logged in §7.

---

## 1. What we are optimising for

Adaptyv measures three properties against trastuzumab's own measured baseline:
binding to HER2 (SPR), epitope preservation (competitive binding), and
thermostability. Tiers count how many a design beats — S-tier 3, A-tier 2,
B-tier 1. Teams rank by S-tier count first.

The brief states: *"Ten designs that each improve one property lose to two
designs that improve all three."*

**Therefore we do not diversify across properties.** Every design we submit is
a credible attempt at all three simultaneously. A candidate predicted to win on
affinity alone is ranked below one predicted to win narrowly on all three.

> `[One sentence: who needs an improved trastuzumab and what they would do with it.]`

---

## 2. Hard gates — abguard decides, we pre-check

`abguard` is authoritative and runs on submission. We run it ourselves first so
no slot is lost to a rejection we could have predicted.

**Shared, every design:**

| # | Gate | Threshold |
|---|------|-----------|
| G1 | Format | VHH, scFv or Fab |
| G2 | Unpaired cysteines | none |
| G3 | Recognisable as an antibody | `immunum` numbering confidence ≥ 0.1 |
| G4 | Total length | ≤ 250 aa (single chain, or VH+VL combined) |

**De novo track only:**

| # | Gate | Threshold |
|---|------|-----------|
| D1 | Identity to nearest known antibody | < 70% |
| D2 | CDRH3 mutations from trastuzumab | ≥ 2 |

**Optimisation track only:**

| # | Gate | Threshold |
|---|------|-----------|
| O1 | CDRH3 identity to trastuzumab | < 70% |
| O2 | Point mutations from any provided data point | strictly > 2 |

**Track is assigned before generation, never after.** A design that misses its
own track's criteria is rejected outright, not reclassified.

Counts: generated `[N]` → passed abguard `[N]`.

---

## 3. Ranking — gate survivors only

Scores give **rank order only**. We do not quote a predicted K_D. For the HER2
exercise RFD4's own pass column stays empty, so we rank on the three scores the
brief names.

| Rank | Score | Source | Use |
|------|-------|--------|-----|
| R1 | `min_ipae` | Apheris `refold-score` | primary ordering |
| R2 | `iptm` | Apheris `refold-score` or `openfold3` | tiers, thresholded per §4 |
| R3 | `binder_plddt` | Apheris `refold-score` | tiers |
| R4 | Binding ΔΔG | Nebius `stab-ddg` | **optimisation track: the affinity axis directly** |
| R5 | Thermostability proxy | ESM-IF1 / AntiFold log-likelihood, germline closeness | **promotes A-tier to S-tier — not a tie-break** |
| R6 | Expression proxy | Aggrescan3D / SoluProt / NetSolP | a design that cannot be expressed is never measured |

**Boltz2 co-folds give no ipTM.** They return confidence and pTM only. We use
Boltz2 for breadth and interface contact counts, never as a source of the ipTM
that §4 thresholds. Calibrating on one score and ranking on another would make
the threshold meaningless.

Model choice: **`her2-novel-2026-11-09`**, because the brief's own benchmark
puts it at 86% passing the novelty filter with 6.8 Å interface error. We
rejected `her2-2026-11-07` conservative and balanced: at 0% and 6% novelty pass
they generate designs abguard discards.

---

## 4. The ipTM threshold carries its false-positive rate

The organisers supply **454 VHH–antigen pairs measured as non-binders, 21 of
which score AlphaFold3 ipTM ≥ 0.7.** Whatever threshold we use, we quote the
false-positive rate that set implies for it, rather than passing 0.7 off as a
bare convention with no error attached.

Threshold used: **ipTM ≥ 0.70** → implies **4.6%** false positives on the
measured negatives (21/454).

**Status: inherited, not yet recalibrated — stated plainly rather than dressed
up.** 21/454 is the only point on the curve the organisers publish, and it is
the 0.70 point. Moving the threshold requires the per-pair ipTM scores for all
454 negatives; we do not have that file at the time of writing. So 0.70 stands,
but with its false-positive rate measured and quoted, which is the half of the
claim we can actually check. We are not pretending to have calibrated what we
have not.

On obtaining the per-pair scores: recompute the FPR at 0.70 / 0.75 / 0.80 /
0.85, take the knee, and log the change in §7 with the new number and the time.
**No threshold change once ranking has begun** — at that point a moved
threshold is exactly the "story told afterwards" this document opens by
refusing.

---

## 5. The ranked list of 20

Every team's **top 4 are tested regardless of placement**; positions 5–20 only
matter if we place. So rank order, not list length, is the deliverable.

| Positions | What goes here |
|-----------|---------------|
| 1–4 | Highest predicted probability of S-tier. These will be tested |
| 5–12 | Next best, same criteria |
| 13–16 | Upside if we win the pitch |
| 17–20 | Spares against duplication between teams |

Track split across the 20: **20** optimisation, **0** de novo — see §7
(2026-11-14 17:20); every design we generated varies only CDRH3 in
trastuzumab's framework, so none is de novo and half the prize pool is not
contested. The prize pool
divides equally between tracks, so entering both doubles our chances of placing.

---

## 6. What would make us abandon this rule

- A gate kills every candidate → relax that gate **once**, logged, new threshold
  stated before re-running.
- Apheris OpenFold unavailable → R1 falls back to `[TOOL]`, stated in the pitch.

Anything else: the rule stands.

---

## 7. Deviation log

| Time | What changed | Why | Who approved |
|------|-------------|-----|--------------|
| 2026-11-14 16:24 | **§3 ranking basis replaced.** R1 `min_ipae`, R2 `iptm`, R3 `binder_plddt` and R4 binding ΔΔG were not used. Ranking is instead a sequence-and-generator model in `score_designs.py`: liability motif count, hydropathy relative to trastuzumab's CDRH3, aromatic fraction, net charge, distance from 6 CDRH3 mutations, plus a per-run generator prior (+0.10 for `her2-novel-2026-11-09`) taken from the organisers' published 6.8 Å vs 9.1 Å interface error. | Only 2 of 400 designs could be folded before the deadline. `iptm`, `min_ipae`, `binder_plddt` and `folding_model` are empty for all 234 ranked rows in `scores.csv` — R1–R4 had no inputs to order on. §6 anticipated an OpenFold fallback for R1; it did not anticipate losing R1–R4 together, so this is a departure from the rule, not an application of it. | `[ ]` |
| 2026-11-14 16:24 | **Gate set differs from §2.** `score_designs.py` applies its own inline gates (total length ≤ 250, even cysteine count, novelty pass, `hotspot_obedience == 1.0`, ≥ 5 CDRH3 mutations) and does not call `gates.py`. Relative to §2 it omits G3 (numbering confidence) and O1 (CDRH3 identity < 70%); it tightens the mutation floor from "strictly > 2" to ≥ 5; it adds `hotspot_obedience == 1.0`; and it applies the novelty filter to the optimisation track, where §2 scopes novelty to de novo. | The 241-eligible / 234-distinct counts in `scores.csv` are therefore **not** the §2 abguard pre-check. `gates.py` (committed 11:33, before ranking) remains the §2 implementation and still reports D1 and O2 as uncheckable without a known-antibody database and the provided data points. abguard stays authoritative. | `[ ]` |
| 2026-11-14 16:24 | **Tie-break added.** §3–§5 specify no tie-break. Composite scores are rounded to the reported 4 decimal places, then ties broken on generator prior, then `abs(mutations − 6)`, then design id. | The penalty model ceilings at `1.0 + prior`, so ties are common: 14 of the top 20 share a score of 1.0000. Unrounded, designs printing as the same value differed at ~1e-16, so the first sort key never tied and order inside a tie group was floating-point noise rather than a stated rule. | `[ ]` |
| 2026-11-14 16:45 | **Submission format corrected: VH-only → scFv.** `submission.csv` carried a bare 120-aa VH and no light chain. `make_scfv.py` pairs every design with trastuzumab's own VL as VH-(GGGGS)₃-VL, 242 aa, written to `submission_scfv.csv`. Rank order is unchanged; only CDRH3 was ever varied and both the VH framework and the VL are trastuzumab's own, so nothing in §3–§5 is re-ranked. | **§2's G1 did not catch this and would not have.** `gates.py:gate_g1_format` reads one H domain with no L as `VHH` and passes it; run on the VH-only submission it reported 20/20 `VHH` 120 aa PASS\*. The brief's own type-aware comparison (p. 3, "a VHH design is compared with known VHHs **or heavy chains**") suggests abguard would likely have passed it too. The failure was never the gate: trastuzumab's HER2 epitope is contacted by both VH and VL CDRs, so a VH-only construct forfeits **epitope preservation** — one of the three properties Adaptyv measures — and contradicts §1's commitment that every design attempt all three at once. Re-run on the scFvs: 20/20 `scFv` 242 aa PASS\*, same unchecked O2. | `[ ]` |
| 2026-11-14 16:45 | **O2 flagged as a live rejection risk, not yet checked.** Unchanged from §2, but recorded: O2 requires *strictly more than 2* point mutations from **any provided data point**, and the organisers provide 524,346 trastuzumab CDRH3 variants of up to 10 mutated positions (brief p. 4). Our designs differ from trastuzumab only inside a 10-aa CDRH3, at ~6 positions. | `gates.py` reports O2 `unchecked` for all 20 because that Zenodo variant set was never pulled onto the VM. This is checkable offline and is the one remaining gate that can still cost real submission slots. It is logged here rather than left as a silent `unchecked`. | `[ ]` |
| 2026-11-14 17:05 | **O2 checked. 13 of the 20 submitted designs fail it.** `check_o2.py` runs O2 against every provided trastuzumab data point — the 524,346 CDRH3 variants plus the 753 distinct SPR-submitted CDRH3s, 525,062 distinct points in all — from the brief's own citation (Chinery et al., DOI 10.5281/zenodo.22310365). The split sizes sum to 367,042 + 78,652 + 78,652 = 524,346, matching the brief exactly, which is how we know it is the set O2 means. Verdict: **7/20 pass, 13/20 fail**, including **rank 1** (2 mutations from `WDEGVFYAFD`) and four designs at a distance of **1**. §5's top 4 are tested regardless of placement, so one of the four guaranteed slots was going to be rejected outright. | The designs vary nothing but a 10-aa CDRH3 in trastuzumab's own framework, and every data point is a CDRH3 in that same framework, so Hamming distance on the 10-mer **is** the point-mutation count — exact, not a proxy. Indel data points (20 nine-mers) go through `gates.py`'s own alignment-based `mutation_count` rather than a second implementation. Reproduced independently: `gates.py --data-points o2_data_points.fasta` returns the same 13 failures with the same distances. 118 of the 234 ranked designs in `scores.csv` do pass O2, so a full 20 can be backfilled from the existing ranking without generating anything new. **That re-ranking is not done here** — it changes §5's list and is a decision, not a correction. | `[ ]` |
| 2026-11-14 17:05 | **Bug fixed in `gates.py` gate O2: it silently passed every scFv.** `gate_o2_point_mutations` compared `max(chains, key=len)` against each data point. For a VH-only candidate that is the VH, and the gate was right. For an scFv it is the whole 242-aa construct, so measuring it against a 120-aa VH data point counted the linker and the entire VL as ~137 extra "mutations" and no scFv could ever fail. It now compares the candidate's **heavy V domain** against each data point's heavy V domain, and data points are reduced to their heavy domain on load. | Caught because the two checks disagreed: `check_o2.py` said 13 failures, `gates.py` on the same scFvs said 20/20 pass. Running `gates.py` on the VH-only form of the same designs reproduced `check_o2.py` exactly, which located the fault in the gate rather than in the data. **Our own §2 fix for the VL finding had introduced a false PASS on the gate most likely to reject us** — the scFv correction and this bug were logged in the same hour, and without the O2 check the submission would have looked cleaner than it was. `--self-test` still passes. | `[ ]` |
| 2026-11-14 17:05 | **`design_id` does not identify a molecule; `run` added to the submission.** `scores.csv` holds 234 rows under only 168 distinct `design_id`s: the `novel` and `exploratory` runs reuse ids for different CDRH3s (`her2_cdrh3_98` is `WGGGRQYTYD` in one, `TTAAVGYAHD` in the other). 7 of the 20 submitted rows carry such an ambiguous id. `submission_scfv.csv` now carries `run`, `weights` and `seed`, joined from `scores.csv` on `(design_id, cdrh3)`. | All 7 resolve uniquely by CDRH3, so **no wrong molecule was submitted** — the sequences in `submission.csv` are correct. But the identifier alone did not name a molecule, and the original `submission.csv` dropped `run`, so a submitted row could not be traced to the job that produced it, against §8's "weights, seeds and settings are recorded per job". The generator prior that fed the ranking (+0.10 for `her2-novel-2026-11-09`) is a per-run quantity, so the run is part of the ranking's own provenance, not a label. | `[ ]` |
| 2026-11-14 17:20 | **Top 20 rebuilt from the O2 survivors. 13 of 20 designs replaced.** `rerank.py` drops the 116 of 234 ranked designs that fail O2 and takes the first 20 that remain, renumbered 1–20 (`submission_o2.csv`), which `make_scfv.py --source` pairs into scFvs. 7 of the previous 20 are retained, 13 replaced by designs previously ranked 21–38. New rank 1 is `her2_cdrh3_150` (was 2). **Verification on the new 20: `check_o2.py` 20/20 pass against all 525,062 data points; `gates.py --data-points` 20/20 clean `PASS` with no gate left unchecked** — G1 `scFv`, G2, G3, G4, O1 and O2 all actually checked, where the original submission had O2 unchecked and 13 latent failures. | **This is a filter, not a re-ranking, and the distinction is the point.** §4 forbids moving a threshold once ranking has begun and §1 refuses stories told afterwards, so no score was recomputed and no weight touched: the order is `scores.csv`'s existing rank order with the §7 tie-break already applied, and `rerank.py` asserts that the survivors' previous ranks stay strictly increasing, so if A outranked B it still does. §3 already scopes ranking to "gate survivors only" — the O2 failures were never eligible, they were ranked while O2 was unchecked. Removing them applies the pre-registered rule rather than departing from it. What genuinely changed is which molecules occupy the 20 slots, and §5's top 4 — the four tested regardless of placement — are now four designs that can actually be submitted. | `[ ]` |
| 2026-11-14 17:20 | **§5 track split filled: 20 optimisation, 0 de novo.** Recorded, not changed. | Every design generated varies only CDRH3 within trastuzumab's framework, so nothing qualifies for the de novo track (D1 needs <70% identity to any known antibody). §5 notes the prize pool divides equally between tracks, so entering one track contests half of it. This was already true of the original submission; it is logged here because the file left it as `[N]`. | `[ ]` |

---

## 8. Provenance

Scores live in `scores.csv`, each with tool and version. Numbering is IMGT via
ANARCII 2.0.8; IMGT positions are mapped explicitly to PDB residue numbers and
never used interchangeably. Weights, seeds and settings are recorded per job.

**Verified, not assumed** (practice run on 1ZVY, `check_1zvy.py`): the
participant guide states residue A1 is unresolved. Chain A is observed 1–125
with no gaps and **A1 is resolved**. Chain B is observed 1–129. IMGT CDR3
105–117 maps to A97–A116, flanked by Cys A96 and Trp A117.

**Verified, not assumed** (light chain, `make_scfv.py --verify-vl`): the VL
paired into every scFv is PDB 1N8Z chain A, "Herceptin Fab — light chain", as
RCSB returns it. The variable-domain boundary is set by ANARCII, not by hand —
numbering the full chain returns `chain_type K`, confidence 30.7026, IMGT 1–127
over the first 107 residues, and the kappa constant domain that follows is
unnumbered and dropped. `--verify-vl` re-derives this and refuses to write if
ANARCII disagrees. Cysteines: 2 in VH, 2 in VL, 4 per scFv, each domain's own
pair.

**Checked, not inherited** (gate O2, `check_o2.py`): the provided data points
are the 524,346 CDRH3 variants and 753 distinct SPR-submitted CDRH3s of
Chinery et al., DOI 10.5281/zenodo.22310365 (BSD-3), fetched from Zenodo
because the hackathon S3 bucket needs the VM credentials. Split sizes sum to
the brief's stated 524,346 exactly. Per-design distances are in
`o2_report.csv`, all 234 ranked designs in `o2_report_all234.csv`, and
`o2_data_points.fasta` holds the nearest data point per design so
`gates.py --data-points` reproduces the verdict without re-reading 525k rows.

---

*Committed before ranking. Judge the method, then the molecules.*
