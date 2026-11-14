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

Track split across the 20: `[N]` optimisation, `[N]` de novo. The prize pool
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

---

## 8. Provenance

Scores live in `scores.csv`, each with tool and version. Numbering is IMGT via
ANARCII 2.0.8; IMGT positions are mapped explicitly to PDB residue numbers and
never used interchangeably. Weights, seeds and settings are recorded per job.

**Verified, not assumed** (practice run on 1ZVY, `check_1zvy.py`): the
participant guide states residue A1 is unresolved. Chain A is observed 1–125
with no gaps and **A1 is resolved**. Chain B is observed 1–129. IMGT CDR3
105–117 maps to A97–A116, flanked by Cys A96 and Trp A117.

---

*Committed before ranking. Judge the method, then the molecules.*
