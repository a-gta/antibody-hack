# antibody-hack — HER2 / trastuzumab optimisation

Oslo AI × Antibodies Hackathon (Apheris, 14 November 2026). Target HER2,
baseline trastuzumab, **optimisation track**. Only CDRH3 is varied; the VH
framework and the VL are trastuzumab's own.

Read [`SELECTION.md`](SELECTION.md) first. It fixes the selection rule and was
committed before anything was ranked, and its §7 logs every departure from it.
*Judge the method, then the molecules.*

## The submission

**[`submission_scfv_fold.csv`](submission_scfv_fold.csv)** — 20 scFvs, ranked by
OpenFold3 interface score.

Every submitted design comes from **one run**: HER2 novel weights
`her2-novel-2026-11-09`, design seed 50. One run, one seed, one command to
reproduce it.

| | |
|---|---|
| Format | scFv, VH-(GGGGS)₃-VL, 242 aa |
| Cysteines | 4 per construct, each domain's own pair, none unpaired |
| Gates | **20/20 clean PASS** — G1 format, G2 cysteines, G3 numbering, G4 length, O1, O2, none left unchecked |
| Ranked on | `min_ipAE` from OpenFold3 (§3's R1), ties broken on ipTM then binder pLDDT |
| Folded | all 20, in the holo complex and apo; OpenFold3 seed 42, one diffusion sample |
| De novo track | not entered; nothing generated qualifies (§5) |

`abguard` is authoritative and runs on submission. We pre-check so no slot is
lost to a rejection we could have predicted.

## How the 20 were reached

| Step | Count | Where |
|---|---|---|
| Designs generated (novel, seed 50) | 200 | [`results/vm/novel-200-per-example-metrics.csv`](results/vm/novel-200-per-example-metrics.csv) |
| Distinct CDRH3s | 177 | 23 of the 200 are repeats |
| Clear the novelty pre-check | 193 | [`results/vm/novel-200-novelty.txt`](results/vm/novel-200-novelty.txt) |
| Survive scoring's inline gates and get ranked | 82 | [`scores.csv`](scores.csv) |
| Pass gate O2 against all 525,063 provided data points | 65 | [`o2_report_all234.csv`](o2_report_all234.csv) |
| Folded by OpenFold3, so R1 can order them | 62 | [`results/vm/novel-200-of3ob-refold-metrics.csv`](results/vm/novel-200-of3ob-refold-metrics.csv) |
| Submitted | 20 | [`submission_scfv_fold.csv`](submission_scfv_fold.csv) |

The three filters are not the same check and the counts are not interchangeable.
The novelty pre-check asks whether a CDRH3 is absent from the provided data *and*
≥3 mutations from trastuzumab. Gate **O2** is stricter and is the one that
rejects: it requires **strictly more than 2** point mutations from *every*
provided data point, and 17 of the 82 ranked novel designs fail it.

Ranking is **§3's own R1**: `min_ipae` from OpenFold3, ascending, with ipTM and
binder pLDDT as tie-breaks only. `score_designs.py`'s sequence-and-generator
composite ranked an earlier version of this submission and is superseded — it
was a stand-in for exactly this score while the fold outputs were thought not to
exist. They were recovered off the VM before shutdown: 142 of the 200 designs,
holo and apo. `rerank_fold.py` asserts the output is monotonic in `min_ipae`, so
the order is the criterion's.

**No predicted K_D is quoted anywhere.** Scores give rank order only, and three
cautions belong with them, all logged in §7:

- ipTM sits in a narrow 0.53–0.69 band, and §6's warning that ipTM cannot
  separate a real antibody–antigen pair from a shuffled one was never
  calibrated. `min_ipae` is the ordering because §3 named it primary, not
  because it is known to be better calibrated.
- Each score is **one diffusion sample at one seed** — a single draw, not a mean.
- **Applying R1 collapsed the diversity.** The submitted 20 hold 39 pairs within
  2 mutations of each other, against 4 in the composite-ranked list, and 16 of
  20 share the motif `ST..V[DGS]YA[FMH]D`. If that family does not bind, most of
  the 20 fail together. No diversity criterion was added, because §3–§5 define
  none and inventing one after seeing the scores is what §1 forbids. The
  previous list is kept as `submission_scfv_novel.csv` so the trade is visible.

## Reproduce the run

```bash
bash reproduce-novel-seed50.sh
```

[`reproduce-novel-seed50.sh`](reproduce-novel-seed50.sh) pins the design image by
sha256 digest rather than by tag, verifies both input files by md5 before
starting, and fingerprints the CDRH3s it produces against the original run
(`90056d11688e5f7e1279b3b4347689a0`, 200 designs). It needs the hackathon VM's
`Data/12-her2/` and one GPU; the design step took ~70 min on an H100 80GB.

It runs the design step only. Nothing in this submission is folded.

## Rebuild the submission from the run outputs

```bash
python3 rerank_fold.py -n 20 --out submission_fold.csv
```

```bash
python3 make_scfv.py --source submission_fold.csv --out submission_scfv_fold.csv
```

```bash
python3 gates.py scfv_gates_input.csv --data-points o2_data_points_fold.fasta
```

`rerank_fold.py` applies §3's R1 to the gate survivors. `rerank.py --run novel`
builds the superseded composite-ranked list, and is a filter rather than a
ranking: it drops rows and renumbers, recomputes no score, and asserts that the
survivors' previous ranks stay strictly increasing.

## What is here

| File | What it is |
|---|---|
| `SELECTION.md` | the selection rule, fixed before ranking, with the §7 deviation log |
| `rerank_fold.py` → `submission_fold.csv` | the ranking actually used: §3's R1 on the OpenFold3 scores |
| `score_designs.py` → `scores.csv` | the superseded composite, 234 rows across both runs; still the eligibility source |
| `gates.py` | the §2 hard gates, with `--self-test` |
| `check_o2.py` → `o2_report*.csv` | gate O2 against all provided data points |
| `rerank.py` | filter the ranking down to gate survivors and one run |
| `make_scfv.py` | pair each VH with trastuzumab's VL; `--verify-vl` re-derives the domain boundary with ANARCII |
| `make_fold_input.py` | prepare the submitted designs for folding |
| `results/vm/` | raw outputs of both VM runs: metrics, novelty reports, run settings |
| `check_1zvy.py` | numbering practice run on PDB 1ZVY, verifying the guide's claims |
| `submission_scfv_novel.csv` | the composite-ranked novel-only 20 — superseded, kept so the diversity trade is inspectable |
| `submission_scfv.csv` | the earlier top 20, mixed runs (16 novel + 4 exploratory) — superseded, kept for comparison |
| `results/vm/novel-200-of3ob-*` | the recovered OpenFold3 metrics and fold config that the ranking rests on |

## Known gaps

- **The folding is thin.** 142 of the 400 designs are folded, all from the novel
  seed-50 run, each as a single diffusion sample at one seed. The 200
  exploratory designs are unfolded, so they cannot be ranked by R1 at all and
  are absent from the submission for that reason as well as the run filter's.
- **ipTM is uncalibrated.** Stage 3 of the hackathon brief exists because ipTM
  often cannot separate a real antibody–antigen pair from a shuffled one. That
  calibration was never run, and these values are narrowly spread.
- **apo/holo RMSD is 3.2–9.8 Å** across the submitted 20 — the designed binders
  move substantially without the antigen. It is reported but not ranked on
  (§7), and it is the clearest unexploited signal in the data.
- **`abguard` and `immunum` were never available**, so `gates.py` reimplements
  G1–G4 and O1–O2. abguard's verdict on submission is the one that counts.
- **No thermostability or solubility prediction.** §1 commits to designs that
  attempt all three measured properties at once; the ranking's proxies for
  stability and expression are sequence features, not models.
- **`her2_cdrh3_N` ids are not unique across runs.** Join on
  `(design_id, cdrh3)`; `run`, `weights` and `seed` are carried in the
  submission for this reason (§7).
