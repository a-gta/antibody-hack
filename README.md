# antibody-hack — HER2 / trastuzumab optimisation

**Team: Totally Claude Spies & the Boys.** Oslo AI × Antibodies Hackathon
(Apheris, 14 November 2026). Target HER2, baseline trastuzumab, **optimisation
track**. Only CDRH3 is varied; the VH framework and the VL are trastuzumab's own.

Read [`SELECTION.md`](SELECTION.md) first. It fixes the selection rule and was
committed before anything was ranked, and its §7 logs every departure from it.
*Judge the method, then the molecules.*

## The submission

**[`submission_scfv_diverse.csv`](submission_scfv_diverse.csv)** — 20 scFvs,
ranked by OpenFold3 interface score, with a minimum pairwise distance between
them.

Every submitted design comes from **one run**: HER2 novel weights
`her2-novel-2026-11-09`, design seed 50. One run, one seed, one command to
reproduce it.

| | |
|---|---|
| Format | scFv, VH-(GGGGS)₃-VL, 242 aa |
| Cysteines | 4 per construct, each domain's own pair, none unpaired |
| Gates | **20/20 clean PASS** — G1 format, G2 cysteines, G3 numbering, G4 length, O1, O2, none left unchecked |
| Ranked on | `min_ipAE` from OpenFold3 (§3's R1), ties broken on ipTM then binder pLDDT |
| Diversity | every pair of CDRH3s differs at **≥3 positions**; closest pair 3, mean 5.87 |
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
| Admitted by the D = 3 diversity rule, walking R1 from the top | 20 of the first 31 | [`submission_scfv_diverse.csv`](submission_scfv_diverse.csv) |

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
- **R1 alone collapsed the diversity, so a diversity rule was added.** Ranking
  on `min_ipAE` by itself put 16 of 20 designs in one motif family,
  `ST..V[DGS]YA[FMH]D`, with 39 pairs within 2 mutations of each other — a
  submission staking most of its 20 slots on one family binding. §3–§5 define
  no diversity criterion, so **D = 3** was pre-registered as a new one (§7, last
  entry) rather than slipped in: the threshold is the challenge's own definition
  of a distinct sequence (O2 needs >2 mutations from any data point), no other
  value was tried, and `rerank_diverse.py` fails rather than relaxing D. It
  costs R1 — reaching rank 31 to fill 20, `min_ipAE` 0.4553–0.5988 instead of
  0.4553–0.5705 — and buys ~20 independent attempts instead of one family's.
  `submission_scfv_fold.csv` keeps the pure-R1 list for comparison.

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
python3 rerank_fold.py -n 62 --out pool_fold_ranked.csv
```

```bash
python3 rerank_diverse.py -n 20 --out submission_diverse.csv
```

```bash
python3 make_scfv.py --source submission_diverse.csv --out submission_scfv_diverse.csv
```

```bash
python3 gates.py scfv_gates_input.csv --data-points o2_data_points_diverse.fasta
```

`rerank_fold.py` applies §3's R1 to the gate survivors; `rerank_diverse.py`
walks that order and admits a design only if it is ≥3 mutations from every design
already admitted, so it filters without reordering — both assert their output is
monotonic in `min_ipae`. `rerank.py --run novel`
builds the superseded composite-ranked list, and is a filter rather than a
ranking: it drops rows and renumbers, recomputes no score, and asserts that the
survivors' previous ranks stay strictly increasing.

## What is here

| File | What it is |
|---|---|
| `SELECTION.md` | the selection rule, fixed before ranking, with the §7 deviation log |
| `rerank_fold.py` → `submission_fold.csv` | §3's R1 on the OpenFold3 scores, pure |
| `rerank_diverse.py` → `submission_diverse.csv` | the submission actually used: R1 plus the D = 3 diversity rule |
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
