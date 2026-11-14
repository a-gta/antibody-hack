# antibody-hack — HER2 / trastuzumab optimisation

Oslo AI × Antibodies Hackathon (Apheris, 14 November 2026). Target HER2,
baseline trastuzumab, **optimisation track**. Only CDRH3 is varied; the VH
framework and the VL are trastuzumab's own.

Read [`SELECTION.md`](SELECTION.md) first. It fixes the selection rule and was
committed before anything was ranked, and its §7 logs every departure from it.
*Judge the method, then the molecules.*

## The submission

**[`submission_scfv_novel.csv`](submission_scfv_novel.csv)** — 20 scFvs, ranked.

Every submitted design comes from **one run**: HER2 novel weights
`her2-novel-2026-11-09`, design seed 50. One run, one seed, one command to
reproduce it.

| | |
|---|---|
| Format | scFv, VH-(GGGGS)₃-VL, 242 aa |
| Cysteines | 4 per construct, each domain's own pair, none unpaired |
| Gates | **20/20 clean PASS** — G1 format, G2 cysteines, G3 numbering, G4 length, O1, O2, none left unchecked |
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
| Submitted | 20 | [`submission_scfv_novel.csv`](submission_scfv_novel.csv) |

The three filters are not the same check and the counts are not interchangeable.
The novelty pre-check asks whether a CDRH3 is absent from the provided data *and*
≥3 mutations from trastuzumab. Gate **O2** is stricter and is the one that
rejects: it requires **strictly more than 2** point mutations from *every*
provided data point, and 17 of the 82 ranked novel designs fail it.

Ranking is `score_designs.py`, a sequence-and-generator model, not a folding
score: only 2 of 400 designs could be folded before the deadline, so `iptm`,
`min_ipae` and `binder_plddt` are empty for every ranked row. That is a
departure from §3's stated ranking basis and is logged in §7, first entry.
**No predicted K_D is quoted anywhere.** Scores give rank order only.

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
python3 rerank.py --run novel -n 20 --out submission_o2_novel.csv
```

```bash
python3 make_scfv.py --source submission_o2_novel.csv --out submission_scfv_novel.csv
```

```bash
python3 gates.py scfv_gates_input.csv --data-points o2_data_points_novel.fasta
```

`rerank.py` is a filter, not a ranking: it drops rows and renumbers, recomputes
no score, and asserts that the survivors' previous ranks stay strictly
increasing, so if A outranked B it still does.

## What is here

| File | What it is |
|---|---|
| `SELECTION.md` | the selection rule, fixed before ranking, with the §7 deviation log |
| `score_designs.py` → `scores.csv` | the ranking actually used, 234 rows across both runs |
| `gates.py` | the §2 hard gates, with `--self-test` |
| `check_o2.py` → `o2_report*.csv` | gate O2 against all provided data points |
| `rerank.py` | filter the ranking down to gate survivors and one run |
| `make_scfv.py` | pair each VH with trastuzumab's VL; `--verify-vl` re-derives the domain boundary with ANARCII |
| `make_fold_input.py` | prepare the submitted designs for folding |
| `results/vm/` | raw outputs of both VM runs: metrics, novelty reports, run settings |
| `check_1zvy.py` | numbering practice run on PDB 1ZVY, verifying the guide's claims |
| `submission_scfv.csv` | the earlier top 20, mixed runs (16 novel + 4 exploratory) — superseded, kept for comparison |

## Known gaps

- **Nothing is folded.** 2 of 400 designs have structures. ipTM, min_ipAE and
  binder_pLDDT play no part in the ranking.
- **`abguard` and `immunum` were never available**, so `gates.py` reimplements
  G1–G4 and O1–O2. abguard's verdict on submission is the one that counts.
- **No thermostability or solubility prediction.** §1 commits to designs that
  attempt all three measured properties at once; the ranking's proxies for
  stability and expression are sequence features, not models.
- **`her2_cdrh3_N` ids are not unique across runs.** Join on
  `(design_id, cdrh3)`; `run`, `weights` and `seed` are carried in the
  submission for this reason (§7).
