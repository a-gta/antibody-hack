#!/usr/bin/env bash
# Reproduce: HER2 novel, 200 designs, seed 50.
#
# One command, no settings to choose. Everything that can change an output is
# pinned here: the image by digest, both seeds, and the two input files by
# checksum. Run it from a folder holding Data/12-her2/.
#
#     bash reproduce-novel-seed50.sh
#
# It refuses to start if an input file differs from the one used on
# 14 November 2026, and at the end it compares the CDRH3s it produced against
# the fingerprint of the original run.
#
# Original run: 2026-11-14T13:43:00+00:00, NVIDIA H100 80GB HBM3,
# driver 580.173.02, Docker 29.9.2. ~70 minutes for the design step.
set -euo pipefail

# ---- pinned ------------------------------------------------------------------
# The tag apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 resolved to this
# digest. Pulling by digest means a retagged image cannot silently change the run.
RFD_DIGEST=quay.io/apheris/oslo-hackathon@sha256:11211fe9678a8365cca1e4ef77eae26af291d9fd0226ad2bac030f33b27f1024

N=200                 # designs
SEED=50               # design seed, also recorded per design in the output CSV
OUT=runs/her2-n200-novel-seed50
SPEC=Data/12-her2/her2_cdrh3.yaml

# HER2 novel takes no --sequence-tsr-k: the weights' own default stands. Passing
# it would make this balanced (2) or exploratory (1) instead, on other weights.

# md5 of the two inputs, as used in the original run.
INPUT_MD5="57c41dfc18e914224f75733252f4d5a5  Data/12-her2/her2_fv_1n8z.cif
36cf2373e9cb6233c1da8b3177132312  Data/12-her2/her2_cdrh3.yaml"

# md5 of the 200 CDRH3s, sorted, newline-joined. The scientific invariant.
EXPECTED_CDRH3_MD5=90056d11688e5f7e1279b3b4347689a0

# ---- check the inputs --------------------------------------------------------
echo "$INPUT_MD5" | md5sum -c - || {
  echo "Inputs differ from the original run. Stopping: the result would not be comparable." >&2
  exit 1
}

# ---- run ---------------------------------------------------------------------
mkdir -p "$OUT"
rm -rf "$OUT/design"

docker run --rm --gpus all --shm-size 16g \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work \
  "$RFD_DIGEST" module-run design \
    --input "$SPEC" \
    --output "$OUT/design" \
    --num-replicates "$N" \
    --seed "$SEED" \
    --recycle-strategy high/front_loaded

# ---- record what produced this -----------------------------------------------
METRICS=$(find "$OUT/design" -name per_example_metrics.csv | head -1)
{
  echo "date          $(date -Is)"
  echo "design image  $RFD_DIGEST"
  echo "N             $N"
  echo "SEED          $SEED"
  echo "TSR           unset (HER2 novel weights default)"
  echo "spec          $SPEC"
  echo "gpu           $(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader)"
  echo "docker        $(docker --version)"
  echo "metrics       $METRICS"
} > "$OUT/run-settings.txt"
cat "$OUT/run-settings.txt"

# ---- verify ------------------------------------------------------------------
# The CDRH3 is cut out of each designed VH by its flanking framework anchors,
# the same way build_designs.py does it.
got=$(python3 - "$METRICS" <<'PY'
import csv, hashlib, sys
PRE, POST = "YYCSR", "WGQGTLVTVSS"
seqs = []
for row in csv.DictReader(open(sys.argv[1], newline="")):
    vh = next((p.split(":", 1)[1].strip()
               for p in row["designed_sequences"].split(";")
               if p.split(":", 1)[0].strip() == "B"), "")
    i, j = vh.find(PRE), vh.find(POST)
    seqs.append(vh[i + len(PRE):j - 1] if i != -1 and j > i else "")
print(len(seqs), hashlib.md5("\n".join(sorted(seqs)).encode()).hexdigest())
PY
)
n=${got%% *}; md5=${got##* }
echo
echo "designs produced   $n (expected $N)"
echo "CDRH3 fingerprint  $md5"
echo "expected           $EXPECTED_CDRH3_MD5"
if [ "$md5" = "$EXPECTED_CDRH3_MD5" ] && [ "$n" = "$N" ]; then
  echo "REPRODUCED: identical CDRH3 set."
else
  echo "DIFFERS. The seed is fixed, so a difference points at the environment," \
       "most likely a different GPU model or CUDA/driver version — floating-point" \
       "reduction order is not guaranteed identical across them. Compare" \
       "$OUT/run-settings.txt against the values recorded above."
fi
