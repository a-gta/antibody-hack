#!/usr/bin/env bash
# Submit the two 1ZVY design jobs to JapanFold's RFD4-Proteina.
# Needs: export JAPANFOLD_API_KEY=...   (japanfold.aiand.com/account/token)
#
# Numbering verified by check_1zvy.py against data/1ZVY.cif:
#   - chain A (VHH D3-L11) observed 1-125, no gaps, A1 IS resolved
#   - chain B (lysozyme)   observed 1-129, no gaps
#   - label_asym_id/label_seq_id == auth numbering for both chains, so the
#     contigs below are valid as written (RFD4 reads the label_* numbering)
#   - IMGT CDR3 105-117 == A97-A116, flanked by Cys A96 and Trp A117
set -euo pipefail

BASE=https://api.japanfold.aiand.com
: "${JAPANFOLD_API_KEY:?set JAPANFOLD_API_KEY first}"
H=(-H "Authorization: Bearer $JAPANFOLD_API_KEY" -H 'Content-Type: application/json')

# RFD4 fetches RCSB entries itself, so no paste is needed: "pdb:1ZVY".
# Swap for  "structure": "$(python3 -c '...read data/1ZVY.cif...')"  to pin the text.

submit () {  # $1=label  $2=json
  echo "--- $1"
  curl -s "${H[@]}" -X POST "$BASE/v1/designs" -d "$2"
  echo
}

# 1. New binder against lysozyme. Keeps all of chain B, designs a fresh
#    110-130 residue chain. 129 + 130 = 259 residues, under the 512 cap.
submit "new binder" '{
  "protocol":"rfd4-binder",
  "name":"1zvy-denovo-binder",
  "structure":"pdb:1ZVY",
  "contig":"B1-129, /0, 110-130",
  "params":{"num_designs":4,"recipe":"thinking_medium/front_loaded","nsteps":60}
}'

# 2. CDR3 loop redesign on the VHH. One continuous chain A: framework up to
#    Cys96, 20 designed residues in place of CDR3, framework from Trp117.
#    Antigen kept as its own chain after the /0.  Use A1-96 to keep Asp1 too.
submit "CDR3 redesign" '{
  "protocol":"rfd4-scaffold",
  "name":"1zvy-cdr3-redesign",
  "structure":"pdb:1ZVY",
  "contig":"B1-129, /0, A2-96, 20, A117-125",
  "params":{"num_designs":4,"recipe":"thinking_high/front_loaded","nsteps":60}
}'

# Poll:      curl -s "${H[@]}" "$BASE/v1/jobs/<id>?wait=60"
# Results:   curl -s "${H[@]}" "$BASE/v1/jobs/<id>/results"
# Each design carries designed_sequence, so no inverse-folding step is needed.
# Check num_chains == 2 and compare motif_rmsd within the job, not against 0.
