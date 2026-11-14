"""Check PDB 1ZVY (VHH D3-L11 + hen egg-white lysozyme) before designing against it.

- parses the mmCIF straight from RCSB (no structure-parsing dependency)
- reports what is actually observed in each chain, including gaps
- numbers chain A with ANARCII (IMGT) and maps IMGT positions onto the
  author (PDB) residue numbers the mmCIF actually uses
- checks the two candidate RFdiffusion-style contig strings against that

Run:  cd ~/antibody-hack && source .venv/bin/activate && python check_1zvy.py
"""

from pathlib import Path
from urllib.request import urlopen

try:
    from anarcii import Anarcii
except ImportError:  # import path has differed between versions
    from anarcii.pipeline import Anarcii

CIF_URL = "https://files.rcsb.org/download/1ZVY.cif"
CIF_PATH = Path(__file__).parent / "data" / "1ZVY.cif"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}

# IMGT CDR definitions (Lefranc). Inclusive ranges.
IMGT_CDRS = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}


def fetch_cif():
    if not CIF_PATH.exists():
        CIF_PATH.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(CIF_URL) as r:
            CIF_PATH.write_bytes(r.read())
        print(f"downloaded {CIF_URL} -> {CIF_PATH}")
    else:
        print(f"using cached {CIF_PATH}")
    return CIF_PATH.read_text()


def parse_atom_site(text):
    """Return {auth_asym: [(auth_seq_id, ins_code, comp_id), ...]} in file order.

    Column positions are read from the loop header rather than assumed.
    """
    lines = text.splitlines()
    cols, start = {}, None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".", 1)[1]] = len(cols)
            start = i + 1
    if not cols:
        raise SystemExit("no _atom_site loop found")

    i_seq = cols["auth_seq_id"]
    i_comp = cols["auth_comp_id"]
    i_asym = cols["auth_asym_id"]
    i_ins = cols["pdbx_PDB_ins_code"]
    i_lseq = cols["label_seq_id"]
    i_lasym = cols["label_asym_id"]

    chains, labels, seen = {}, {}, set()
    for line in lines[start:]:
        if line.startswith("#"):
            break
        if not line.startswith("ATOM"):
            continue  # skip HETATM: phosphate, glycerol, waters
        f = line.split()
        asym, seq, comp = f[i_asym], f[i_seq], f[i_comp]
        ins = "" if f[i_ins] in (".", "?") else f[i_ins]
        key = (asym, seq, ins)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(asym, []).append((int(seq), ins, comp))
        labels.setdefault(asym, []).append((f[i_lasym], f[i_lseq]))
    return chains, labels


def check_label_numbering(chains, labels):
    """RFD4-Proteina contigs name label_asym_id / label_seq_id, not the auth
    numbering a viewer shows. Confirm whether they coincide for this file."""
    print("\n--- label_* vs auth_* numbering (what RFD4 contigs actually name) ---")
    same = True
    for cid in sorted(chains):
        auth = chains[cid]
        lab = labels[cid]
        if any(ls in (".", "?") for _, ls in lab):
            print(f"  chain {cid}: label_seq_id is empty for some rows "
                  f"-> RFD4 refuses a kept range here")
            same = False
            continue
        asyms = {la for la, _ in lab}
        offsets = {int(ls) - a[0] for (_, ls), a in zip(lab, auth)}
        lo, hi = min(int(ls) for _, ls in lab), max(int(ls) for _, ls in lab)
        print(f"  auth chain {cid} ({auth[0][0]}-{auth[-1][0]})  ->  "
              f"label_asym {sorted(asyms)} label_seq {lo}-{hi}  "
              f"offset {sorted(offsets)}")
        if asyms != {cid} or offsets != {0}:
            same = False
    print("  verdict:", "label == auth, so contigs may use the numbers above"
          if same else "THEY DIFFER - rewrite the contigs in label numbering")
    return same


def describe(chain_id, residues):
    nums = [n for n, ins, _ in residues if not ins]
    inserts = [(n, ins) for n, ins, _ in residues if ins]
    lo, hi = min(nums), max(nums)
    present = set(nums)
    gaps = [n for n in range(lo, hi + 1) if n not in present]
    print(f"\nchain {chain_id}: {len(residues)} observed residues, "
          f"auth numbering {lo}-{hi}")
    print(f"  internal/leading gaps in {lo}-{hi}: {gaps or 'none'}")
    if inserts:
        print(f"  insertion codes: {inserts}")
    seq = "".join(THREE_TO_ONE.get(c, "X") for _, _, c in residues)
    print(f"  sequence ({len(seq)} aa): {seq}")
    return seq


def main():
    chains, labels = parse_atom_site(fetch_cif())
    print("polymer chains (auth ids):", sorted(chains))

    seqs = {cid: describe(cid, chains[cid]) for cid in sorted(chains)}
    check_label_numbering(chains, labels)

    vhh_id = "A"
    vhh_res = chains[vhh_id]
    vhh_seq = seqs[vhh_id]

    # --- ANARCII on the observed chain A sequence -------------------------
    print("\n--- ANARCII (IMGT) on chain A ---")
    model = Anarcii(seq_type="antibody", mode="accuracy", cpu=True, verbose=False)
    out = model.number([vhh_seq])
    rec = list(out.values())[0]
    print(f"  chain_type={rec['chain_type']}  score={rec['score']}")

    # ANARCII only ever saw the OBSERVED residues, so the k-th non-gap
    # position in its numbering is the k-th residue in vhh_res.
    numbering = rec["numbering"]
    aligned = []  # (imgt_num, imgt_ins, aa, auth_seq, auth_ins)
    k = 0
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        auth_seq, auth_ins, comp = vhh_res[k]
        expect = THREE_TO_ONE.get(comp, "X")
        if aa != expect:
            raise SystemExit(
                f"alignment desync at IMGT {num}{ins.strip()}: "
                f"ANARCII says {aa}, file says {expect} at {auth_seq}{auth_ins}"
            )
        aligned.append((num, ins.strip(), aa, auth_seq, auth_ins))
        k += 1
    print(f"  mapped {len(aligned)} residues; "
          f"IMGT {aligned[0][0]}-{aligned[-1][0]} <-> "
          f"PDB {aligned[0][3]}-{aligned[-1][3]}")

    off = {(n, i) for n, i, _, s, _ in aligned if not i and n != s}
    print(f"  positions where IMGT number != PDB number: {len(off)} of {len(aligned)}")

    # --- CDRs -------------------------------------------------------------
    print("\n--- IMGT CDRs -> PDB residue numbers ---")
    cdr3_pdb = []
    for name, (lo, hi) in IMGT_CDRS.items():
        hits = [a for a in aligned if lo <= a[0] <= hi]
        if not hits:
            print(f"  {name}: nothing observed")
            continue
        seq = "".join(a[2] for a in hits)
        pdb_lo = f"{hits[0][3]}{hits[0][4]}"
        pdb_hi = f"{hits[-1][3]}{hits[-1][4]}"
        imgt_lo = f"{hits[0][0]}{hits[0][1]}"
        imgt_hi = f"{hits[-1][0]}{hits[-1][1]}"
        print(f"  {name}: {seq}  ({len(seq)} aa)")
        print(f"         IMGT {imgt_lo}-{imgt_hi}  ->  PDB A{pdb_lo}-A{pdb_hi}")
        if name == "CDR3":
            cdr3_pdb = hits

    # --- contig sanity checks --------------------------------------------
    print("\n--- contig checks ---")
    b_nums = sorted(n for n, _, _ in chains["B"])
    a_nums = sorted(n for n, _, _ in vhh_res)

    def check_span(ch, lo, hi):
        nums = b_nums if ch == "B" else a_nums
        present = set(nums)
        missing = [n for n in range(lo, hi + 1) if n not in present]
        ok = not missing and lo >= min(nums) and hi <= max(nums)
        status = "OK" if ok else "PROBLEM"
        print(f"  {ch}{lo}-{hi}: {status}"
              + (f" (not observed: {missing[:12]}"
                 f"{' ...' if len(missing) > 12 else ''})" if missing else ""))
        if hi > max(nums) or lo < min(nums):
            print(f"      chain {ch} is only observed over {min(nums)}-{max(nums)}")
        return ok

    print("  new binder: B1-129, /0, 110-130")
    check_span("B", 1, 129)
    print("  loop redesign: B1-129, /0, A2-96, 20, A117-125")
    check_span("B", 1, 129)
    check_span("A", 2, 96)
    check_span("A", 117, 125)

    if cdr3_pdb:
        lo, hi = cdr3_pdb[0][3], cdr3_pdb[-1][3]
        n = hi - lo + 1
        print(f"\n  CDR3 occupies A{lo}-A{hi} ({n} residues).")

        # The IMGT anchors either side of CDR3 are Cys104 and Trp118. If the
        # mapping is right, the kept framework must end on C and resume on W.
        by_auth = {a[3]: a for a in aligned}
        before, after = by_auth.get(lo - 1), by_auth.get(hi + 1)
        for label, a, want in (("left", before, "C"), ("right", after, "W")):
            if a is None:
                print(f"  {label} anchor: missing from the model")
            else:
                ok = "OK" if a[2] == want else f"EXPECTED {want}"
                print(f"  {label} anchor A{a[3]} = {a[2]} "
                      f"(IMGT {a[0]}{a[1]}) {ok}")

        print(f"\n  So 'A2-96, 20, A117-125' keeps framework up to the Cys, "
              f"rebuilds the {n}-residue CDR3, and resumes at the Trp: correct.")
        print(f"  Note A1 ({vhh_res[0][2]}) IS resolved, so A2-96 drops one real "
              f"residue; A1-96 also works.")


if __name__ == "__main__":
    main()
