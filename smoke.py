"""Smoke test: confirms anarcii works end to end and prints the data shape.

Run:  cd ~/antibody-hack && source .venv/bin/activate && python smoke.py
"""

# Trastuzumab (Herceptin) variable regions.
# VERIFY THESE against Thera-SAbDab before using them for calibration.
TRAST_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRAST_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

# IMGT region boundaries, inclusive. Fixed by the IMGT scheme -- this is the
# whole reason for numbering a sequence in the first place.
IMGT_REGIONS = {
    "FR1":  (1, 26),
    "CDR1": (27, 38),
    "FR2":  (39, 55),
    "CDR2": (56, 65),
    "FR3":  (66, 104),
    "CDR3": (105, 117),
    "FR4":  (118, 128),
}


def regions(numbering):
    """Split an anarcii numbering list into IMGT regions.

    `numbering` is a list of ((position, insertion_code), amino_acid) pairs.
    Gap characters '-' are dropped.
    """
    out = {name: "" for name in IMGT_REGIONS}
    for (pos, _ins), aa in numbering:
        if aa == "-":
            continue
        for name, (lo, hi) in IMGT_REGIONS.items():
            if lo <= pos <= hi:
                out[name] += aa
                break
    return out


def main():
    import numpy, pandas, Bio, torch
    print(f"numpy {numpy.__version__} | pandas {pandas.__version__} "
          f"| biopython {Bio.__version__} | torch {torch.__version__}")

    import anarcii
    try:
        from anarcii import Anarcii
    except ImportError:  # import path has differed between versions
        from anarcii.pipeline import Anarcii
    print(f"anarcii {anarcii.__version__}\n")

    model = Anarcii(seq_type="antibody", mode="accuracy",
                    batch_size=2, cpu=True, verbose=False)

    result = model.number({"trast_VH": TRAST_VH, "trast_VL": TRAST_VL})

    print(f"RETURN TYPE: {type(result)}")
    print(f"KEYS: {list(result)}\n")

    for name, rec in result.items():
        print(f"--- {name}")
        print(f"    record keys : {list(rec)}")
        print(f"    chain_type  : {rec['chain_type']}   (H=heavy, K/L=light, F=failed)")
        print(f"    score       : {rec['score']:.3f}")
        print(f"    scheme      : {rec['scheme']}")
        print(f"    query range : {rec['query_start']}-{rec['query_end']}")

        if rec["numbering"] is None:
            print(f"    FAILED      : {rec.get('error')}")
            continue

        print(f"    first 3 numbering entries: {rec['numbering'][:3]}")
        for region, seq in regions(rec["numbering"]).items():
            print(f"    {region:5s} ({len(seq):2d}) {seq}")
        print()

    cdr_h3 = regions(result["trast_VH"]["numbering"])["CDR3"]
    print(f"\nTrastuzumab CDR-H3 = {cdr_h3}")
    print("Should contain WGGDGFYAMDY. If it does, the whole chain works.")


if __name__ == "__main__":
    main()
