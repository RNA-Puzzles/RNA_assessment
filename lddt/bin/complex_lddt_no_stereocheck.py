import argparse
import os
import json

from ost import io, seq
from ost.mol.alg.lddt import lDDTScorer
import pandas as pd
import ast

def _parse_args():
    desc = ("Computes global and per-residue lDDT score on monomer models "
            "without performing a stereo-chemistry check")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("mdl_file", help="model file in PDB format")
    parser.add_argument("trg_file", help="target file in PDB format")
    parser.add_argument("chain_mapping", help="the chain mapping between mdl and trg file")
    return parser.parse_args()

def main():
    args = _parse_args()
    mdl = io.LoadPDB(args.mdl_file)
    trg = io.LoadPDB(args.trg_file)
    chain_mapping=json.loads(args.chain_mapping)
 
    lddt_scorer = lDDTScorer(trg)
    lddt_score, _ = lddt_scorer.lDDT(mdl, local_lddt_prop="lddt",chain_mapping=chain_mapping)

    local_scores = dict()
    for r in mdl.residues:
        if r.HasProp("lddt"):
            local_scores[r.GetNumber().GetNum()] = r.GetFloatProp("lddt")
        else:
            local_scores[r.GetNumber().GetNum()] = None

    # create and dump output
    json_summary = {
        "mdl_file": args.mdl_file,
        "trg_file": args.trg_file,
        "lDDT": lddt_score,
        "local_lDDT": local_scores
    }
   
    print(json_summary['lDDT'])

    

if __name__ == '__main__':
    main()

