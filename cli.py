#!/usr/bin/env python3

import argparse
from rdkit import Chem
from hfrpkg.runner import run_reaction

def parse_smiles_list(arg):
    if not arg:
        return None
    if isinstance(arg, list):
        raw = ",".join(arg)
    else:
        raw = arg
    smiles_list = raw.split(",")
    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise argparse.ArgumentTypeError(f"Invalid SMILES: {smi}")
        mols.append(mol)
    return mols

def main():
    parser = argparse.ArgumentParser(description="Run reaction balancing (HFR)")
    parser.add_argument("action_type", choices=["write", "view", "count"], help="Action to perform")
    parser.add_argument("reaction_type", choices=["isogyric", "isodesmic", "homodesmotic", "hypohomodesmotic"], help="Reaction type")
    parser.add_argument("input", type=str, help="Input molecule SMILES")
    parser.add_argument("--lhs", nargs="+", help="Required LHS molecules (comma-separated or list)")
    parser.add_argument("--rhs", nargs="+", help="Required RHS molecules (comma-separated or list)")
    parser.add_argument("--substruct", nargs="+", help="Substructures to replace")
    parser.add_argument("--replacement", nargs="+", help="Replacement structures")
    parser.add_argument("--outfolder", help="Output outfolder (only required for write)")
    parser.add_argument("--m", help="[METHOD]")
    parser.add_argument("--b", help="[BASIS]")
    parser.add_argument("--s", help="Software to use ['g'-Gaussian, 'o'-Orca, 'p'-Psi4] Default:Gaussian")

    args = parser.parse_args()

    lhs = parse_smiles_list(args.lhs)
    rhs = parse_smiles_list(args.rhs)
    substruct = parse_smiles_list(args.substruct)
    replacement = parse_smiles_list(args.replacement)

    
    run_reaction(
        action_type=args.action_type,
        reaction_type=args.reaction_type,
        input_smiles=args.input,
        lhs=lhs,
        rhs=rhs,
        substruct=substruct,
        replacement=replacement,
        outfolder=args.outfolder,
        method = args.m,
        basis = args.b,
        extension = args.s
    )

if __name__ == "__main__":
    main()

