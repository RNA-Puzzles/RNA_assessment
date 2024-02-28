"""
-*- coding: utf-8 -*-

@Author : Bu_Fan
@Time : 2023/12/23 17:01
@File : metrics_lddt.py
@aim:extract the lddt part from metrics
"""

from structure_processed import is_valid_pdb
from structure_processed import run_bash_command
from rna_tools.rna_tools_lib import *
from Bio import PDB
import os






class LDDT:
   

    def __init__(self):
        current_script_path = os.path.abspath(__file__)
        
        # Get the directory where the current script is located
        self.binary_path = os.path.dirname(current_script_path)
       

    def calculate(self, reference, model, chain_mapping, force=False):
        reference = os.path.abspath(reference)
        model = os.path.abspath(model)

        assert is_valid_pdb(reference) and is_valid_pdb(model), "reference and model PDB must be valid."

        pdb_path = os.path.dirname(reference)
        c = f"{os.path.dirname(self.binary_path)}/ema --mount {pdb_path} {self.binary_path}/complex_lddt_no_stereocheck.py {model} {reference} '{chain_mapping}'"
        lddt_out = float(run_bash_command(c).strip())


        return lddt_out

