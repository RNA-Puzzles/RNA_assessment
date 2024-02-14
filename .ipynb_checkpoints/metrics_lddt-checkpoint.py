"""
-*- coding: utf-8 -*-

@Author : Bu_Fan
@Time : 2023/12/23 17:01
@File : metrics_lddt.py
@aim:extract the lddt part from metrics
"""
from normal import is_valid_pdb
from normal import rel_path
from normal import run_bash_command
from normal import make_list_of_pdbs
from normal import get_lines
from normal import write_lines
from rna_tools.rna_tools_lib import *
import os
import shutil
import subprocess
import glob
import json
import csv
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

def custom_output_function(*args):
    try:
        print(*args)
    except Exception as e:
        pass  # 忽略异常，不打印任何东西

ic.configureOutput(outputFunction=custom_output_function)

# TODO: Consider using decorations instead
class Metric:
    project_path = os.path.abspath(os.getcwd())

    def __init__(self, name):
        self.name = name
        self.binary_path = None

        os.makedirs("scores", exist_ok=True)

    def chdir_run_folder_ref_mod(self, run_name, r, m):
        run_folder_ref_mod = os.path.abspath(f"runs/{run_name}/{r[0]}/{r[0]}.{m[0]}/{self.name}")

        # Mkdir
        os.makedirs(run_folder_ref_mod, exist_ok=True)
        os.chdir(run_folder_ref_mod)

    def chdir_run_folder_ref(self, run_name, r):
        run_folder_ref = os.path.abspath(f"runs/{run_name}/{r[0]}/{self.name}")

        # Mkdir
        os.makedirs(run_folder_ref, exist_ok=True)
        os.chdir(run_folder_ref)

    def calculate(self, reference, model, run_name="default"):
        reference = os.path.abspath(reference)
        model = os.path.abspath(model)

        assert is_valid_pdb(reference) and is_valid_pdb(model), "reference and model PDB must be valid."

        self.cleanup_stack = []
        #os.path.basename(reference)形成文件名
        r = os.path.splitext(os.path.basename(reference))#形成包含文件名和扩展名的元祖
        m = os.path.splitext(os.path.basename(model))

        self.chdir_run_folder_ref_mod(run_name, r, m)#更改当前目录到新创建的目录下

        # LGA does not take paths outside the current working directory, so workaround needed
        # Create a folder called references and models
        #"".join(r) 和 "".join(m) 分别将文件名和扩展名重新组合成一个字符串，作为目标文件名。例如，如果 r 是 ("file", ".xyz")，"".join(r) 将返回 "file.xyz"。
        shutil.copy(rel_path(reference), "".join(r))
        shutil.copy(rel_path(model), "".join(m))
        #清理临时文件
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(r))))
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(m))))

    def cleanup(self):
        """
        清除创建的临时文件并返回到项目根目录。
        """
        while self.cleanup_stack:
            os.remove(self.cleanup_stack[0])
            self.cleanup_stack.pop(0)

        os.chdir(self.project_path)

    def calc_bulk(self, reference, models, force=False, run_name="default"):
        reference = os.path.abspath(reference)
        r = os.path.splitext(os.path.basename(reference))

        models_abs = [os.path.abspath(model) for model in models]

        self.chdir_run_folder_ref(run_name, r)

        assert is_valid_pdb(reference), f"reference PDB {reference} must be valid."

        self.cleanup_stack = []
        r = os.path.splitext(os.path.basename(reference))
        shutil.copy(rel_path(reference), "".join(r))
        self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), "".join(r))))

        for model in models_abs:
            m = os.path.splitext(os.path.basename(model))

            shutil.copy(rel_path(model), "".join(m))
            self.cleanup_stack.append(os.path.abspath(os.path.join(os.getcwd(), ''.join(m))))

    def consol_bulk(self, target_name):
        raise NotImplementedError("This functionality has not yet been implemented.")

    def consolidate(self, target_name):
        raise NotImplementedError("This functionality has not yet been implemented.")

    def process_out_file(self, out_file):
        raise NotImplementedError("This functionality has not yet been implemented.")



class LDDT(Metric):

    def __init__(self):
        super().__init__("lddt")

        lddt_path = os.path.abspath("/home/bufan/RNA_assessment/lddt")
        self.binary_path = os.path.abspath(lddt_path)

    def consolidate(self, target_name):
        # Define the headers for the CSV file
        headers = ['reference', 'model', 'lddt']

        # Open the CSV file for writing
        with open(f"scores/{self.name}.{target_name}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)

            # Iterate through all the JSON files in the directory
            for json_file in glob.glob(f"runs/*/*/*/{self.name}/{self.name}.*.txt"):
                print(f"current file:{self.name}.*.txt")
                # Open the JSON file and extract the relevant data
                with open(json_file) as jf:
                    data = json.load(jf)
                    reference = data['trg_file']
                    model = data['mdl_file']
                    lddt = data['lDDT']

                    # Write the data to the CSV file
                    writer.writerow([reference, model, lddt])

    def calculate(self, reference, model, force=False):
        super().calculate(reference, model)

        r = os.path.splitext(os.path.basename(reference))
        m = os.path.splitext(os.path.basename(model))
        output = f"{self.name}.{r[0]}.{m[0]}.txt"

        # if force == False and os.path.exists(output) and os.path.getsize(output) > 0:
        #     print(f"{output} exists. Skipping calculations...")

        # else:
            #print('processing test :2023-12-23')
        c = f"{self.binary_path}/ema --mount {os.getcwd()} {self.binary_path}/bin/complex_lddt_no_stereocheck.py {''.join(m)} {''.join(r)} {output}"
        
        lddt_out=float(run_bash_command(c).strip())
        
        #print(f"::::{lddt_out}")
        ## Cleanup
        self.cleanup()
        return lddt_out
        
        
       
