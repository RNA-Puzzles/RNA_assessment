"""
-*- coding: utf-8 -*-

@Author : Bu_Fan
@Time : 2023/12/16 11:49
@File : normal.py.py
@aim:
"""
import ast

import pandas as pd
from rna_tools.rna_tools_lib import *
from Bio import PDB
import os
import subprocess
import glob
import itertools
import csv


def create_dict_from_csv(csv_file):
    df=pd.read_csv(csv_file)
    result_dict={}
    for index,row in df.iterrows():
        key=row['new_model_name_pdb']
        chain_mp_dict=ast.literal_eval((row['chain_mp']))
        value=[row['new_extract_name_pdb'],chain_mp_dict]
        result_dict[key]=value
    return result_dict

def rel_path(path):
    """
    Get the relative path of a file.

    Parameters:
    path (str): The path to the file.

    Returns:
    str: The relative path of the file.
    """
    return os.path.relpath(path, os.getcwd())

def preprocess_rec(target_path):
    """
    Preprocess all PDB files in a directory.

    Parameters:
    target_path (str): Path to directory containing PDB files.

    Returns:
    None
    """
    for type in ["references", "models"]:
        subdir = os.path.join(target_path, type)

        if not os.path.isdir(subdir):
            print(f"{subdir} does not exist in {target_path}.")
            return

        pdb_files = [f for f in os.listdir(subdir) if f.endswith(".pdb")]
        for pdb_file in pdb_files:
            pdb_file_path = os.path.join(subdir, pdb_file)
            preprocess(pdb_file_path)

def is_valid_pdb(fn):
    """
    Check if a file is a valid PDB file.

    Parameters:
    fn (str): The path to the file to check.

    Returns:
    bool: True if the file is a valid PDB file, False otherwise.
    """

    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)

    try:
        # Attempt to parse the PDB file
        structure = parser.get_structure('X', fn)

        # If the structure has no models, it is invalid
        if len(structure) == 0:
            print("a")
            return False

        # If the structure has no chains, it is invalid
        elif len(list(structure.get_chains())) == 0:
            return False

        # If the structure has no residues, it is invalid
        elif len(list(structure.get_residues())) == 0:
            return False

        else:
            return True

    except Exception:
        print(f"{fn} is not a valid PDB.")
        return False


def suffix_and_write(suffix, fn, subdir=""):
    """
    Helper function to write a file to the same directory as fn. fn will have an additional suffix name

    Parameters:
    suffix (str): suffix that will append the filename
    fn (str): the path to a specific file

    Returns:
    The fn suffixed with the suffix parameter
    """

    # Split the file path into its directory and filename components
    dir_name, base_name = os.path.split(fn)

    # Create a new filename with a suffix
    new_base_name = f'{os.path.splitext(base_name)[0]}_{suffix}{os.path.splitext(base_name)[1]}'

    # Create new subdir
    new_dir_name = os.path.join(dir_name, subdir)

    try:
        os.mkdir(new_dir_name)
    except FileExistsError:
        #print("Did not create new folder.")
        #print(f"{new_dir_name} already exists.")
        print('')

    output = os.path.join(new_dir_name, new_base_name)
    return output



def preprocess(input):
    """
    Preprocess a PDB using rna-tools and clean up nonstandard nucletides.

    Parameters:
    fn (str): Path to .pdb file

    Returns:
    bool: True if the PDB file is valid, False otherwise.
    """
    try:
        print(f"Processing {input}")
        # Checks if PDB is valid
        assert is_valid_pdb(input) == True, "pdb_structure must be valid"

        s = RNAStructure(input)
        rnatools_remarks = s.get_rnapuzzle_ready()
        print('the processed result is',rnatools_remarks)

        casp15rna_remarks = "REMARK 250 Model processed with casp-rna\n" # TODO: Think about text
        rnatools_remarks = [line + "\n" for line in rnatools_remarks]

        lines = [line + "\n" for line in s.lines]

        output_fn = suffix_and_write("processed", input, subdir="processed")

        with open(output_fn, 'w') as f:
            f.write(casp15rna_remarks)
            f.writelines(rnatools_remarks)
            f.writelines(lines)
    except Exception as e:
        print(f"Error processing file : {e}")

def make_list_of_pdbs(input_dir):
    """
    Write a file containing PDB files in a directory. Necessary for calculating GDT-TS.

    Parameters:
    input_dir (str): Path to directory containing PDB files.

    Returns:
    None
    """
    prev_path = os.path.abspath(os.getcwd())

    os.chdir(input_dir)
    print(f"Generating list for path: {os.path.abspath(os.getcwd())}")
    run_bash_command("ls -1 *.pdb | sed -e 's/\.pdb$//' > list")

    os.chdir(prev_path)

def get_lines(file_pattern, m, n):
    """
    Get lines from files that match the specified pattern. m and n are the first and last line numbers to include.

    Parameters:
    file_pattern (str): A glob pattern to match files.
    m (int): The first line number to include.
    n (int): The last line number to include. If -1, all lines after m will be included.

    Returns:
    list: A list of lines from the files that match the pattern.
    """
    lines = []
    # iterate through each file that matches the pattern
    for filename in glob.glob(file_pattern):
        with open(filename, 'r') as f:
            # iterate through each line in the file
            for i, line in enumerate(f):
                # check if the line number is within the specified range
                if i+1 >= m and (n == -1 or i+1 <= n):
                    lines.append(line.strip())
    return lines

def write_lines(file, lines):
    """
    Write the lines to the specified file.

    Parameters:
    file (str): The path to the file to write.
    lines (list): A list of lines to write.

    Returns:
    None
    """
    print(f"lines: {lines}")

    with open(file, "w") as outfile:
        outfile.flush()

        outfile.write("reference,model,inf_all,inf_stack,inf_WC,inf_nWC,sns_WC,ppv_WC,sns_nWC,ppv_nWC")

        # Create a csv writer object
        writer = csv.writer(outfile)

        # Write each second line to the output file
        for line in lines:
            print(f"line: {line}")
            outfile.write("\n")
            outfile.write(line)


def run_bash_command(cmd):
    """
    Run a bash command. Print the command and the output.

    Parameters:
    cmd (str): The bash command to run.

    Returns:
    subprocess.CompletedProcess: The result of the command.
    """
    #print(f"command: {cmd}")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    #print(f"stdout: {result.stdout}")
    #return result.stdout
    #print(f"stderr: {result.stderr}\n")

    return result.stdout

if __name__=="__main__":
    #检查路径是否存在
    # directory_path = '/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/Puzzles/original/deleteG'
    # if os.path.exists(directory_path):
    #     predicted_list = os.listdir(directory_path)
    #     for i in predicted_list:
    #         if i.endswith('.pdb'):
    #             file_path = os.path.join(directory_path, i)
    #
    #             preprocess(file_path)
    # else:
    #      print(f"Directory does not exist: {directory_path}")
     input_pdb="/home/bufan/RNA_assessment/example/14_ChenPostExp_2.pdb"
     preprocess(input_pdb)