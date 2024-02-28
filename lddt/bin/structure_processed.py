"""
-*- coding: utf-8 -*-

@Author : Bu_Fan
@Time : 2023/12/16 11:49
@File : strutcture_processed.py
@aim:
"""
import ast
from rna_tools.rna_tools_lib import *
from Bio import PDB
import os
import subprocess
import glob




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
        print("Did not create new folder.")
        print(f"{new_dir_name} already exists.")
        
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
        #print('the processed result is',rnatools_remarks)

        casp15rna_remarks = "REMARK 250 Model processed with RNA Puzzles\n" # TODO: Think about text
        rnatools_remarks = [line + "\n" for line in rnatools_remarks]

        lines = [line + "\n" for line in s.lines]

        output_fn = suffix_and_write("processed", input, subdir="processed")

        with open(output_fn, 'w') as f:
            f.write(casp15rna_remarks)
            f.writelines(rnatools_remarks)
            f.writelines(lines)
    except Exception as e:
        print(f"Error processing file : {e}")


def run_bash_command(cmd):
    """
    Run a bash command. Print the command and the output.

    Parameters:
    cmd (str): The bash command to run.

    Returns:
    subprocess.CompletedProcess: The result of the command.
    """
    print(f"command: {cmd}")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    # print(f"stdout: {result.stdout}")
    # print(f"stderr: {result.stderr}\n")
    

    return result.stdout

if __name__=="__main__":
    #检查路径是否存在
    directory_path = '/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/Puzzles/original/deleteG'
    if os.path.exists(directory_path):
        predicted_list = os.listdir(directory_path)
        for i in predicted_list:
            if i.endswith('.pdb'):
                file_path = os.path.join(directory_path, i)

                preprocess(file_path)
    else:
         print(f"Directory does not exist: {directory_path}")
    # input_pdb="/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/Puzzles/original/PZ38_Dfold_1_deleteA.pdb"
    # preprocess(input_pdb)
