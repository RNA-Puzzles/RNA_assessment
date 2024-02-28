"""
-*- coding: utf-8 -*-

@Author : Bu_Fan
@Time : 2023/12/23 17:01
@File : metrics_lddt.py
@aim:extract the lddt part from metrics
"""

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


def summary_to_csv():
    # Define the headers for the CSV file
    headers = ['reference', 'model', 'lddt']

        # Open the CSV file for writing
    with open(f"scores/20240120summary_all.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)

            # Iterate through all the JSON files in the directory
        for json_file in glob.glob(f"runs/*/*/*/lddt/lddt.*.txt"):
            #print(f"current file:lddt.*.txt")
                # Open the JSON file and extract the relevant data
            with open(json_file) as jf:
                data = json.load(jf)
                reference = data['trg_file']
                model = data['mdl_file']
                lddt = data['lDDT']

                # Write the data to the CSV file
                writer.writerow([reference, model, lddt])
