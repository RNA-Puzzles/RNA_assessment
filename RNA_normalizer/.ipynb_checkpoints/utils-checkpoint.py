#  utils.py
#  
#  Copyright 2017 Chichau Miau <zmiao@ebi.ac.uk>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  util functions
from operator import attrgetter
import os

from .msgs import *

def command( cmd ):
    ret_code = os.system( cmd )
    if( ret_code != 0 ):
        print ("FATAL", "Command '%s' ended with error code '%s'" %(cmd, ret_code) )
    
class Result:
    def __init__(self, problem, original_file, lab, result ):
        self.problem = problem
        self.file_original = original_file
        self.lab = lab
        self.result = result
        self.file = "%d_%s_%s.pdb" %(problem, lab, result)

def read_results_list( problem, fname ):
    results = {}
    
    rows = filter( lambda row: len(row) == 3, map( lambda row: row.split(), open( fname ).read().strip().split( "\n" ) ) )
    
    for row in rows:
        results[row[0]] = Result( problem, row[0], row[1], int(row[2]) )
    
    return( results )

def get_index_file( pdb_file, pdb_result_file="" ):
    index_file = "%s.%s.index" %(pdb_file.replace( ".pdb", "" ), pdb_result_file.replace( ".pdb", "" ))

    if (pdb_result_file is "") or (not os.path.isfile( index_file )):
        index_file = pdb_file.replace( ".pdb", ".index" )
    
    if( not os.path.isfile( index_file ) ):
        show( "INFO", "INDEX SKIPPED! '%s' skipped for '%s'." %(index_file, pdb_file) )
        index_file= None
    else:
        show( "INFO", "INDEX FOUND! '%s' for '%s'." %(index_file, pdb_file) )

    return( index_file )

def molprobity_parse(f, evals):
	for line in file(f):
		line = line.strip()
		
		if( (line != "") and ("#" not in line) ):
			data = line.split( ":" )
			
			id = data[0].strip( "FH.pdb" ).split( "_" )
			print(data[0]) 
			
			if( len(id) == 3 and id[1] != 'solution'):
				eval = find_eval( evals, int(id[0]), id[1], int(id[2]) )
				if eval is None:
					print ('Cannot find %s'%data[0])
				else:
					eval.clashscore = float(data[8])
			else:
				print ("skip: ", line)

class Eval:
    def __init__(self, problem=0, original="", lab="", result="", result_fit=""):
        self.ok = False
        
        self.problem = problem
        self.original = original
        self.lab = lab
        self.result = result
        #self.result_fit = result_fit
        
        self.rmsd = 1e100
        #~ self.rmsb = 1e100
        self.pvalue = 1e100
        self.DI_ALL = 1e100
        self.INF_ALL = 0.0
        self.INF_WC = 0.0
        self.INF_NWC = 0.0
        self.INF_STACK = 0.0
        self.clashscore = 1e100
        self.best_sol_ndx = -1
        self.mcq=1e100
        self.gdt=1e100
    
    def parse(self, row):
        result = True
        data = row.split()
        
        if( len(data) == 15 ):
            self.problem = int(data[0])
            self.original = data[1]
            self.lab = data[2]
            self.result = int(data[3])
            
            self.rmsd = float(data[4])
            self.pvalue = float(data[5])
            self.DI_ALL = float(data[6])
            self.INF_ALL = float(data[7])
            self.INF_WC = float(data[8])
            self.INF_NWC = float(data[9])
            self.INF_STACK = float(data[10])
            self.clashscore = float(data[11])
            self.mcq=float(data[12])
            self.gdt=float(data[13])
            self.best_sol_ndx = int(data[14])
            
            self.rmsd_rank = 0
            self.pvalue_rank = 0
            self.DI_ALL_rank = 0
            self.INF_ALL_rank = 0
            self.INF_WC_rank = 0
            self.INF_NWC_rank = 0
            self.INF_STACK_rank = 0
            self.clashscore_rank = 0
            self.mcq_rank=0
            self.gdt_rank=0
        else:
            result = False

        return result
    
    def file_name(self):
        return( "%d_%s_%d" %(self.problem, self.lab, self.result) )
    
    def set_rank(self, attr, rank):
        if( attr == "rmsd" ):
            self.rmsd_rank = rank
        elif( attr == "pvalue" ):
            self.pvalue_rank = rank
        elif( attr == "DI_ALL" ):
            self.DI_ALL_rank = rank
        elif( attr == "INF_ALL" ):
            self.INF_ALL_rank = rank
        elif( attr == "INF_WC" ):
            self.INF_WC_rank = rank
        elif( attr == "INF_NWC" ):
            self.INF_NWC_rank = rank
        elif( attr == "INF_STACK" ):
            self.INF_STACK_rank = rank
        elif( attr == "clashscore" ):
            self.clashscore_rank = rank
        elif( attr == "mcq" ):
            self.mcq_rank = rank
        elif( attr == "gdt" ):
            self.gdt_rank = rank
        else:
            show( "FATAL", "Can't set attribute '%s' in class 'Eval'" %(attr) )
    
    def __str__(self):
        s = []
        s += [str(self.problem)]
        s += [str(self.original)]
        s += [self.lab]
        s += ["%d" %self.result]
        s += ["%7.3f" %self.rmsd]
        s += ["%.3e" %self.pvalue]
        s += ["%7.3f" %self.DI_ALL]
        s += ["%7.3f" %self.INF_ALL]
        s += ["%7.3f" %self.INF_WC]
        s += ["%7.3f" %self.INF_NWC]
        s += ["%7.3f" %self.INF_STACK]
        s += ["%7.3f" %self.clashscore]
        s += ["%7.3f" %self.mcq]
        s += ["%7.3f" %self.gdt]
        s += ["%d" %self.best_sol_ndx]
        return " ".join( s )

    def rankline(self):
        print (self.lab,self.result)
        s = []
        s += [self.lab.replace('PostExp','_PostExp').replace('PreExp','_PreExp')]
        s += ["%d" %self.result]
        s += ["%.3f" %self.rmsd]
        s += ["%d" %self.rmsd_rank]
        s += ["%.3f" %self.DI_ALL]
        s += ["%d" %self.DI_ALL_rank]
        s += ["%.3f" %self.INF_ALL]
        s += ["%d" %self.INF_ALL_rank]
        s += ["%.3f" %self.INF_WC]
        s += ["%d" %self.INF_WC_rank]
        s += ["%.3f" %self.INF_NWC]
        s += ["%d" %self.INF_NWC_rank]
        s += ["%.3f" %self.INF_STACK]
        s += ["%d" %self.INF_STACK_rank]
        s += ["%.3f" %self.clashscore]
        s += ["%d" %self.clashscore_rank]
        s += ["%.3f" %self.mcq]
        s += ["%d" %self.mcq_rank]
        s += ["%.3f" %self.gdt]
        s += ["%d" %self.gdt_rank]
        s += ["%.3e" %self.pvalue]
        s += ["%d" %self.best_sol_ndx]
        return "\t".join( s )

def find_eval( evals, problem, lab, result ):
    for eval in evals:
        if( eval.problem == problem and eval.lab == lab and eval.result == result ):
            return eval

def load_evals_list( fname ):
    evals = []
    
    rows = open( fname ).read().strip().split( "\n" )
    
    for row in rows:
        eval = Eval()
        if( not eval.parse( row ) ):
            print ("FATAL", "Syntax error in evals file '%s' row '%s'" %(fname, row) )
        
        eval.ok = True
        evals.append( eval )
    
    return( evals )

def save_evals_list( evals, fname ):
    fo = open( fname, "w" )
    for eval in filter( lambda e: e.ok, evals ):
        fo.write( "%s\n" %eval )
    fo.close()
    
def compute_evals_ranks( evals ):
    for (attr, reverse) in [("rmsd", False), ("pvalue", False), ("DI_ALL", False), ("INF_ALL", True), ("INF_WC", True), ("INF_NWC", True), ("INF_STACK", True), ("clashscore", False), ("mcq",False), ("gdt",False)]:
        evals_sorted = sorted( evals, key=attrgetter(attr), reverse=reverse )
        
        for (i, eval) in enumerate( evals_sorted ):
            eval.set_rank( attr, i+1 )

def sort_evals( evals, attr ):
    evals.sort( key=attrgetter(attr) )
