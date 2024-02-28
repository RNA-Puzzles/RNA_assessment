#!/usr/bin/env python

#
# Fits two or more molecules
#
import sys
import os

from Bio.PDB import *

class MySelect(Select):
    def config( self, res_list ):
        self.res_list = []

        for res_data in res_list:
            chain = res_data[0]
            res_id = int(res_data[1])
            count = int(res_data[2])

            for i in range( 0, count ):
                self.res_list.append( "%s|%s" %(chain, str(res_id+i)) )

    def accept_residue(self, residue):
        key = "%s|%s" %(residue.parent.id.strip(), residue.get_id()[1])

        save_it = key in self.res_list

        #if( save_it ):
        #    print key

        return( save_it )


def WritePDB( struct, select_class, file ):
    io = PDBIO()
    io.set_structure( struct )
    io.save( file, select_class )

        
def parse_res_list( s ):
    res_list = []
    pieces = s.split( "," )
    
    for piece in pieces:
        data = piece.split( ":" )

        if( len(data) != 3 ):
            print("Wrong data: %s!" %piece) 
            quit()

        res_list.append( data )

    return( res_list )

def extract_PDB(p1,p2,p3):
	
	parser = PDBParser()
	# Open and parse the structure (the first parameter is arbitrary)
	sinput = parser.get_structure( "SI", p1 )

	# prepares the select class
	select_class =  MySelect()
	res_list = parse_res_list( p2 )
	select_class.config( res_list )

	WritePDB( sinput, select_class, p3 )
#
# Main
#
# Create the PDB file parser
if __name__ == '__main__':
    model_list = ['PZ22_Xiao_1', 'PZ22_Dokholyan_5', 'PZ22_Xiao_4', 'PZ22_Dokholyan_4', 'PZ22_SimRNA_3',
                  'PZ22_Bujnicki_5', 'PZ22_SimRNA_4', 'PZ22_Bujnicki_4', 'PZ22_SimRNA_2', 'PZ22_SimRNA_1',
                  'PZ22_Bujnicki_3', 'PZ22_Dokholyan_1', 'PZ22_Dokholyan_2']
    path_dir=r'/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/Puzzles'
    coords = "A:1:13,A:15:67,B:69:13"
    pdb_files=[os.path.join(path_dir,f"{i}_splitchainA.pdb") for i in model_list]
    for i in pdb_files:

        new_pdb = os.path.basename(i).split('.pdb')[0] + '_deleteU.pdb'
        extract_PDB(i, coords, os.path.join(path_dir,'pz22deleteu',new_pdb))

    # original_pdb = "/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/solution/PZ22_extract_1.pdb"
    # coords_new = "A:69:13,B:1:13,B:15:67"
    # new_pdb = "/Users/bf/Documents/1doctor/RoundV_all/aFigure_preprocess/lddt/step4_all/solution/PZ22_extract_1_deleteC.pdb"
    # extract_PDB(original_pdb, coords_new, new_pdb)

	

