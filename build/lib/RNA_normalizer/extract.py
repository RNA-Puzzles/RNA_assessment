#
# Fits two or more molecules
#
import sys

from Bio.PDB import *

class MySelect(Select):
    def config( self, res_list ):
        self.res_list = []

        for res_data in res_list:
            chain = res_data[0]
            res_id = int(res_data[1])
            count = int(res_data[2])

            for i in xrange( 0, count ):
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
            print ("Wrong data: %s!" %piece)
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
	if( len(sys.argv) < 4 ):
		print("\n\n%s" %("- " * 40))
		print("extract.py - extracts the selected residues from the input PDB file")
		print("%s\n" %("- " * 40))
		print("Usage:")
		print("$ python extract.py <input model> <residue list> <output model>\n")
		print("<input model> - Initial model in PDB format.")
		print("<residue list> - Residues to extract")
		print("<output model> - Name of the new file\n")
		print("Residue lists should be in the following format:\n'chain:res_id:count,...,chain:res_id:count'\n")
		print("Examples:")
		print("$ python extract.py INPUT.pdb A:1:10,A:21:5,B:2:9 OUTPUT.pdb\b\b")
		quit()

	extract_PDB(sys.argv[1],sys.argv[2],sys.argv[3])
	

