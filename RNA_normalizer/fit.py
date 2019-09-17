#
# Fits two or more molecules
#
import os

from Bio.Clustalw import *
from Bio.PDB import *

BACKBONE = ["C1'", "C1*", "C2'", "C2*", "C3'", "C3*", "C4'", "C4*", "C5'", "C5*", "O2'", "O2*", "O3'", "O3*", "O4'", "O4*", "O5'", "O5*", "P"]
FULL_ATOMS = ["C1'", "C1*", "C2'", "C2", "C2*", "C3'", "C3*", "C4'", "C4", "C4*", "C5'", "C5", "C5*", "C6", "C8", "N1", "N2", "N3", "N7", "N9", "O2'", "O2*", "O3'", "O3*", "O4'", "O4*", "O5'", "O5*", "O6", "P"]

#ATOM_LIST = BACKBONE
ATOM_LIST = FULL_ATOMS

def WritePDB( struct, file ):
    io = PDBIO()
    io.set_structure( struct )
    io.save( file )

def ResiduesFromModel( model, res_list ):
    residues = []

    for res_data in res_list:
        chain = res_data[0]
        res_id = int(res_data[1])
        count = int(res_data[2])

        if( res_id < 0 ):
            residues.extend( [None] * count )
        else:
            res_id = int(res_id)
            
            all_residues = [r for r in model[chain]]
            sub_residues = []
    
            start = False
            for res in all_residues:
                if( (res.get_id()[1] == res_id) ):
                    start = True
    
                if( start and (len(sub_residues) < count) ):
                    sub_residues.append( res )
    
            if( len(sub_residues) < count ):
                print("!!! Error: less than '%d' residues in chain '%s' starting from nt '%d'" %(count, chain, res_id))
                quit()
    
            residues.extend( sub_residues )

    return( residues )

def GetAtomsFromResiduesAux( ref_res, cmp_res ):
    global ATOM_LIST
    
    rr_name = ref_res.get_resname().strip()
    cr_name = cmp_res.get_resname().strip()

    ref_atom_list = []
    cmp_atom_list = []
    count = 0

    ref_atom_list_tmp = [a for a in ref_res if (a.get_name() in ATOM_LIST)]
    cmp_atom_list_tmp = [a for a in cmp_res if (a.get_name() in ATOM_LIST)]

    # for each atom in reference
    found = False
    for ra in ref_atom_list_tmp:
        ra_name = ra.get_full_id()[4][0].replace( "*", "'" )
        
        # search for an atom with the same name in the comparison
        for ca in cmp_atom_list_tmp:
            ca_name = ca.get_full_id()[4][0].replace( "*", "'" )

            # if the atom was found keep it and jump to the next
            if( (ra_name == ca_name) and ((rr_name == cr_name) or (ra_name in ATOM_LIST)) ):
                ref_atom_list.append( ra )
                cmp_atom_list.append( ca )
                break

    return( ref_atom_list, cmp_atom_list, count )

def GetAtomsFromResidues( ref_residues, cmp_residues ):
    ref_atoms = []
    cmp_atoms = []
    count = 0
    i = 0
    j = 0

    if( len(ref_residues) != len(cmp_residues) ):
        print("!! Different number of residues!")

    for i in xrange( 0, min(len(ref_residues), len(cmp_residues)) ):
        rr = ref_residues[i]
        cr = cmp_residues[i]
        
        if( (rr is None) or (cr is None) ):
            continue

        (r_a, c_a, n) = GetAtomsFromResiduesAux( rr, cr )

        ref_atoms.extend( r_a )
        cmp_atoms.extend( c_a )

    return( ref_atoms, cmp_atoms )

def Fitter( ref_struct, cmp_struct, ref_res_list, cmp_res_list ):
    min_sup = None
    min_ref_model = None
    min_cmp_model = None
    
    # for each model in the reference
    ref_residues = ResiduesFromModel( ref_struct[0], ref_res_list )
    cmp_residues = ResiduesFromModel( cmp_struct[0], cmp_res_list )
    
    (ref_atoms, cmp_atoms) = GetAtomsFromResidues( ref_residues, cmp_residues )
        
    sup = Superimposer()
    sup.set_atoms( ref_atoms, cmp_atoms )
    sup.apply( cmp_struct.get_atoms() )
    
    #for (a1, a2) in zip(ref_atoms,cmp_atoms):
    #    print a1.parent.parent.id, a1.parent.id[1], a1.parent.resname, a1.id, a2.parent.parent.id, a2.parent.id[1], a2.parent.resname, a2.id, a1-a2

    return( cmp_struct, sup.rms )
        
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

def go_fit( pdb_ref, pdb_cmp, res_ref, res_cmp ):
    # Open and parse the structure (the first parameter is arbitrary)
    parser = PDBParser()
    
    s1 = parser.get_structure( "S1", pdb_ref )
    s2 = parser.get_structure( "S2", pdb_cmp )
    res_list_1 = parse_res_list( res_ref )
    res_list_2 = parse_res_list( res_cmp )
    
    # Fits s2 into s1
    return(  Fitter( s1, s2, res_list_1, res_list_2 ) )
    
#
# Main
#
if __name__ == '__main__':
    if( len(sys.argv) < 6 ):
        print("""\n\n%s
fit.py - fits to PDB files minimizing the RMSD of selected residues
%s\n
Usage:
$ python fit.py <ref. model> <cmp. model> <ref. res. list> <cmp. res. list> <out file>\n
<ref. model> - Reference model in PDB format.
<cmp. model> - Comparing model in PDB format.
<ref. res. list> - Residues in the ref. model.
<cmp. res. list> - Residues in the comp. model.
<out file> - Name of the fitted model file\n
Residue lists should be in the following format:\n'chain:res_id:count,...,chain:res_id:count'\n
Examples:
$ python fit.py AAAA.pdb BBBB.pdb A:1:10,A:21:5,B:2:9 X:1:10 result.pdb\b\b"""%("- " * 40, "- " * 40))
        
        quit()
    
    (sf, rmsd) = go_fit( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )
    
    if( sys.argv[5] != "-" ):
        WritePDB( sf, sys.argv[5] )
    
    print("RMSD: %s vs %s = %f" %(os.path.basename( sys.argv[1]), os.path.basename( sys.argv[2]), rmsd ))
