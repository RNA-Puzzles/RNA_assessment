#  pdb_utils.py
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
#  core functions for rna structure manipulation and comparison

import copy
import math
import os

from Bio.PDB import *

from .msgs import *
from .mcannotate import *
#'from .utils import *
from .extract import *

BIN_DIR=os.getcwd()
# from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
# Implements the Gauss error function.
#   erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z)
#
# fractional error in math formula less than 1.2 * 10 ^ -7.
# although subject to catastrophic cancellation when z in very close to 0
# from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
def erf(z):
		t = 1.0 / (1.0 + 0.5 * abs(z))
		# use Horner's method
		ans = 1 - t * math.exp( -z*z -  1.26551223 +
								  t * ( 1.00002368 +
								  t * ( 0.37409196 + 
								  t * ( 0.09678418 + 
								  t * (-0.18628806 + 
								  t * ( 0.27886807 + 
								  t * (-1.13520398 + 
								  t * ( 1.48851587 + 
								  t * (-0.82215223 + 
								  t * ( 0.17087277))))))))))
		if z >= 0.0:
				return ans
		else:
				return -ans

class PDBNormalizer:
	MAX_ERRORS = 5
	
	def __init__(self, fres_list, fatoms_list):
		self._load_res_list( fres_list )
		self._load_atom_list( fatoms_list )
		
	def parse( self, finput, foutput ):
		# state variables for the parse process
		self._in_model = False
		self._in_atom = False
		
		self._chain_found = False
		self._row_count = 0
		self._ok = True

		out_txt = ""
		
		fi = open( finput )
	
		for row in fi:
			self._row_count += 1
			
			row = row.strip()

			rec_name = row[:6]

			if( rec_name == "MODEL " ):
				row = self.parse_model( row )
				row = ""
			elif( rec_name == "ENDMDL" ):
				row = self.parse_endmdl( row )
				row = ""
			elif( rec_name[:3] == "TER" ):
				row = self.parse_ter( row )
			elif( rec_name in ("ATOM  ", "HETATM") ):
				row = self.parse_atom( row )
			else:
				continue
			
			if( row != "" ):
				out_txt += row + "\n"
		
		fi.close()
		
		if( self._in_atom ):
			out_txt += "TER\n"

		if( self._ok ):
			open( foutput, "w" ).write( out_txt )

		return( self._ok )
	
	def parse_model(self, row ):
		if( self._in_model ):
			self.show_err( "'ENDMDL' not found." )
		
		if( self._in_atom ):
			self.show_err( "Missing 'MODEL' before 'ATOM' declaration." )

		self._in_model = True
		return( row )

	def parse_endmdl(self, row):
		
		if( not self._in_model ):
			msgs.show( "Warning", "Missing 'MODEL' declaration." )
			#~ self.show_err( "Missing 'MODEL' declaration." )
		
		if( self._in_atom ):
			show( "Warning", "Missing 'TER' declaration." )
			#~ self.show_err( "Missing 'TER' declaration." )
		self._in_model = False
		self._in_atom = False
		return( "ENDMDL" )

	def parse_ter(self, row):
		result = ""
		
		if( self._in_atom ):
			result = "TER"
		
		self._in_atom = False
		return( result )

	def parse_atom(self, row):
		# get all the fields from the line
		serial = row[6:11]
		name = row[12:16].strip()
		altLoc = row[16]
		resName = row[17:20].strip()
		chainID = row[21]
		resSeq = int(row[22:26])
		iCode = row[26]
		x = row[30:38]
		y = row[38:46]
		z = row[46:54]
		occupancy = row[54:60]
		tempFactor = row[60:66]
		element = row[76:78]
		charge = row[78:80]

		# check residue name
		resName_norm = self._res_list.get( resName, None )
		if( resName_norm is None ):
			self.show_err( "Unknown residue name: '%s'." %resName )
			return ""
		elif( resName_norm == "-" ):
			return ""

		resName = resName_norm

		# check atom name
		name_norm = self._atom_list.get( name, None )
		if( name_norm is None ):
			self.show_err( "Unknown atom name: '%s' in residue'%s'" %(name, resName) )
			return ""
		elif( name_norm == "-" ):
			return ""
		
		name = name_norm.ljust(3)
		
		# check chainID
		if( chainID == " " ):
			if( self._chain_found ):
				self.show_err( "One of the chains is missing!" )
			else:
				chainID = "A"
		else:
			self._chain_found = True
		
		# check occupancy
		#if( occupancy == "" ):
		
		# This is the only that works with MolProbity
		occupancy = "  1.00"

		# check tempFactor
		if( tempFactor == "" ):
			tempFactor = "  0.00"
			
		# check element
		if( element == "" ):
			element = name[0]

		self._in_atom = True
		return "ATOM  %5s  %3s%s%3s %s%4d%s   %8s%8s%8s%6s%6s		  %2s%2s" %(serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge)

	def show_err( self, msg ):
		show( "ERROR", "Line %d: %s\n" %(self._row_count, msg) )
		self._ok = False

	def _load_res_list(self, fres_list):
		# read residues list
		pairs = map( lambda x: x.split(), filter( lambda row: not row.startswith( "#" ), open( fres_list ).read().strip().split( "\n" ) ) )

		self._res_list = {}		 
		for (name, nt) in pairs:
			self._res_list[name] = nt

	def _load_atom_list(self, fatoms_list):
		# read residues list
		pairs = map( lambda x: x.split(), filter( lambda row: not row.startswith( "#" ), open( fatoms_list ).read().strip().split( "\n" ) ) )

		self._atom_list = {}		 
		for (name, name_norm) in pairs:
			self._atom_list[name] = name_norm

#		
# get the sequence list from a pdb either raw or indexed
#
class Residue:
	def __init__(self, chain, pos, nt, res):
		self.chain = chain
		self.pos = pos
		self.nt = nt
		self.res = res
		
	def key(self):
		return "%s:%s" %(self.chain, self.pos)
	
	def __str__(self):
		return "%s:%s:%s > %s" %(self.chain, self.pos, self.nt, self.res)

class PDBStruct(object):
	def __init__(self):
		self._pdb_file = None
		self._struct = None
		self._res_list = []
		self._res_seq = []
		self._res_index = {}
		self._interactions = []
		#self._brackets = []
		#self._wcpairs = []
	
	def load(self, pdb_file, index_name=None ):
		self._pdb_file = pdb_file
		
		ok = self._load_struct()
		
		if( ok and not index_name is None ):
			ok = self._load_index( index_name )
		else:
			ok = self._load_index2()
			
		if( ok ):
			ok = self._load_annotations_3D()
		#~ print pdb_file,self._interactions,index_name
		#if( ok and not fbrackets is None ):
		#	ok = self._load_brackets( fbrackets )
		
		return( ok )
	
	def raw_sequence(self):
		seq = ""
		for ndx in self._res_seq:
			seq += self._res_list[ndx].nt
	
		return seq

	def res_sequence(self):
		result = []
		for ndx in self._res_seq:
			result.append( self._res_list[ndx].res )
	
		return result

	def get_interactions(self, type="ALL"):
		if( type == "ALL" ):
			# "ALL": returns all interactions
			return self._interactions
		elif( type in ("PAIR") ):
			# "PAIR": returns all pairs irrespective of their type
			return list(filter( lambda x: x[0] in ("PAIR_2D", "PAIR_3D"), self._interactions ))
		elif( type in ("PAIR_2D", "PAIR_3D", "STACK") ):
			# "PAIR_2D", "PAIR_3D", "STAK": returns the interactions of the specified type
			return list(filter( lambda x: x[0] == type, self._interactions ))
		else:
			show( "FATAL", "Wrong interaction type '%s' expected: 'ALL', 'PAIR', 'PAIR_2D', 'PAIR_3D' or 'STACK'" %type)
	
	# --- properties ---
	def struct_get(self):
		return self._struct
	
	def res_seq_get(self):
		return self._res_seq

	def res_list_get(self):
		return self._res_list
	
	def pdb_file_get(self):
		return self._pdb_file
	
	def rad_gir(self):
		rmean = Vector([0.0, 0.0, 0.0])
		count = 0
		for res in self._res_list:
			for a in res.res:
				rmean += a.get_vector()
				count += 1
		
		rmean = rmean/float(count)

		rsum = 0.0
		for res in self._res_list:
			for a in res.res:
				rsum += (a.get_vector() - rmean) * (a.get_vector() - rmean)
		
		return math.sqrt( rsum / count )
		
	
	#def brackets_get(self):
	#	return self._brackets
	
	struct = property( struct_get )
	res_seq = property( res_seq_get )
	res_list = property( res_list_get )
	pdb_file = property( pdb_file_get )
	#brackets = property( brackets_get )
	# ---
	
	def _load_struct(self):
		parser = PDBParser()
		self._struct = parser.get_structure( "struct", self._pdb_file )
		
		if( len(self._struct) > 1 ):
			show( "WARNING", "%d models found. Only the first will be used!" %(len(self._struct)) )

		self._res_list = []
		self._res_seq = []
		self._res_index = {}
		
		# gets only the first model
		model = self._struct[0]
		count = 0
		for chain in model.child_list:
			for res in chain.child_list:
				new_residue = Residue(chain.id, res.id[1], res.resname.strip(), res)
				
				self._res_list.append( new_residue )
				self._res_seq.append( count )
				self._res_index[new_residue.key()] = [count, None]
				
				count += 1

		return( True )
			
	def _load_index(self, index_name):
		self._res_seq = []
		entries = []
		for row in open( index_name ).read().split( "\n" ):
			row = row.strip()
			if( (not row.startswith( "#" )) and (row != "") ):
				entries.extend( map( lambda row: row.split( ":" ), row.split( "," ) ) )
		
		for entry in entries:
			if( len(entry) != 3 ):
				show( "ERROR", "Bad index entry: '%s'" %entry)
				return( False )
			
			chain = entry[0]
			pos = int(entry[1])
			count = int(entry[2])
		
			# get the index position
			ndx = self._get_index( chain, pos, 0 )
			
			if( ndx is None ):
				return( False )
			
			# get the positions
			for i in range( ndx, ndx + count ):
				if( i >= len(self._res_list) ):
					show( "ERROR", "Bad count %d in index entry: '%s'" %(count, entry) )
					return( False )
				
				if( self._res_list[i].chain != chain ):
					show( "ERROR", "Position %d in index entry: '%s' is outside the chain" %(i, entry) )
					return( False )
				self._res_seq.append( i )
				
				# update the index with the rank of the residue
				self._res_index[self._res_list[i].key()][1] = (len(self._res_seq) - 1)
		return( True )
	
	def _load_index2(self):
		self._res_seq = []
		for i in xrange( 0, len(self._res_list) ):
			self._res_seq.append( i )
			self._res_index[self._res_list[i].key()][1] = (len(self._res_seq) - 1)
		return True 

	def _load_annotations_3D(self):
		self._interactions = []
		mca = MCAnnotate()
		mca.load( self._pdb_file, os.path.dirname( self._pdb_file ) )
		#~ print mca.interactions
		for (type, chain_a, pos_a, nt_a, chain_b, pos_b, nt_b, extra1, extra2, extra3) in mca.interactions:
			# get the rank of the first position of the pair
			rank_a = self._get_index( chain_a, pos_a, 1 )
			rank_b = self._get_index( chain_b, pos_b, 1 )
			
			if( (rank_a is None) or (rank_b is None) ):
				continue
				#~ return False
			
			if( type == "STACK" ):
				extra = extra1
			else:
				extra = "%s%s" %(extra1, extra2)
			self._interactions.append( (type, min( rank_a, rank_b ), max( rank_a, rank_b ), extra ))
		 
	def _get_index(self, chain, pos, field):
		key = "%s:%s" %(chain, pos)
		data = self._res_index.get( key, None )[field]
		if data is None and field == 0:
			sys.stderr.write("ERROR	Bad index key: '%s'\n" %key)
		
		return data 
			
class PDBComparer:
	BACKBONE_ATOMS = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'","O4'", "O5'", "OP1", "OP2", "P"]
	HEAVY_ATOMS =	["C2", "C4", "C5", "C6", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]
	ALL_ATOMS = BACKBONE_ATOMS + HEAVY_ATOMS

	RMSDD_ATOMS =	["C4", "C8", "P", "C1'"]


	def __init__(self):
		pass
	
	def mcq(self, f1,f2):
		cmd='java -cp %s/mcq.ws.client-0.0.1-SNAPSHOT-jar-with-dependencies.jar pl.poznan.put.mcq.ws.client.Global -m %s -t %s >mcq.log'%(BIN_DIR,f1,f2)
		os.system(cmd)
		try:
			v=float(open('mcq.log').read().strip())
		except Exception:
			v=0
		return v
	
	def gdt(self, f1, f2):
		cmd='java -jar %s/gdt.jar %s %s >gdt.log'%(BIN_DIR,f2,f1)
		#~ print cmd
		os.system(cmd)
		try:
			x=open('gdt.log').read().strip().split('\n')[1].split(',')[-1]
			if x == 'NaN':return 0
			v=float(x)
		except Exception:
			v=0
		return v
	
	def rmsd( self, src_struct, trg_struct, fit_pdb=None ):
		# for each model in the reference
		src_residues = src_struct.res_sequence()
		trg_residues = trg_struct.res_sequence()
		
		atoms = self._get_atoms_struct( PDBComparer.ALL_ATOMS, src_residues, trg_residues )
		
		if( not atoms is None ):
			(src_atoms, trg_atoms) = atoms
		else:
			return None

		# compute the rmsd value and apply it to the target structure			
		sup = Superimposer()
		sup.set_atoms( src_atoms, trg_atoms )
		
		# we copy the fit_struct to leave the target struct unmodified for posterior processing
		fit_struct = copy.deepcopy( trg_struct.struct )
		sup.apply( fit_struct.get_atoms() )
				
		# save the fitted structure
		if( not fit_pdb is None ):
			io = PDBIO()
			io.set_structure( fit_struct )
			io.save( fit_pdb )

		return sup.rms
		
	
	# From Hajdin et al., RNA (7) 16, 2010 
	def pvalue( self, m, N, param ):
		if( param == "+" ):
			a = 5.1
			b = 15.8
		elif( param == "-" ):
			a = 6.4
			b = 12.7
		else:
			show( "FATAL", "Wrong p-value parameter '%s'. Expected '+' or '-'" %param )
			
		RMSD = a * (N ** 0.41) - b
		
		Z = (m - RMSD) / 1.8
		
		pv = (1.0 + erf( Z / (2**0.5) )) / 2.0 
		
		return( pv )

	def INF(self, src_struct, trg_struct, type):
		(P, TP, FP, FN) = (0, 0, 0, 0)
		
		for (stype, sb1, sb2, sextra) in src_struct.get_interactions( type ):
			P += 1
			found = False
			for (ttype, tb1, tb2, textra) in trg_struct.get_interactions( type ):
				if( (stype==ttype) and (sb1==tb1) and (sb2==tb2) and (sextra==textra) ):
					found = True
					break
			
			if( found ):
				#print "TP>", (stype, sb1, sb2, sextra)
				TP += 1
			else:
				#print "FN>", (stype, sb1, sb2, sextra)
				FN += 1

		for (ttype, tb1, tb2, textra) in trg_struct.get_interactions( type ):
			found = False
			for (stype, sb1, sb2, sextra) in src_struct.get_interactions( type ):
				if( (stype==ttype) and (sb1==tb1) and (sb2==tb2) and (sextra==textra) ):
					found = True
					break
			
			if( not found ):
				FP += 1
				#print "FP>", (ttype, tb1, tb2, textra)
		
		if( TP == 0 and (FP == 0 or FN == 0) ):
			INF = -1.0
		else:
			PPV = float(TP) / (float(TP) + float(FP))
			STY = float(TP) / (float(TP) + float(FN))
			INF = (PPV * STY) ** 0.5
		
		#print "##>", INF, P, TP, FP, FN
		
		return( INF )
	
	def DP(self, src_struct, trg_struct, template_txt, dname, dp_script):
		# prepare the config file
		txt = ""
		txt += "matrix=True\n"
		txt += "quiet_err = True\n"
		txt += "out_dir = '%s'\n" %dname
		txt += "ref_model = ('%s', 0)\n" %src_struct.pdb_file
		txt += "cmp_model = [('%s', 0)]\n" %trg_struct.pdb_file
		
		aligns = self._build_dp_alignments(src_struct, trg_struct)
		aligns_txt = []
		
		for align in aligns:
			aligns_txt.append( "('%s', %s, '%s', %s, %s)" %(align[0], align[1], align[2], align[3], align[4]) )
			
		txt += "aligns = [%s]\n" %(", ".join( aligns_txt) )
		txt += template_txt
		
		fname_cfg = "%s.cfg" %(trg_struct.pdb_file)
		fname_log = "%s.log" %(trg_struct.pdb_file) 
		open( fname_cfg, "w" ).write( txt )
		
		# runs the DP generator
		os.system( "python %s -c %s > %s" %(dp_script, fname_cfg, fname_log) )
	
	def VARNA(self, src_struct, trg_struct, algorithm="radiate"):
		edges = {"W":"wc", "S":"s", "H":"h"}
		
		data = {}
		data["sequenceDBN"] = src_struct.raw_sequence()
		data["structureDBN"] = "." * len(src_struct.raw_sequence())

		aux_bps = []
		
		for (stype, sb1, sb2, sextra) in src_struct.get_interactions( "PAIR" ):
			color = "#FF0000"
			for (ttype, tb1, tb2, textra) in trg_struct.get_interactions( "PAIR" ):
				if( (stype==ttype) and (sb1==tb1) and (sb2==tb2) and (sextra==textra) ):
					color = "#00FF00"
					break
			aux_bps.append( "(%d,%d):color=%s,edge5=%s,edge3=%s,stericity=%s" %(sb1+1, sb2+1, color, edges[sextra[0]], edges[sextra[1]], sextra[2:]) )
			
		data["auxBPs"] = ";".join( aux_bps )
		data["algorithm"] = algorithm
		
		return( data )
	
	def _get_atoms_residue( self, atom_list, src_res, trg_res ):
		src_atom_list = []
		trg_atom_list = []
		
		src_atom_list_tmp = list(filter( lambda a: a.get_name() in atom_list, src_res ))
		trg_atom_list_tmp = list(filter( lambda a: a.get_name() in atom_list, trg_res ))
		
		# for each atom in reference
		for src_atom in src_atom_list_tmp:
			found = False
			src_name = src_atom.get_full_id()[4][0]
			
			# search for an atom with the same name in the comparison
			for trg_atom in trg_atom_list_tmp:
				trg_name = trg_atom.get_full_id()[4][0]
	
				# if the atom was found keep it and jump to the next
				if( src_name == trg_name ):
					src_atom_list.append( src_atom )
					trg_atom_list.append( trg_atom )
					found = True
					break
			
			if( not found ):
				show( "WARNING", "Atom %s from residue %s not found in target atom list" %(src_name, src_res.id))
		return( src_atom_list, trg_atom_list )
	
	def _get_atoms_struct( self, atom_list, src_residues, trg_residues ):
		src_atoms = []
		trg_atoms = []
	
		if( len(src_residues) != len(trg_residues) ):
			show( "ERROR", "Different number of residues!" )
			return( None )
	
		for (src_res, trg_res) in zip(src_residues, trg_residues):
			(sa, ta) = self._get_atoms_residue( atom_list, src_res, trg_res )
			
			src_atoms.extend( sa )
			trg_atoms.extend( ta )
		#'print('%d %d'%(len(src_atoms),len(trg_atoms)))
		return( src_atoms, trg_atoms )
	
	def _build_dp_alignments(self, src_struct, trg_struct):
		aligns = []
		
		(schain, tchain) = ("", "")
		(spos, tpos) = (-1, -1)
		count = 0
		item = None
		
		for (i, j) in zip(src_struct.res_seq, trg_struct.res_seq):
			(sres, tres) = (src_struct.res_list[i], trg_struct.res_list[j])

			# if the numbering or the chain change			
			if( (sres.pos != (spos+1)) or (tres.pos != (tpos+1)) or (sres.chain != schain) or (tres.chain != tchain) ):
				if( count > 0 ):
					item[4] = count
					aligns.append( item )
				
				(schain, tchain) = (sres.chain, tres.chain)
				item = [sres.chain, sres.pos, tres.chain, tres.pos, None]
				count = 1
			else:
				count += 1
			
			(spos, tpos) = (sres.pos, tres.pos)

		# adds the last alignment
		if( count > 0 ):			
			item[4] = count
			aligns.append( item )
		
		return( aligns )
