from Bio import *
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser
def Cal_Sasa(filename):
    '''
    Function to calculate SASA score for each Amino Acid residue in a PDB file
    '''
    p = PDBParser(QUIET=1)
    struct = p.get_structure("protein","../data/alphafoldpten.pdb")
    sr = ShrakeRupley()
    sr.compute(struct[0], level="R")
    sasa_list = []
    # print(len(struct))
    for chain in struct[0]:
        for residue in chain:
            sasa_list.append((residue.get_resname(), residue.sasa))
    return sasa_list

if __name__ == "__main__":
    RASA = Cal_Sasa("../data/alphafoldpten.pdb")
    print(RASA)