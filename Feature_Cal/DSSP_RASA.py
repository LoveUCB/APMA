from Bio.PDB import PDBParser
from Bio.PDB import DSSP

def DSSP_RASA(position, filename):
    print("RASA Calculating...", end=" ")
    protein_rasa = []
    p = PDBParser()
    structure = p.get_structure("protein", filename)
    model = structure[0]
    dssp = DSSP(model,filename)
    for i in range(len(dssp.keys())):
        a_key = list(dssp.keys())[i]
        protein_rasa.append(dssp[a_key][3])
    new_rasa = []
    for pos in position:
        selected_rows = protein_rasa[pos - 1]
        new_rasa.append(selected_rows * 100)
    print("Done")
    return new_rasa

def all_DSSP_RASA(filename):
    print("RASA Calculating...", end=" ")
    protein_rasa = []
    p = PDBParser()
    structure = p.get_structure("protein", filename)
    model = structure[0]
    dssp = DSSP(model,filename)
    for i in range(len(dssp.keys())):
        a_key = list(dssp.keys())[i]
        protein_rasa.append(dssp[a_key][3])
    print("Done")
    return protein_rasa

if __name__ == "__main__":
    protein_rasa = all_DSSP_RASA("/home/wangjingran/APMA/data/alphafoldpten.pdb")
    print(protein_rasa)