from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
def run_blast(seq):
    result_handle = NCBIWWW.qblast("blastp", "nr", seq)
    blast_result = open("my_blast.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()
    E_VALUE_THRESH = 1e-50
    for record in NCBIXML.parse(open("my_blast.xml")):
        if record.alignments : #skip queries with no matches
            print ("QUERY: %s" % record.query[:60])
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print("MATCH: %s " % align.title[:60])
                        print(hsp.expect)

def extract_sequence_from_pdb(pdb_file):
    sequences = []
    with open(pdb_file, "r") as handle:
        for record in SeqIO.parse(handle, "pdb-seqres"):
            sequences.append(str(record.seq))
    return sequences

if __name__ == "__main__":
    sequence = extract_sequence_from_pdb("/Users/wangjingran/Desktop/APMA Server/files_ALPL/alpl.pdb")
    run_blast(sequence)
