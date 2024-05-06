from Bio import SeqIO
import time

def extract_sequence_from_pdb(pdb_file):
    sequences = []
    with open(pdb_file, "r") as handle:
        for record in SeqIO.parse(handle, "pdb-seqres"):
            sequences.append(str(record.seq))
    return sequences


from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def blast_search(sequence, output_file, blast_program="blastp", database="swissprot", evalue=1e-15, num_alignments=1000): # 搜索多少个，然后是不是使用人的数据库
    # Perform BLAST search
    print("BLAST search started...")
    result_handle = NCBIWWW.qblast(blast_program, database, sequence, expect=evalue, hitlist_size=num_alignments,alignments=1000)
    # Parse BLAST result
    blast_records = NCBIXML.parse(result_handle)
    with open(output_file, "w") as out_fasta:
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    out_fasta.write(">{} E-Value:{}\n".format(alignment.title.split("|")[1], hsp.expect))
                    out_fasta.write("{}\n".format(alignment.hsps[0].sbjct))
    
    print("BLAST search complete. Results saved to", output_file)

from Bio.Align.Applications import ClustalOmegaCommandline

def run_clustal(input_fasta, output_fasta):# 处理blast结果，生成的序列需要全部都是独立的
    # 设置Clustal Omega的命令行参数
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta, outfile=output_fasta, verbose=True, auto=True,dealign=True)
    
    # 执行Clustal Omega
    stdout, stderr = clustalomega_cline()


import subprocess
def run_blast(input, output, dbname, evalue, alignments):
    command = f'''blastp -query {input}
    -out {output}
    -db {dbname}
    -outfmt 5 -evalue {evalue}
    -num_threads 2 -num_alignments {alignments} 
    -max_target_seqs 200'''
    subprocess.run(command, shell=True)

if __name__ == "__main__":
    sequence = extract_sequence_from_pdb("/Users/wangjingran/Desktop/APMA Server/files_PIK3CA/PIK3CA.pdb")
    sequence = ''.join(sequence)
    print(sequence)
    max_try_for_blast = 3
    current_try_for_blast = 0
    while current_try_for_blast < max_try_for_blast:
        current_try_for_blast += 1
        try:
            output_file = "/Users/wangjingran/Desktop/APMA/data/PIK3CA.fasta"
            print(f"BLAST Search Started {current_try_for_blast} time")
            blast_search(sequence, output_file)
            print(f"BLAST Search success")
            break
        except Exception as e:
            print(f"Blast search failed {current_try_for_blast} times, {3 - current_try_for_blast} remaining")
            print(f"Error: {e}")
            time.sleep(30)
    else:
        print("Error: BLAST search failed after multiple tries.")

