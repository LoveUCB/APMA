import requests
from Bio import SeqIO

def extract_sequence_from_pdb(pdb_file):
    sequences = []
    with open(pdb_file, "r") as handle:
        for record in SeqIO.parse(handle, "pdb-seqres"):
            sequences.append(str(record.seq))
    return sequences

pdb_sequences = extract_sequence_from_pdb('data/alphafoldpten.pdb')
protein_sequence = ''.join(pdb_sequences)

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def blast_search(sequence, output_file, blast_program="blastp", database="nr", evalue=1e-50, num_alignments=5000): # 搜索多少个，然后是不是使用人的数据库
    # Perform BLAST search
    result_handle = NCBIWWW.qblast(blast_program, database, sequence, expect=evalue, hitlist_size=num_alignments,alignments=5000)
    
    # Parse BLAST result
    blast_records = NCBIXML.parse(result_handle)
    with open(output_file, "w") as out_fasta:
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    out_fasta.write(">{} {}\n".format(alignment.title.split()[1], alignment.title))
                    out_fasta.write("{}\n".format(alignment.hsps[0].sbjct))
    
    print("BLAST search complete. Results saved to", output_file)

from Bio.Align.Applications import ClustalOmegaCommandline

def run_clustal(input_fasta, output_fasta):# 处理blast结果，生成的序列需要全部都是独立的
    # 设置Clustal Omega的命令行参数
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta, outfile=output_fasta, verbose=True, auto=True,dealign=True)
    
    # 执行Clustal Omega
    stdout, stderr = clustalomega_cline()




# Example usage
if __name__ == "__main__":
    # Specify your sequence
    sequence = protein_sequence
    f = open("data/sequence.txt","w")
    f.write(sequence)
    f.close()
    # Perform BLAST search and save results to a file
    output_file = "/home/wangjingran/APMA/data/blast_results.fasta"
    blast_search(sequence, output_file)
        # 输入的FASTA文件，这里假设你已经有了一些同源序列的FASTA文件
    input_fasta = "data/blast_resultss.fasta"
    # 输出的FASTA文件，用于保存比对结果
    output_fasta = "data/msaa.fasta"
    # 运行多序列比对
    run_clustal(input_fasta, output_fasta)
