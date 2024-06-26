import warnings
warnings.filterwarnings('ignore')
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
import time
import os


def extract_sequence_from_pdb(pdb_file):
    sequences = []
    with open(pdb_file, "r") as handle:
        for record in SeqIO.parse(handle, "pdb-seqres"):
            sequences.append(str(record.seq))
    return sequences


def blast_search(sequence, fasta_file, output_file, blast_program="blastp", evalue=1e-50, num_alignments=3000):
    # Create a temporary input file for the query sequence
    # os.environ['PATH'] = "/Users/wangjingran/Downloads/ncbi-blast-2.15.0+/bin:" + os.environ['PATH']
    query_file = "/home/wangjingran/APMA/data/temp_query.fasta"
    with open(query_file, "w") as seq_file:
        seq_file.write(">Query\n")
        seq_file.write(sequence)
    

    
    # Define the BLAST command
    blast_cline = NcbiblastpCommandline(
        cmd=blast_program,
        query=query_file,
        db=fasta_file,
        evalue=evalue,
        outfmt=5,
        out="temp_blast.xml",
        max_target_seqs=num_alignments,
        num_threads = 4
    )
    
    # Execute the BLAST command
    print("[INFO]...Searching...")
    stdout, stderr = blast_cline()
    print(stdout)
    print(stderr)
    
    # Parse the BLAST results
    with open("/home/wangjingran/APMA/data/temp_blast.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        with open(output_file, "w") as out_fasta:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        out_fasta.write(">{} E-Value:{}\n".format(alignment.title.split("|")[1], hsp.expect))
                        out_fasta.write("{}\n".format(alignment.hsps[0].sbjct))
    
    # Cleanup temporary files
    os.remove(query_file)
    os.remove("/home/wangjingran/APMA/data/temp_blast.xml")
    
    print("BLAST search complete. Results saved to", output_file)


def run_clustal(input_fasta, output_fasta):# 处理blast结果，生成的序列需要全部都是独立的
    # 设置Clustal Omega的命令行参数
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta, outfile=output_fasta, verbose=True, auto=True,dealign=True)
    
    # 执行Clustal Omega
    stdout, stderr = clustalomega_cline()


if __name__ == "__main__":
    sequence = extract_sequence_from_pdb("/Users/wangjingran/Desktop/SpencerW-APMA/APMA Connect/files_ALPL/alpl.pdb")
    sequence = ''.join(sequence)
    print(sequence)
    fasta_file = '/Users/wangjingran/Desktop/prdatabase/uniref50.fasta'
    max_try_for_blast = 3
    current_try_for_blast = 0
    while current_try_for_blast < max_try_for_blast:
        current_try_for_blast += 1
        try:
            output_file = "/Users/wangjingran/Desktop/APMA/data/alpl.fasta"
            print(f"...BLAST Search Started {current_try_for_blast} time...")
            blast_search(sequence, fasta_file ,output_file)
            print(f"BLAST Search success")
            break
        except Exception as e:
            print(f"Blast search failed {current_try_for_blast} times, {3 - current_try_for_blast} remaining")
            print(f"Error: {e}")
    else:
        print("Error: BLAST search failed after multiple tries.")
