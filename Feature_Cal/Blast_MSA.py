import warnings
warnings.filterwarnings('ignore')
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import PDB
import time
import os


def extract_sequence_from_pdb(pdb_file):
    """
    Extracts the sequence of residues from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file to be parsed.

    Returns:
        list: A list of strings representing the sequences of residues for each chain in the PDB file. 
              Each residue is represented by its one-letter code.
    
    Notes:
        The function uses a dictionary to map three-letter residue names to their corresponding one-letter codes.
        It processes the PDB file using the PDBParser from the Biopython library to extract the residue sequence.
        
    Example:
        >>> sequences = extract_sequence_from_pdb('example.pdb')
        >>> print(sequences)
        ['C', 'A', 'D', 'G', ...]  # The output will be a list of sequences for each chain
    """

    # Dictionary to map three-letter residue names to one-letter codes
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    # Initialize the PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)    

    # List to store sequences of residues for each chain
    sequences = []

    # Iterate over each model in the structure
    for model in structure:
        # Iterate over each chain in the model
        for chain in model:
            seq = []
            # Iterate over each residue in the chain
            for residue in chain:
                # Append the one-letter code of the residue to the sequence
                seq.append(d3to1[residue.resname])
            # Add the sequence of this chain to the list of sequences
            sequences.append(seq)

    return sequences


# def extract_sequence_from_pdb(pdb_file):
#     sequences = []
#     with open(pdb_file, "r") as handle:
#         for record in SeqIO.parse(handle, "pdb-seqres"):
#             sequences.append(str(record.seq))
#     return sequences


def blast_search(sequence, fasta_file, output_file, blast_program="blastp", evalue=1e-50, num_alignments=3000):
    """
    Performs a BLAST search against a sequence database using a given query sequence and saves the results to an output file.

    Args:
        sequence (str): The query sequence to be searched.
        fasta_file (str): Path to the FASTA file containing the sequence database for BLAST search.
        output_file (str): Path to the output file where the BLAST results will be saved.
        blast_program (str, optional): The BLAST program to use (default is "blastp"). Other options include "blastn", "blastx", etc.
        evalue (float, optional): The E-value threshold for reporting matches (default is 1e-50).
        num_alignments (int, optional): The maximum number of alignments to report (default is 3000).

    Returns:
        None

    Notes:
        - The function creates a temporary query file for the sequence and performs the BLAST search using the specified BLAST program.
        - Results are saved in XML format and parsed to extract relevant information.
        - Temporary files are cleaned up after processing.
        - The function prints status messages and any errors encountered during the BLAST search.

    Example:
        >>> blast_search("MKTIIALSYIFCLVFA", "database.fasta", "results.fasta")
        # BLAST search will be performed and results saved to 'results.fasta'
    """

    # Create a temporary input file for the query sequence
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
        num_threads=4
    )
    
    # Execute the BLAST command
    print("[INFO]...Searching...")
    stdout, stderr = blast_cline()
    print(stdout)
    print(stderr)
    
    # Parse the BLAST results
    with open("temp_blast.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        with open(output_file, "w") as out_fasta:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        out_fasta.write(">{} E-Value:{}\n".format(alignment.title.split("|")[1], hsp.expect))
                        out_fasta.write("{}\n".format(hsp.sbjct))
    
    # Cleanup temporary files
    os.remove(query_file)
    os.remove("temp_blast.xml")
    
    print("BLAST search complete. Results saved to", output_file)



def run_clustal(input_fasta, output_fasta):
    """
    Runs Clustal Omega to align sequences from an input FASTA file and save the aligned sequences to an output FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file containing sequences to be aligned.
        output_fasta (str): Path to the output FASTA file where the aligned sequences will be saved.

    Returns:
        None

    Notes:
        The function uses the Clustal Omega command-line tool to perform multiple sequence alignment.
        The following command-line options are set:
        - `verbose=True`: Provides detailed output.
        - `auto=True`: Automatically determines the alignment strategy.
        - `dealign=True`: Removes gaps from the alignment.

    Example:
        >>> run_clustal('input_sequences.fasta', 'aligned_sequences.fasta')
        # The aligned sequences will be saved to 'aligned_sequences.fasta'
    """
    
    # Set up the Clustal Omega command-line parameters
    clustalomega_cline = ClustalOmegaCommandline(
        infile=input_fasta,
        outfile=output_fasta,
        verbose=True,
        auto=True,
        dealign=True
    )
    
    # Execute Clustal Omega
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
