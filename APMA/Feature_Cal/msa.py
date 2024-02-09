from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline

def pdb_to_msa_fasta(pdb_file, output_fasta):
    # Step 1: 使用PDB文件中的蛋白质序列作为查询序列
    # 解析PDB文件并提取蛋白质序列
    print("FASTA searching...",end=' ')
    from Bio.PDB import PDBParser
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]
    chain = model
    query_sequence = ''.join([residue.resname for residue in chain.get_residues() if residue.get_id()[0] == ' '])

    # Step 2: 使用BLAST进行远程同源序列搜索
    # 提交BLAST搜索请求
    result_handle = NCBIWWW.qblast('blastp', 'nr', query_sequence)

    # 解析BLAST搜索结果
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)

    # 获取BLAST搜索结果中的同源序列
    homologous_sequences = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            homologous_sequences.append(hsp.sbjct)
    print("Done")
    # Step 3: 使用Clustal Omega进行多序列比对
    # 将查询序列和同源序列保存为FASTA文件
    with open('query.fasta', 'w') as f:
        f.write(f'>query\n{query_sequence}\n')
    with open('homologous.fasta', 'w') as f:
        for i, seq in enumerate(homologous_sequences):
            f.write(f'>homologous_{i+1}\n{seq}\n')

    # 运行Clustal Omega进行多序列比对
    clustalomega_exe = '/path/to/clustalo'  # 替换成Clustal Omega的可执行文件路径
    clustalomega_cline = ClustalOmegaCommandline(clustalomega_exe, infile='homologous.fasta', outfile='aligned.fasta', verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()

    # Step 4: 将MSA保存为FASTA文件
    # 读取多序列比对结果并保存为FASTA文件
    with open('aligned.fasta', 'r') as f:
        with open(output_fasta, 'w') as out_f:
            for record in SeqIO.parse(f, 'fasta'):
                SeqIO.write(record, out_f, 'fasta')



pdb_file = '../data/alphafoldpten.pdb'
output_fasta = 'output.fasta'
pdb_to_msa_fasta(pdb_file, output_fasta)

