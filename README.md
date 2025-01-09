![APMA](Figure/LOGO.png)
# protPheMut ---- ML & Graph based Tool
## Introduction
![APMA](Figure/graph_abstract.png)
- **FROM MUTATIONS TO PHENOTYPES**
- The deePheMut Server is at: **[protPheMut Server](http://106.54.2.54/protPheMut)**
- **deePheMut** is an online server tool for distinguishing and predicting inseparable single gene missense mutation multiphenotypic diseases. We will automatically calculate the phenotypic mutation characteristics and use the machine learning framework for model training and interpretation. The final results will be presented in a visual result. We will provide biological indicators that distinguish different disease phenotypes caused by single gene missense mutations and mutation score results to measure the biological significance of the current mutation in this phenotype.
- More information can be found on **[protPheMut Server](http://106.54.2.54/protPheMut)**

## Data
- Protein structure data can be downloaded from the [AlphaFold Database](https://alphafold.ebi.ac.uk) or using your structure file
- Mutation data can be downloaded from public databases: [gnomAD](https://gnomad.broadinstitute.org) and [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) etc.

## Links
- ProDy: http://prody.csb.pitt.edu
- FoldX: https://foldxsuite.crg.eu
- NACEN: http://sysbio.suda.edu.cn/NACEN
- Clustal Omega: http://www.clustal.org/omega
- Rate4site: https://www.tau.ac.il/~itaymay/cp/rate4site.html
- Blast: https://blast.ncbi.nlm.nih.gov
- SHAP: https://shap.readthedocs.io/en/latest
- Uniref50: https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50

## Features
1. Entropy: measures the mutation frequency on the position
2. Coevolution: coevolution on the position
3. Conservation: conservation score based on rate4site
4. ddG: energy change based on FoldX
5. RASA: single residue's relative exposed area
6. Polarity: mutation's global effect on polarity sum(MT_Polarity - WT_Polarity)
7. Hydrophobicity: mutation's global effect on hydrophobicity sum(MT_Hydrophobicity - WT_Hydrophobicity)
8. Betweenness: mutations' combined global effect on betweenness on a residue rowmeans(matrix(MT_betweenness) - WT_betweenness)[position]
9. Closeness: mutations' combined global effect on closeness on a residue rowmeans(matrix(MT_closeness) - WT_closeness)[position]
10. Eigenvector: mutations' combined global effect on eigenvector on a residue rowmeans(matrix(MT_eigenvector) - WT_eigenvector)[position]
11. Effectiveness: dynamic network features based on ProDy
12. Sensitivity: dynamic network features based on ProDy
13. DFI: dynamic network features based on ProDy
14. MSF: dynamic network features based on ProDy
15. Sensitivity: dynamic network features based on ProDy


# Install
- The deePheMut is linux-based tool, do not use Windows or MacOS
## Installation of Blast tool
1. In order to do blastp, please download ncbi-blast+
```sh
sudo apt get install ncbi-blast+
```
2. Download your local blast database
> In my research, I used uniref50, be sure your server have enough space.
```sh
wget -h https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
```
```sh
unzip uniref50.fasta.gz
```
3. Construct local blast database
```sh
makeblastdb -in uniref50.fasta -dbtype prot
```
## Installation of Clustal Omega and Rate4site
1. Download Clustal Omega
```sh
sudo apt get install clustalo
```
2. Download Rate4site
```sh
sudo apt get install rate4site
```
## Installation of FoldX
1. The FoldX is in https://foldxsuite.crg.eu
2. Please move all the files into deePheMut/mutation folder

## Installation of deePheMut
1. To get the tool, run the following code
```
git clone https://github.com/Spencer-JRWang/deePheMut
```
2. To install the python dependency, run the following code
```
pip install .
```
## Installation of R dependency
1. To install the R base, run the following code
```sh
sudo apt-get install r-base
```
2. Open R console
```sh
sudo R
```
3. Install bio3d
```R
install.packages("bio3d")
```
4. Install igraph
```R
install.pacakges("igraph")
```
5. Install NACEN
```R
install.packages("/your/route/to/NACEN", repos = NULL, type = "source")
```
> You can get NACEN source code at NACEN website: http://sysbio.suda.edu.cn/NACEN

> Or you can get NACEN from [NACEN package](data/NACEN_0.1.0.tar.gz)


# Message

> 📧: spencer-jrwang@foxmail.com
>
> Department of Bioinformatics, Medical School of Soochow University
>
> protPheMut is free for everyone to use, if you have used our tools in your research, please cite: 
> Jingran Wang, Miao Yang, Chang Zong, Gennady Verkhivker, Fei Xiao, Guang Hu. protPheMut: An Interpretable Machine Learning Tool for Classification of Cancer and Neurodevelopmental Disorders in Human Missense Variants. bioRxiv 2025.01.06.631365; doi: https://doi.org/10.1101/2025.01.06.631365

