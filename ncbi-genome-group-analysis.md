## Analysis of NCBI's *Staphylococcus aureus* Genome Groups
We downloaded 25 representative genomes from [NCBI's *S. aureus* Genome Group](http://www.ncbi.nlm.nih.gov/genome/genomegroups/154?). 
The multi-locus sequence type was determined for each genome and varying coverages of FASTQ sequences were simulated using ART.

### List of genomes and corresponding MLST
| *S. aureus Genome*                                            | Sequence Type |
|:---------------------------------------------------------------|---------------:|
| GCF_000011265.1_ASM1126v1_genomic.fna                         | 1             |
| GCF_000009645.1_ASM964v1_genomic.fna                          | 5             |
| GCF_000013425.1_ASM1342v1_genomic.fna                         | 8             |
| GCF_000013465.1_ASM1346v1_genomic.fna                         | 8             |
| GCF_000215425.1_ASM21542v2_genomic.fna                        | 22            |
| GCF_000205385.1_ASM20538v2_genomic.fna                        | 25            |
| GCF_000011505.1_ASM1150v1_genomic.fna                         | 36            |
| GCF_000242475.1_ASM24247v2_genomic.fna                        | 45            |
| GCF_000237125.1_ASM23712v1_genomic.fna                        | 59            |
| GCA_001469105.1_ASM146910v1_genomic.fna                       | 126           |
| GCF_000210315.1_ASM21031v1_genomic.fna                        | 133           |
| GCF_000162635.1_ASM16263v1_genomic.fna                        | 145           |
| GCF_000009005.1_ASM900v1_genomic.fna                          | 151           |
| GCF_900004855.1_BB155_genomic.fna                             | 152           |
| GCF_000027045.1_ASM2704v1_genomic.fna                         | 239           |
| GCF_000705395.1_Staphylococcus_aureus_SA3-LAU_WGS_genomic.fna | 291           |
| GCF_000239655.1_ASM23965v2_genomic.fna                        | 395           |
| GCF_000009585.1_ASM958v1_genomic.fna                          | 398           |
| GCF_001049595.1_Staphylococcus_aureus_AH2_genomic.fna         | 630           |
| GCF_000189435.1_ASM18943v2_genomic.fna                        | 700           |
| GCF_000153665.1_ASM15366v1_genomic.fna                        | 1159          |
| GCF_000363005.1_Stap_aure_M1216_V1_genomic.fna                | 1613          |
| GCF_000174575.1_ASM17457v1_genomic.fna                        | 1892          |
| GCF_000751255.1_Sa_FSA084_genomic.fna                         | 2022          |
| GCF_000606605.1_Stap_aure_M21126_V1_genomic.fna               | 2250          |

#### Determination of MLST
A in house script was written to BLAST each genome against each of the MLST loci. Only the top BLAST hit was retained. MLST 
was only determined if the top hits were identical hits across the full length of the gene.

##### Example
```
mkdir -p data/mlst
cd data/mlst

# Use 'getmlst.py' from SRST2, to download current MLST info
getmlst.py --species "Staphylococcus aureus"

# Make BALST databases for each loci
mkdir blastdb
ls *.tfa | sed 's/.tfa//' | xargs -I {} makeblastdb -in {}.tfa -out blastdb/{} -dbtype nucl -title "{} database"

# Run the in house script 'determine-mlst.py'
cd ../../
mkdir -p data/blast-results/
python ../scripts/determine-mlst.py data/completed-genomes/ data/blast-results/ \
                                    ~/bin/blastn data/mlst/blastdb data/mlst/saureus.txt > data/genome-mlst.txt

```

#### Simulation of FASTQ files
We used ART to simulate varying coverages of HiSeq 2000 reads for each of the 25 genomes.

##### Example
```
mkdir -p data/simulated-reads

# Simulate using ART
scripts/simulate-reads.sh data/completed-genomes/ data/simulated-reads/ scripts/art/art_illumina
```
