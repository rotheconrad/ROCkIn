# ROCkIn

The preparation pipeline for building ROCker models with [ROCkOut](https://github.com/KGerhardt/ROCkOut)

This pipeline walks the researcher through the process of collecting the necessary sequence information needed to build and refine ROCker models for any functional gene group of interest. The steps involve a combination of Python, Bash, and bioinformatics tools. Examples are provided for PBS and Sbatch job schedulers for users with access to a compute cluster.

To use the provided PBS or Sbatch scripts replace "PATH/to/GitHub/repo" with the path to your local copy of this Github repo, and replace all the "YOUR_PROMPTs" with the relevant information.

The steps are left separately so the user can more easily follow the workflow, and so individual steps can be more efficiently parallelized depending on the users system.

## Dependencies

- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Clustal Omega](http://www.clustal.org/omega/)
- [TrimAl](http://trimal.cgenomics.org/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
- [Python](https://www.python.org/)

 #### References

 1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC bioinformatics. 2009 Dec;10(1):1-9.
 1. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology. 2017 Nov;35(11):1026-8.
 1. Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75
 1. Salvador Capella-Gutiérrez, José M. Silla-Martínez, Toni Gabaldón, trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses, Bioinformatics, Volume 25, Issue 15, 1 August 2009, Pages 1972–1973
 1. Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490.
 1. A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014, open access.
 1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

#### Python Packages

- [xlmtramp2](https://github.com/tBaxter/xmltramp2)
- [pandas](https://pandas.pydata.org/) 
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [scipy](https://scipy.org/)
- [sklearn](https://scikit-learn.org/stable/)
- [Bio](https://biopython.org/)
- [hdbscan](https://hdbscan.readthedocs.io/)
- [pycircos](https://github.com/ponnhide/pyCircos)

#### References

1. Aaron Swartz, Kristian Glass, T. Carter Baxter. https://github.com/tBaxter/xmltramp2, 2002
1. McKinney W, others. Data structures for statistical computing in python. In: Proceedings of the 9th Python in Science Conference. 2010. p. 51–6.
1. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature. 2020;585:357–62.
1. Hunter JD. Matplotlib: A 2D graphics environment. Computing in science & engineering. 2007;9(3):90–5.
1. Waskom ML. Seaborn: statistical data visualization. Journal of Open Source Software. 2021 Apr 6;6(60):3021.
1. Pauli Virtanen, et. al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
1. Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.
1. Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
1. L. McInnes, J. Healy, S. Astels, hdbscan: Hierarchical density based clustering In: Journal of Open Source Software, The Open Journal, volume 2, number 11. 2017

 *Python, it's packages, and all program above can be installed with [Conda](https://docs.conda.io/en/latest/miniconda.html) except for pycircos needs to use "pip install python-circos".*


# PART 00: Curate reference sequences

ROCker model building starts with a set of curated reference sequences. Curation of these sequences is up to the researcher building the model and they should be either experimentally verified or highly specific. **You must start with a reliable set of gene sequences for gene function you wish to study.**

Specialized databases are a good starting resources such as the [NCBI ref gene database](https://www.ncbi.nlm.nih.gov/pathogens/refgene) for antibiotic resistance genes.

When curating sequences, it is insightful to look at a multiple sequence alignment, a phylogenetic tree such as a quick neighbor joining tree, and/or a clustered heatmap of sequence similarity. There are many approaches to this. We will outline a quick and easy one here utilizing EBI's website and we provide a Python script to build a sequence similarity heatmap from results.

 1. [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Select Pearson/FASTA as the output format in step 2.
 2. Download the alignment file and view it with your favorite multiple sequence alignment tool such as [AliView](https://ormbunkar.se/aliview/).
 3. Under the "Results Viewers" tab, select "Send to Simple Phylogeny" at the bottom of the options. In STEP 2 of Simple Phylogeny, turn on the distance correction, exclude gaps, and P.I.M. (percent identity matrix) options.
 4. At this point you should see a phylogram of the results. You can scroll down and "View Phylogenetic Tree File" which is in newick format (.nwk). You can save this file and view/edit it with tools such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTol](https://itol.embl.de/).
 5. For the sequence similarity heatmap, look under the "Results Summary" tab at the top of the Simple Phylogeny results page. Download the "Percent Identity Matrix" file (.pim) and use the 00a_PIM_clustered_heatmap.py Python script included in the 02_Python directory of this GitHub repo to create a heatmap figure.

 ```bash
 python /Path/to/GitHub/repo/02_Python/00a_PIM_clustered_heatmap.py -i your_files_name.pim -o name_your_output_file.pdf

 # for script options and info
 python /Path/to/GitHub/repo/02_Python/00a_PIM_clustered_heatmap.py -h
 ```

 ![Example Figure of sequence similarity heatmap](https://github.com/rotheconrad/ROCkIn/blob/main/05_Example_Figs/00_Example-A.png)

Once you have made your selections, create a fasta formatted file with the amino acid sequences and short meaningful defline names. We will refer this file of curated sequences as RefSeqs.faa. Example file can be found in the 06_Example_Files directory of this repo.

# PART 01: UniProt sequence search

Since ROCker models are used with metagenomics data, we want to account for broad sequence diversity around the reference sequences. To do this, we will perform a BLAST search of UniProt's SwissProt and TrEMBL databases.

If there are multiple reference sequences with ≥90% or ≥95% sequence similarity it is suggested to select one representative sequence to use in the BLAST search as it is unlikely they will yield different search results. This will save computational time and make the figures easier to interpret.

EBI webservices provides code to access and Blast search the UniProt database programmatically. Similar to searching through the UniProt website, this will query the database remotely and use EBI resources to perform the Blast search. The results are returned as a separate text files contianing one UniProt ID per match for each sequence in the input fasta file (RefSeqs.faa). For more information see [https://www.ebi.ac.uk/Tools/webservices/](https://www.ebi.ac.uk/Tools/webservices/).

```bash
# path to ebi blast script
ebi_blast='Path/to/GitHub/repo/02_Python/01a_ebi_ncbiblast.py'

# executes the command Run 30 sequences at a time from 1 fasta file returning 1000 matches (alignments) per query sequence.
python ${ebi_blast} --email YOUR_EMAIL --program blastp --stype protein \
--sequence RefSeqs.faa --database uniprotkb --multifasta --useSeqId --maxJobs 30 --pollFreq 60 \
--outformat ids --exp 10 --alignments 1000
```

PBS example:
```bash
qsub -v fasta=RefSeqs.faa /Path/to/GitHub/repo/01a_PBS/01a_ebi_blast.pbs
```

Sbatch example:
```bash
sbatch --export fasta=RefSeqs.faa /Path/to/GitHub/repo/01b_Sbatch/01a_ebi_blast.sbatch
```

Once we have the Blast results, there are a few additional steps to get the corresponding fasta sequences that we need. First, we will use dbfetch to retrieve the .dat file for each UniProt ID returned from the blast search.

```bash
# file house keeping
mkdir 01a_ebi_blast_ids 01b_ebi_dat 01c_ebi_fasta
mv *ids.txt 01a_ebi_blast_ids

# path to dbfetch script
dbfetch='Path/to/GitHub/repo/02_Python/01b_ebi_dbfetch.py'
# loop over ids.txt files and run dbetch fetchData on each UniProt ID.
for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; while read p; do name=`echo $p | cut -d: -f2`; python ${dbfetch} fetchData $p > ${odir}/${gene}/${name}.dat; echo $name; done < $f; done
```

PBS example:
```bash
# Download the *.dat files from EBI with dbfetch
for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; qsub -v input=${f},odir=${odir},gene=${gene} /Path/to/GitHub/repo/01a_PBS/01b_ebi_dbfetch.pbs; done
```

Sbatch example:
```bash
# Download the *.dat files from EBI with dbfetch
for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; sbatch --export input=${f},odir=${odir},gene=${gene}  /Path/to/GitHub/repo/01b_Sbatch/01b_ebi_dbfetch.sbatch; done
```

Now we will parse the .dat file to retrieve the corresponding fasta sequences, and we will use the fasta deflines to store additional relevant information.

```bash
for d in 01b_ebi_dat/*; do n=`basename $d`; echo $n; python /Path/to/GitHub/repo/02_Python/01d_parse_dat_file.py -i $d -o 01c_ebi_fasta/${n}.fasta; done

# concatenate to single fasta file
cat 01c_ebi_fasta/MCR-* >> 01c_ebi_fasta/ALL_EBI_BLAST_MATCHES.faa
```

At this point we have a single file (ALL_EBI_BLAST_MATCHES.faa) containing all the fasta sequences from the UniProt database that returned a match to our curated RefSeqs.faa.

# PART 02: Select sequences for model training

Part 2 is divided into the following series of steps:

1. Remove duplicate sequences from ALL_EBI_BLAST_MATCHES.faa.
1. Filter the sequences to remove poor matches.
1. Cluster the sequences and select a representative sequence for each cluster.
1. Generate a multiple sequence alignment (MSA).
1. Trim/clean the multiple sequence alignment.
1. Build phylogenetic tree.
1. Compute branch distance and cluster sequences into clades.
1. Generate annotation file and tree figure labeled by cluster.

Our goal is to build a phylogenetic tree with our curated sequences (RefSeqs.faa) and known surrounding sequence diversity (ALL_EBI_BLAST_MATCHES.faa) that we will use to make decisions about positive and negative sequence sets for training the ROCker model.

#### Step 1: Remove duplicates

Since we have multiple verified sequences that we searched to UniProt, we likely have found overlapping search results. I wrote a Python script to deduplicate the concatenated fasta.

```bash
# setup directories
mkdir 02a_tree_prep
# Remove duplicates in fasta file
python /Path/to/GitHub/repo/02_Python/02a_Remove_Duplicate_Fasta_Entries.py -f 01c_ebi_fasta/ALL_EBI_BLAST_MATCHES.faa -o 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa
```

#### Step 2: Filter

The initial sequence search does not have great options to filter the matches that are returned. As such, it usually returns many spurious matches, and many nearly identical matches to our curated sequences. We don't need either of these types of sequences to build our ROCker models.

To filter our search results, I run Blastp locally using the RefSeqs.faa as the reference database and the Blast search results as the query. Then I use that to filter which sequences to keep.

I wrote a Python script that applies 2 filters:
	1 Removes matches <= 30% identity to verified sequences
	1 Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

This script outputs 100% sequence matches separately as well in tabular blast, fasta, and in a text list of Reference sequences with their corresponding list of uniprot IDs with 100% sequence identity. Ideally, we find at least 1 uniprot ID for each reference sequence. If not, we can look through the filtered blast output and find the closest matches to our reference sequences available in the uniprot database. This important because we'll want to include 1 exact match (or closest match) uniprot ID in our final input list to ROCkOut. During the similar sequence clustering these ID's can sometimes be dereplicated as they aren't guaranteed selection as a cluster representative.

The script also plots histograms of each parameter post filtering and can be rerun with different filtering options if neccessary.

```bash
# make blast database from RefSeqs.faa
makeblastdb -dbtype prot -in RefSeqs.faa
# Run blast search
blastp -num_threads 2 -max_target_seqs 10 \
-db RefSeqs.faa -query 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa -out 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.blast  -subject_besthit -outfmt \
'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
# Set python script directory
scripts=/Path/to/GitHub/repo/02_Python
# filter blast results
python ${scripts}/02b_Blastp_filter_hist.py -i 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.blast -o 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.fltrd.blast
# retrieve fasta sequences for the filtered blast results
python ${scripts}/02c_Get_Fasta_from_Filtered_Blast.py -b 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.fltrd.blast -q 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa -o 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa
```

![Example histogram figures for Blast sequence alignments.](https://github.com/rotheconrad/ROCkIn/blob/main/05_Example_Figs/02_Example-B.png)

PBS example:
```bash
qsub -v ref=RefSeqs.faa,qry=02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa,out=02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa /Path/to/GitHub/repo/01a_PBS/02b_Blastp.pbs
```

Sbatch Example:
```bash
sbatch --export ref=RefSeqs.faa,qry=02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa,out=02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa /Path/to/GitHub/repo/01b_Sbatch/00b_BlastP.sbatch
```

#### Step 3: Cluster and select representative sequences

If you have fewer than hundreds of searched sequences (FILTER_EBI_BLAST_MATCHES.faa) at this point you can skip this step. This step reduces the number of sequences by clustering them at 90% amino acid sequence similarity and choosing one reference sequence for each cluster.

Cluster with MMSeqs2
```bash
# setup directories
mkdir 02b_mmseqs
mkdir 02b_mmseqs/temp
# create mmseqs database
mmseqs createdb 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa 02b_mmseqs/MyDB01 --dbtype 1
# Cluster input sequences at 90% amino acid sequence identity
mmseqs cluster 02b_mmseqs/MyDB01 02b_mmseqs/MyDB02 02b_mmseqs/temp --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign
# convert to tsv format
mmseqs createtsv 02b_mmseqs/MyDB01 02b_mmseqs/MyDB01 02b_mmseqs/MyDB02 02a_tree_prep//mmseqs_90.tsv
# Get the representative sequence for each cluster in fasta format
mmseqs createsubdb 02b_mmseqs/MyDB02 02b_mmseqs/MyDB01 02b_mmseqs/MyDB03
mmseqs convert2fasta 02b_mmseqs/MyDB03 02a_tree_prep/mmseqs_reps.fasta
# Count the number of representative sequences
grep -c '>' 02a_tree_prep/mmseqs_reps.fasta
```

PBS Example:
```bash
qsub -v infile=02a_tree_prep/01_MCR_fltr_ebi_blast_fltrd.fasta,odir=02a_tree_prep ../00b_PBS/02c_mmseqs.pbs
``` 

Sbatch Example:
```bash
sbatch --export infile=02a_tree_prep/01_MCR_fltr_ebi_blast_fltrd.fasta,odir=02a_tree_prep /Path/to/GitHub/repo/01b_Sbatch/27c_MMSeqs2_Cluster.sbatch
```

# PART 03: Create and trim multiple align. Build Phylogenetic tree and create clades.

Part 3 is divided into the following series of steps:

1. Generate a multiple sequence alignment (MSA).
1. Trim/clean the multiple sequence alignment.
1. Build phylogenetic tree.
1. Compute branch distance and cluster sequences into clades.
1. Generate annotation file and tree figure labeled by cluster.

Our goal is to build a phylogenetic tree with our curated sequences (RefSeqs.faa) and known surrounding sequence diversity (ALL_EBI_BLAST_MATCHES.faa) that we will use to make decisions about positive and negative sequence sets for training the ROCker model.


#### Step 1: Multiple Sequence Alignment

We need a good multiple sequence alignment to build the phylogenetic tree. We'll start with an alignment using only the curated sequences (RefSeqs.faa). We can check this alignment with something like [AliView](https://ormbunkar.se/aliview/) and make manual adjustments if necessary, then we'll align the searched sequences (mmseqs_reps.fasta) to this alignment.

```bash
# First align the curated sequences
clustalo -i RefSeqs.faa -o 02a_tree_prep/RefSeqs.aln
# then align the searched sequences to that alignment
clustalo -i 02a_tree_prep/mmseqs_reps.fasta -o 02a_tree_prep/my_MSA.aln --profile1 02a_tree_prep/RefSeqs.aln
```

PBS Example:
```bash
qsub -v verified=RefSeqs.faa,newseqs=02a_tree_prep/mmseqs_reps.fasta /Path/to/GitHub/repo/01a_PBS/02d_seq_alignment.pbs
``` 

Sbatch Example:
```bash
sbatch --export verified=RefSeqs.faa,newseqs=02a_tree_prep/mmseqs_reps.fasta /Path/to/GitHub/repo/01b_Sbatch/02d_seq_alignment.sbatch
```

#### Step 2: Trim/Clean Multiple Sequence Alignment

Before building a phylogenetic tree, MSAs should be cleaned to removed spurious sequences and columns consisting mainly of gaps. TrimAl's algorithm does a pretty good and quick job, but it is still recommended to inspect the alignment manually as well with something like [AliView](https://ormbunkar.se/aliview/).

TrimAl sometimes removes reference sequences or does other strange things. Please check the mmseqs_reps.fasta.aln file before TrimAl and the trimmed TrimAl output before proceeding. It is quite common for step 4 and 5 to require manual work to achieve a good alignment file. A good alignment file is necessary to get a reliable tree for step 6, 7, & 8.

Sometimes you may want to consider going back to step 2 and increasing the filtering parameters or to step 3 and cluster at 85 or 95 instead of 90 with mmseqs.

We're using the following two parameters for TrimAl ([User Guide](http://trimal.cgenomics.org/)):
    - resoverlap: Minimum overlap of a positions with other positions in the column to be considered a "good position". Range: [0 - 1].
    - seqoverlap: Minimum percentage of "good positions" that a sequence must have in order to be conserved. Range: [0 - 100].

```bash
trimal -in 02a_tree_prep/my_MSA.aln -out 02a_tree_prep/my_MSA_trimmed.aln -resoverlap 0.75 -seqoverlap 80 -automated1

# count sequences after trimming
grep -c '>' 02a_tree_prep/my_MSA_trimmed.aln
```

At his point the sequence names look something like this:
>A0A8E0FMZ4_ECOLX//Arylsulfatase//Escherichia

Which is really hard to look at on the tips of the tree. This script splits at the '//' leaving only the uniprot id to be display at the tree leaves.
>A0A8E0FMZ4_ECOLX

```bash
# clean up sequence names for the tree.
# this script alters the file in place. no new file created.
python /Path/to/GitHub/repo/02_Python/02e_clean_seq_names.py -i 02a_tree_prep/my_MSA_trimmed.aln
```

PBS Example:
```bash
qsub -v input=02a_tree_prep/my_MSA.aln,output=02a_tree_prep/my_MSA_trimmed.aln /Path/to/GitHub/repo/01a_PBS/02e_seq_trim.pbs
``` 

Sbatch Example:
```bash
sbatch --export input=02a_tree_prep/my_MSA.aln,output=02a_tree_prep/my_MSA_trimmed.aln /Path/to/GitHub/repo/01b_Sbatch/02e_seq_trim.sbatch
```

#### Step 3: Build Phylogenetic Tree

Now that we have a good multiple sequence alignment, we are ready to build the tree. Here is an example with FastTree and RAxML. FastTree is recomended unless you are using a cluster. RAxML has very long run times. This step will provide you with Newick formated files than can be viewed/edited with tools such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTol](https://itol.embl.de/).


```bash
FastTree 02a_tree_prep/my_MSA_trimmed.aln > 02a_tree_prep/Tree_FastTree.nwk
``` 

PBS Example:
```bash
# bootstrapped ML tree w/ RAxML
qsub -v input=02a_tree_prep/my_MSA_trimmed.aln,output=02a_tree_prep/Tree_Bootstrap /Path/to/GitHub/repo/01a_PBS/02f_RAxML_AminoAcid-Bootstrap.pbs

# Fasttree - approximate maximum likelihood tree
qsub -v input=02a_tree_prep/my_MSA_trimmed.aln,output=02a_tree_prep/Tree_FastTree.nwk /Path/to/GitHub/repo/01a_PBS/02f_FastTree.pbs 
```

Sbatch Example:
```bash
sbatch --export input=02a_tree_prep/my_MSA_trimmed.aln,output=02a_tree_prep/Tree_Bootstrap /Path/to/GitHub/repo/01b_Sbatch/002f_RAxML_AminoAcid-Bootstrap.sbatch

sbatch --export input=02a_tree_prep/my_MSA_trimmed.aln,output=02a_tree_prep/Tree_FastTree.nwk /Path/to/GitHub/repo/01b_Sbatch/02f_FastTree.sbatch
```

#### Step 4: Compute branch distance and cluster sequences into clades.
#### Step 5: Generate annotation file and tree figure labeled by cluster.

*Step 4 & 5 are combined.*

Create an organized and annotated tsv file to explore clustered clades and make decisions for ROCkOut inputs. The annotation information is pulled from the fasta deflines (where we stored the information from the .dat files back at the begining) and organized into columns. The branch distances are calculated between all sequences and the HDBSCAN algorithm is used to cluster the sequences into initial clades.

This step creates a .tsv (Easily opened in Excel) that contains the UniProt IDs, assigned cluster, functional annotation, and taxonomic classification for each sequence. It aslo plots a tree figure with the leaf nodes color labelled by assigned cluster. The .tsv file follows the order of the tree figure. The sequences in the .tsv file can be reodered/reassigned based on visual inspection of the tree figure, and the researcher can decide clades for positive and/or negative reference sequences based on these results. Once the researcher has made their selections, simply copy and paste the UniProt IDs into a positive.txt and a negative.txt file (one ID per line). Then insert the text files into the [ROCkOut](https://github.com/KGerhardt/ROCkOut) pipeline.

Essentially, it is up to the researcher to make good choices based on the results. The HDBSCAN clustering algorithm is just to provide initial suggestions for plausible clades/clusters. The researcher should double check these results and make changes/decisions accordingly.

*The HDBSCAN algorithm assigns many -1's to sequences that fall outside what the algorithm has determined to be the cluster boundaries. The algorithm isn't perfect at building the clades, it is just a quick way to get started. The researcher can quickly reassign -1's into the appropriate clades while making their selections*

*DBSCAN is a clustering algorithm that works by finding data points that are all within a "step" of at least one other point and putting them into a cluster. A cluster is therfoew made up of points that can all be walked to by taking steps on points within the cluster. When no new points can be reached by taking a single step, a new cluster is started.*

```bash
# concatenate the curated sequences and the dereplicated and filtered searched sequences
cat RefSeqs.faa 02a_tree_prep/mmseqs_reps.fasta >> 02a_tree_prep/final_sequence_set.fasta

# setup directory
mkdir 02c_Annoted_Tree

# Set python script directory
scripts=/Path/to/GitHub/repo/02_Python

# convert nwk to distance matrix, cluster the matrix and add annotations

### 02g_Tree_Distance_Cluster.py can take a while (several hours) - especially for larger trees.

# For RAxML Bootstrap
python ${scripts}/02g_Tree_Distance_Cluster.py -i 02a_tree_prep/Tree_Bootstrap_bestTree.MCL_bootstrap.RAxML -f 02a_tree_prep/final_sequence_set.fasta -o 02c_Annoted_Tree/Cluster_annotations

python ${scripts}/02i_Plot_Annotated_Tree_v2.py -n 02a_tree_prep/Tree_Bootstrap_bestTree.MCL_bootstrap.RAxML -a 02c_Annoted_Tree/Cluster_annotations_annotated.tsv -o ROCker_prep_tree.pdf

# or for Fasttree
python ${scripts}/02g_Tree_Distance_Cluster.py -i 02a_tree_prep/Tree_FastTree.nwk -f 02a_tree_prep/final_sequence_set.fasta -o 02c_Annoted_Tree/Cluster_annotations

python ${scripts}/02i_Plot_Annotated_Tree_v2.py -n 02a_tree_prep/02a_tree_prep/Tree_FastTree.nwk -a 02c_Annoted_Tree/Cluster_annotations_annotated.tsv -o ROCker_prep_tree.pdf
```

![Example phylogenetic tree labelled by assigned cluster/clade.](https://github.com/rotheconrad/ROCkIn/blob/main/05_Example_Figs/07_Example-C.png)

# PART 04: Make UniProt ID selections.

ROCkOut takes as input a list of UniProt IDs. pos.txt and option neg.txt depending on your reference sequence, functional target, and phylogenetic tree. Some gene targets may only need positive target sequences (pos.txt) while other genes clearly will have a closely related but off target clade appear in the phylogenetic and sequence similarity analyses above. In which case they can be used as neagetive targets (neg.txt). pos.txt and optional neg.txt are simply single column text files with one uniprot entry per line.

Once you've made your selection, proceed to ROCkOut!

