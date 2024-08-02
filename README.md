# Metapresence
Calculation of metrics for evaluating the distribution of aligned reads onto nucleotidic sequences.  
While Metapresence was developed for the purpose of identifying present species in shotgun metagenomic sequencing data, it can potentially be used to assess the evenness of read distribution onto any kind of DNA/RNA sequences starting from the result of a sequence alignment, provided that the sequences are in FASTA format and the alignment is reported as a sorted and indexed BAM file.  

The associated article can be found at: https://doi.org/10.1128/msystems.00213-24   
If you use Metapresence, please cite via: 
*Sanguineti D, Zampieri G, Treu L, Campanaro S.0.Metapresence: a tool for accurate species detection in metagenomics based on the genome-wide distribution of mapping reads. mSystems0:e00213-24.https://doi.org/10.1128/msystems.00213-24*

## Installation
The Python packages numpy and pysam are needed to use Metapresence. You can either create a conda environment and install them in the following way:
```
# create conda environment
conda create metapresence_env

# activate conda environment
conda activate metapresence_env

# install packages
conda install anaconda::numpy
conda install bioconda::pysam
```
Or you can install them with pip:
```
pip install numpy
pip install pysam
```
Once you have installed the two dependencies, simply download the file < metapresence.py > from this repository, and run it with Python:
```
python3 metapresence.py -h
```

## Usage
The toy_example will be used to illustrate the basic usage of Metapresence. In this example dataset, a reference database comprising seven genomes is present. 250,000 sequencing reads were generated with CAMISIM (https://github.com/CAMI-challenge/CAMISIM/) from three of these genomes (those without the suffix "_added").  
We start with a set of reference genomes and shotgun metagenomic sequencing reads. The first thing to do is to align the sequencing reads on the reference genomes to generate a sorted bam file, which can be then indexed. Bowtie2 (https://github.com/BenLangmead/bowtie2) and Samtools (https://github.com/samtools/samtools) will be used:
```
# generate a multifasta from reference genomes
cat ref_genomes/* > allgenomes.fasta

# generate bowtie2 index
bowtie2-build --threads 4 allgenomes.fasta allgenomes_index

# align reads with bowtie2, convert to bam and sort with samtools
bowtie2 -1 toy_R1.fastq.gz -2 toy_R2.fastq.gz -x allgenomes_index -p 4 | samtools view -b -@ 4 | samtools sort -@ 4 > alignment.sorted.bam

# generate aligment index (alignment.sorted.bam.bai) with samtools
samtools index -@ 4 alignment.sorted.bam 
```
Metapresence can now be launched giving it as inputs the folder containing the reference genome files in fasta format, and the sorted bam file:
```
python3 metapresence.py ref_genomes alignment.sorted.bam -p 4 -o example
```
< -p > option defines the number of processes to be used for parallelization. < -o > is the prefix of output files.  
Two output files are generated:  
 - **example_metrics.tsv**: this file stores, for each genome, coverage, breadth, breadth-expected breadth ratio (BER), read-distance metric (FUG) for both group of mates, and the number of mapping reads. Only the reference genomes with at least one mapped read will be in this file.
```
genome  length  coverage        breadth BER     FUG1    FUG2    read_count
Abiotrophia_2_added.fasta       2015532 0.01265174653639833     0.00963666168535156     0.8674391645602579      0.3533318744284423      0.3531612227095277 170
Abiotrophia_defectiva.fasta     2046826 0.9947596913465043      0.5889850920400659      1.0076016570018371      0.6319770867879021      0.6315087011082681 13574
Acidithiobacillus_173_added.fasta       2086976 0.00014374865834585544  0.00012122803520500475  0.9551380029939979      0.15377361609639206     0.15384764661996755 2
Acidithiobacillus_45_added.fasta        3675789 0.016078180766088586    0.009586241212430855    0.6800329266380903      0.12016570884526044     0.12016434170080421 394
Acidithiobacillus_ferridurans.fasta     2921399 2.647156379529123       0.9024282544082476      0.9988963071103649      0.6278898204513724      0.6266337387192086  51556
Xanthobacter_10_added.fasta     5216997 0.026567007801614608    0.011320880575549497    0.4882709529452287      0.0254046698676058      0.025422533966759663924 924
Xanthobacter_dioxanivorans.fasta        6650762 4.135420873578095       0.9731994619563894      0.9991269898195253      0.6356478161084415      0.6355469590475508  183358
```
 - **example_abundances.tsv**: this file stores the relative abundance (%) of each of the genomes identified as present using BER and FUG values. The relative abundance is defined as the coverage of a given genome divided by the sum of coverages of all present genomes.
```
genome  relative_abundance_%
Abiotrophia_defectiva.fasta     12.79049240699158
Acidithiobacillus_ferridurans.fasta     34.03679689378635
Xanthobacter_dioxanivorans.fasta        53.17271069922207
```
## Advanced usage
By default, Metapresence uses certain values of coverage, BER, FUG and number of mapped reads to define a sequence as "present" and add it to the abundances file. To change the default values, just change the corresponding flag. In particular, the default criteria to define a sequence as present are the following:  
- If less than 80 reads map on a given sequence, it is considered as absent. Corresponding flag: ```-min_reads```
- If a sequence has BER value lower than 0.8, it is considered as absent. Corresponding flag: ```-ber_threshold```
- If a sequence has a coverage lower than 0.1 and FUG values for both the group of mates lower than 0.5, it is considered as absent. To change the fug threshold, use the flag ```-fug_threshold``` . To change the maximum coverage above which fug metric is no longer used (the "0.1" just mentioned), use the flag ```-max_for_fug```. To change whether both FUG values (both group of mates), only one FUG value or the mean FUG value must be above the threshold, use the flag ```-fug_criterion```.
For example:
```
python3 metapresence.py ref_genomes alignment.sorted.bam -min_reads 30 -ber_threshold 0.7 -fug_threshold 0.4 -max_for_fug 0.2 -fug_criterion mean
```
If the reads are not paired, the flag ```--unpaired``` should be set. In this case, the flag ```-fug_criterion``` is irrelevant.  

The input sequences can be either a folder of fasta files, or a fasta file with one or more contigs. If the case is the latter, a text file should be provided with the flag ```-input_bins```. This file is necessary to cluster different contigs as belonging to the same sequence (for instance, a genome). In this text file, each line should list the names of a contig and of the sequence to which it belongs, tab-separated.  
If each contig should be treated as an independent sequence, then each line of this text file should list the name of a contig and, tab-separated, again the name of the contig. More simply, the same result can be obtained by setting the flag ```--all_contigs```. 

### Parsing metric output
The time consuming step of metapresence.py is that of parsing the bam file to calculate the metric values. On the contrary, inferring relative abundances from the metric values is straightforward and fast.  
It may be useful to test different metrics thresholds depending on the specific instance under investigation. In this case, it is useless and time-consuming to re-calculate the metric values by running metapresence.py. In this repository, the script parse_metrics.py is available in the utils subdirectory. This script takes as input the metric output of metapresence.py (*_metrics.tsv) and outputs an abundance file based on the metric threshold parameters given by the user.
```
python3 parse_metrics.py -metrics example_metrics.tsv -output_abundances new_abundances.tsv -min_reads 30 -ber_threshold 0.7 -fug_threshold 0.4 -max_for_fug 0.2 -fug_criterion mean
```
### Merging abundance outputs
In the utils subdirectory of this repository, the script merge_abundances.py is available. This script takes as input a folder containing (only) abundance outputs of metapresence.py, and it outputs a tab-separated table where the different samples (abundance outputs) are merged together.
```
python3 merge_abundances.py -abundance_folder samples/ -output_table merged_table.tsv
```
