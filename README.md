# Metapresence
Calculation of metrics for evaluating the distribution of aligned reads onto nucleotidic sequences.  
While Metapresence was developed for the purpose of identifying present species in shotgun metagenomic sequencing data, it can potentially be used to assess the evenness of read distribution onto any kind of DNA/RNA sequences starting from the result of a sequence alignment, provided that the sequences are in FASTA format and the alignment is reported as a sorted and indexed BAM file.  

The associated article can be found at: ...   
If you use Metapresence, please cite via: ...

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
We start with a set of reference genomes and shotgun metagenomic sequencing reads. With any sequence aligner, the sequencing reads can be aligned on the reference genomes. The resulting SAM file has to be converted into a sorted and indexed BAM file, which can be done with samtools (https://github.com/samtools/samtools):
```
samtools view alignment.sam -b -o alignment.bam
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
```
Metapresence can now be launched giving it as inputs the folder containing the reference genome files in fasta format, and the sorted bam file:
```
python3 metapresence.py refgen_folder alignment_sorted.bam -p 4 -o example
```
< -p > option defines the number of processes to be used for parallelization. < -o > is the prefix of output files.  
Two output files are generated:  
 - example_metrics.tsv: this file stores, for each genome, coverage, breadth, breadth-expected breadth ratio (BER), read-distance metric (FUG) for both group of mates, and the number of mapping reads. Only the reference genomes with at least one mapped read will be in this file.
 - example_abundances.tsv: this file stores the relative abundance (%) of each of the genomes identified as present using BER and FUG values. The relative abundance is defined as the coverage of a given genome divided by the sum of coverage of all present genomes.
## Advanced usage
By default, Metapresence uses certain values of coverage, BER, FUG and number of mapped reads to define a sequence as "present" and add it to the abundances file. To change the default values, just change the corresponding flag. In particular, the default criteria to define a sequence as present are the following:  
- If less than 80 reads map on a given sequence, it is considered as absent. Corresponding flag: ```-min_reads```
- If a sequence has BER value lower than 0.8, it is considered as absent. Corresponding flag: ```-ber_threshold```
- If a sequence has a coverage lower than 0.1 and FUG values for both the group of mates lower than 0.5, it is considered as absent. To change the fug threshold, use the flag ```-fug_threshold``` . To change the maximum coverage above which fug metric is no longer used (the "0.1" just mentioned), use the flag ```-max_for_fug```. To change whether both FUG values (both group of mates), only one FUG value or the mean FUG value must be above the threshold, use the flag ```-fug_criterion```.
For example:
```
python3 metapresence.py refgen_folder alignment_sorted.bam -min_reads 30 -ber_threshold 0.7 -fug_threshold 0.4 -max_for_fug 0.2 -fug_criterion mean
```
If the reads are not paired, the flag ```--unpaired``` should be set. In this case, the flag ```-fug_criterion``` is irrilevant.  

The input sequences can be either a folder of fasta files, or a fasta file with one or more contigs. If the case is the latter, a text file should be provided with the flag ```-input_bins```. This file is necessary to cluster different contigs as belonging to the same sequence (for instance, a genome). In this text file, each line should list the names of a contig and of the sequence to which it belongs, tab-separated.  
If each contig should be treated as an independent sequence, then each line of this text file should list the name of a contig and, tab-separated, again the name of the contig. More simply, the same result can be obtained by setting the flag ```--all_contigs```. 




