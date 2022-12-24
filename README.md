# metapresence
Calculation of metrics for evaluating the distribution of aligned reads on a genome.

Starting from an indexed sorted bam file (with the index in the same folder) and a list of genomes, metapresence evaluates the randomness of the distribution of the reads by calculating different metrics:
- coverage and breadth per genome using coverm (version 0.6.1 or above) launched via inStrain quick_profile (version 1.5.7). Both the actual breadth and the ratio between breadth and expected breadth are returned. The expected breadth is calculated using the formula reported in https://instrain.readthedocs.io/en/latest/important_concepts.html ,section 6: detecting organisms in metagenomic data.
- average distance between the mapping position of two consecutive reads among all the possible pairs (in case of paired-end reads, only the first mate encountered in the sorted bam file is considered). This metric is returned as a ratio real/expected, where the expected value is given by the length of the scaffold divided by the number of reads considered.
- average distance between the mapping positions of two reads of randomly sampled pairs, divided by the length of the scaffold (in case of paired-end reads, only the first mate encountered in the sorted bam file is considered). This metric is returned as a ratio real/expected, where the expected value is 1/3.

The two latter metrics are calculated at the level of a single scaffold. The whole-genome values are given by an average of the values of each scaffold weighted for the ratio between the length of the scaffold and the length of the genome.

If needed, instead of giving as input a sorted bam file, it is possible to launch bowtie2 (version 2.3.5.1) via metapresence.py. The output of bowtie2 will be piped into samtools (version 1.9 or above) to be converted to bam and sorted. This sorted bam file will be used as input of metapresence.py .
## Output

The output consists both in the output of inStrain quick_profile, that is a folder marked by the suffix "coverm", and in two tsv files, one containing the metrics at the level of single scaffolds, the other at the genome level, marked by the suffix "scaffolds" or "genomes" respectively.

## Usage
**Help page:**
  `python3 metapresence.py -h`

**With bowtie2:**
  `python3 metapresence.py bowtie2 alignment launch ! [options] -s coverm_input.txt all_sequences.fasta`
  The exclamation mark separates the command for bowtie2 and the one for metapresence.py.
  
**Without bowtie2:**
  `python3 metapresence.py (indexed)sorted.bam ! [options] -s coverm_input.txt all_sequences.fasta`
  The exclamation mark separates the input sorted bam file and the command for metapresence.py.

### Mandatory:
- all_sequences.fasta: a fasta file containing all the genomes that has to be analyzed.
- sorted_bam/bowtie2: either an indexed sorted bam file of the alignment on the fasta file, or the command line of bowtie2 to do it. In case, the sorted bam file must
  be indexed and the index (.bai) must be in the same folder (the index is needed by pysam). If bowtie2 is launched, the index will be automatically created.
- -s coverm_input.txt: input file for coverm. A file with each line listing a scaffold and a bin name, tab-seperated. 

### Options:
- -p [int], number of processes to be executed in parallel. Default: 1 .
- -l [int] minimum length of a contig to be included in the distance-metrics calculation. Default: 0 .
  Scaffolds that are shorter than the assigned value will not be considered. They will not weight in the genomic average and the genome length will be recalculated
  subtracting theirs.
- -v [float] value of the distance metrics to be given to a scaffold with less than three reads mapping on it. Default: 0.5 .
  For scaffolds with less than three unique mates mapping on it the distance metrics will not be calculated, but they will weight in the genomic average with the
  assigned value
- -o the prefix of the coverm output and of the .tsv output files. Default: metrics

## Dependencies
- inStrain 1.5.7
- coverM 0.6.1
- samtools 1.14
- bowtie 2 2.3.5.1

