# metapresence
Calculation of metrics for evaluating the distribution of aligned reads on a genome.

Starting from an indexed sorted bam file (with the index in the same folder) and a list of genomes, metapresence evaluates the randomness of the distribution of the reads by calculating different metrics, some of which are returned relative to their expected value:
- coverage per genome and per scaffold. The coverage is calculated as the mean of the depth at each position. When calculating the per-base depth, the only cigar operation that is potentially considered is the clipping, the others are ignored.
- breadth per genome and per scaffold.  Both the actual breadth and the ratio between breadth and expected breadth are returned. The expected breadth is calculated using the formula reported in https://instrain.readthedocs.io/en/latest/important_concepts.html ,section 6: detecting organisms in metagenomic data.
- average distance between the mapping position of two consecutive reads among all the possible pairs (in case of paired-end reads, only the first mate encountered in the sorted bam file is considered). This metric is returned as a ratio real/expected, where the expected value is given by the length of the scaffold divided by the number of reads considered.
- average distance between the mapping positions of two reads of randomly sampled pairs, divided by the length of the scaffold (in case of paired-end reads, only the first mate encountered in the sorted bam file is considered). This metric is returned as a ratio real/expected, where the expected value is 1/3.

The two latter metrics are calculated at the level of a single scaffold. The whole-genome values are given by an average of the values of each scaffold weighted for the ratio between the length of the scaffold and the length of the genome.

## Output

The output consists in two tsv files, one containing the metrics at the level of single scaffolds, the other at the genome level, marked by the suffix "scaffolds" or "genomes" respectively.

## Usage
**Help page:**
  `python3 metapresence.py -h`
  
**Command line**
   `python3 metapresence.py [options] scaffold_to_bin.txt indexed_sorted.bam all_sequences.fasta`

### Mandatory:
- all_sequences.fasta: a fasta file containing all the genomes that are to be analyzed.
- indexed_sorted_bam: indexed sorted bam file of the alignment on the fasta file. The sorted bam file must be indexed and the index (.bai) must be in the same
  folder (the index is needed by pysam). 
- scaffold_to_bin.txt: A file with each line listing a scaffold and a bin name, tab-seperated. 

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
- numpy 1.23.3 (https://www.numpy.org)
- pysam 0.19.1 (https://github.com/pysam-developers/pysam)

