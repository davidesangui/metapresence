# metapresence
Calculation of metrics for evaluating the distribution of aligned reads on a genome.

Starting from an indexed sorted bam file (with the index in the same folder) and a list of genomes, metapresence evaluates the evenness of the distribution of the reads by calculating various metrics:
- coverage per genome. The coverage is calculated as the mean of the depth at each position. When calculating the per-base depth, cigar operations and quality scores are ignored.
- breadth per genome and BER metric.  Both the actual breadth and the ratio between breadth and expected breadth (BER) are returned. The expected breadth is calculated using the formula reported in https://instrain.readthedocs.io/en/latest/important_concepts.html ,section 6: detecting organisms in metagenomic data.
- FUG metric. The FUG metric is calculated as the as normalized value of the cumulative distribution function of the frequencies of distances between consecutive reads for values greater than the expected distance. The output reports FUG metric for both group of mates when the reads are paired-end. 

## Output

The output consists in a tsv file containing the metrics at the genome level, and in a png file containing the scatterplot between BER and average FUG per each genome.

## Usage
**Help page:**
  `python3 metapresence.py -h`
  
**Command line**
   `python3 metapresence.py [options] contig_genomes indexed_sorted.bam all_sequences.fasta`

### Mandatory:
- all_sequences.fasta: a fasta file containing all the genomes that are to be analyzed. 
- indexed_sorted_bam: indexed sorted bam file of the alignment on the fasta file. The sorted bam file must be indexed and the index (.bai) must be in the same
  folder (the index is needed by pysam). 
- contig_genomes: Either a file with each line listing a scaffold and a bin name, tab-seperated, or the path to a folder containing all (and only) the fasta files of the genomes to analyze. 

### Options:
- -p [int], number of processes to be executed in parallel. Default: 1 .
- -o the prefix of the coverm output and of the output files. Default: metrics
- --unpaired set this flag if the aligned reads are not paired. Default: FALSE
## Dependencies
Metapresence requires only the Python modules Pysam and Numpy. It was tested using the following versions of these packages, but it should work fine with any more recent version.
- numpy 1.23.3 (https://www.numpy.org)
- pysam 0.19.1 (https://github.com/pysam-developers/pysam)



