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

# metapresence modificato
- La cosa diversa è che ho cambiato le metriche di distanza. Ho lasciato solo la distanza tra reads successive, ma fatta in maniera tale che distanze inattese pesino proporzionalmente alla loro dimensione.
- La metrica è calcolata per ciascuna mate. La distanza attesa è la lunghezza della sequenza diviso il numero di reads (una o l'altra mate) che vi mappano.
- Per prima cosa unisco tutte le contig in un'unica sequenza continua, salvo tutte le posizioni in cui le reads mappano su questa sequenza in una lista, e quindi metto una read inventata all'inizio e alla fine di questa sequenza.
- Quindi definisco una finestra ampia tanto quanto è la distanza attesa fra due reads più uno. mi muovo con questa finestra una base alla volta e per ogni base vedo se c'è o non c'è una read dentro alla finestra. Se c'è non faccio niente, se non c'è tolgo 1 ad una variabile che parte da un valore che è uguale alla lunghezza della sequenza. Questo per ogni base fino alla fine della sequenza continua (meno la dimensione della finestra).
- Lo score finale è la variabile a cui ho tolto i vari 1 diviso la lunghezza della sequenza. Se tutte le read mappano entro la distanza attesa, lo score finale sarà uno.
- Maggiore è il gap tra la distanza attesa e la distanza effettiva, più si abbasserà lo score. Così, mentre per genomi presenti lo score sta intorno a 0.65, per genomi assenti tende a essere più basso perchè, anche se gruppi di reads possono mappare uniformemente in regioni conservate, tenderanno a mappare a macchie e quindi i gap inaspettatamente grandi saranno visti e pesati proporzionalmente.
- Se le reads mappassero uniformemente in un'unica regione, la metrica avrebbe un valore corretto. Per questo metto una read all'inizio e alla fine della sequenza.
- In realtà non faccio veramente questa cosa base per base, ma per ogni posizione registrata tolgo alla variabile: distanza con successiva - dimensione finestra. è la stessa cosa ma più veloce.



# functional redundancy (fred.py)
Per adesso calcolo tutte le metriche che ho raccolto e le mando a un file tsv, sia quelle "community-level" sia quelle per la singola funzione. le singole funzioni per ora sono i ko e i moduli (anche se i moduli per ora li considero sempre presenti, poi bisognerà definire se ci sono o no in base alla completezza)
## Script
* Nella prima parte dello script calcolo le relative abundances a partire dal bam, dal fasta e dal file di input che serve a binnare in genomi le varie contig presenti nel fasta. la relative abundance la calcolo normalizzando le reads che mappano su un genoma cosi: Rnorm = R * (sum_all_genome_length / genome_length ). La relative abundance è quindi (Rnorm / sum_all_Rnorm). Alla fine le abbondanze sono ritornate con un dizionario (genoma : abbondanza).
  * La velocità della prima parte dipende dal numero di reads mappate nel bam (a meno che non sia un fasta file molto grande che diventa significativo nel tempo della prima parte). con le 100 milioni di reads di cami_medium ci mette un paio di minuti (su ion)
- Nella seconda parte invece analizzo le funzioni dal file emapper.annotations.
  - considerando ad esempio i ko, innanzitutto assegno a ciascun ko che incontro nel file un indice in un dizionario (ko : indice) che sarà la sua posizione in una lista che corrisponde al profilo funzionale del genoma
  - quindi assegno a ciascun genoma il profilo funzionale in un dizionario (genoma : profilo). il profilo funzionale è una lista di n posizioni (n = numero di ko nel file), dove ogni posizione corrisponde a un ko definito dall'indice assegnato prima. Ogni posizione della lista è 1 se il genoma ha il ko, 0 se non ce l'ha. Per assegnare ciascun ko al genoma a cui appartiene bisogna per forza che il nome di ciascun gene nel file di annotazioni contenga il nome del genoma da cui proviene (o al limite della contig, il genoma della quale può essere ripescato dal file di input contig-genoma).
  - quindi, a partire da questo dizionario di profili funzionali, calcolo l'overlap funzionale fra tutte le coppie di genomi e lo salvo in un altro dizionario (genoma : lista di overlap); nella lista di overlap ogni posizione corrisponde a un genoma, avendo preventivamente indexato i genomi analogamente a quanto fatto coi ko.
  - l'overlap funzionale lo calcolo come jaccard distance (1 - intersezione/unione, la jaccard distance perchè lo fa Tian per lo stesso scopo https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7719190/ ) iterando parallelamente sopra i profili funzionali della coppia di genomi.
  - con le relative abundances, i profili funzionali e la jaccard distance di ogni coppia di genomi posso calcolare tutte le metriche, sia per la singola funzione che community-level.
- La velocità della seconda parte dipende dal numero di genomi e dal numero di funzioni considerate. coi 132 genomi di cami_medium e coi soli ko ci mette meno di un minuto (su ion).
- lo stesso procedimento può essere ripetuto per tutte le funzioni presenti. Per ora l'ho fatto solo per i moduli, ma per i moduli bisogna appunto aggiustare in base alla loro completezza, cosa che ancora non ho messo.
- l'output per ora è un file tsv (uno per funzione) con le prime due linee che sono i valori community-level, e poi i valori di single function per ciascuna funzione.
## Dependencies
- pysam 0.19.1

