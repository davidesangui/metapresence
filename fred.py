import sys

if len(sys.argv)==1:
    print('syntax: [options] annotations bam fasta binning.txt') # input.txt is the coverm-like input for contig binning
    print('mandatory: (in this order)')
    print('\tannotations: emapper annotation file\n\tbam: sorted and indexed bam file with the index in the same folder\n\tfasta: fasta of the alignment in the bam file\n\tbinning.txt: a text file with each line listing a scaffold and a bin name, tab-seperated')
    print('options:\n\t-f [ko,module] the functions to analyze, comma-separated (no space). default: all (ko,module)')
    print('\t-o output prefix, default: fred')
    exit()


import pysam
import math
import statistics

# options setting
default={'ko':True,'module':True,'o':'fred'}
options=' '.join(sys.argv[1:-4]).split('-')
solved=False
while not solved:
    found=False
    for i in range(len(options)):
        if len(options[i].split())==1:
            found=True
            options[i-1]=options[i-1]+'-'+options[i]
            options.pop(i)
            break
    if not found:
        solved=True
for opt in options:
    a=opt.split()
    if len(a)==0:
        continue
    if opt[0]=='f':
        default['ko']=False
        default['module']=False
        functions=a[1].split(',')
        for funct in functions:
            default[funct]=True
        continue
    try:
        default[a[0]]
        default[a[0]]=a[1]
    except KeyError:
        print('-'+a[0], 'is not an available option')
        print('help:',sys.argv[0])
        exit()
#########################################################################################

def rel_ab(bam,fasta,genome_contigs):
    print('Calculating relative abundances')
    print('\treading fasta file...') # calculating contig lengths from fasta file and saving them into a dictionary
    f=open(fasta)
    contigs_lengths_d={}
    contig_name=''
    for line in f:
        if line[0]=='>':
            header=line[1:].strip().split()
            contig_name=header[0]
            contigs_lengths_d[contig_name]=0
        else:
            contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()
    print('\tdone')

    genomes=[]                
    total_length_allgenomes=0 # calculation of each genome length for normalization purposes
    genome_length={} 
    for genome in genome_contigs:
        genomes.append(genome)
        genome_length[genome]=0
        for contig in genome_contigs[genome]:
            genome_length[genome]+=contigs_lengths_d[contig]
        total_length_allgenomes+=genome_length[genome]

    print('\treading bam file...') # reading bam file (with pysam) for read numbers on each genome 
    sam=pysam.AlignmentFile(bam,"rb") 
    genome_readnumbers={} # read numbers are to be saved into a dictionary
    total_reads=0
    for genome in genomes:
        reads_number=0
        for contig in genome_contigs[genome]:
            for line in sam.fetch(contig=contig): #iterating over pysam.fetch to count the reads mapping on each contig of the genome
                reads_number+=1
        genome_readnumbers[genome]=reads_number
        total_reads+=reads_number
    print('\tdone')

    norm_reads_total=0       # relative abundance calculation (RA = genome_normalized_reads / total_normalized_reads)
    for genome in genomes:
        norm_reads=genome_readnumbers[genome]*(total_length_allgenomes/genome_length[genome]) 
        genome_readnumbers[genome]=norm_reads
        norm_reads_total+=norm_reads
    genome_ra={}

    for genome in genome_readnumbers:
        genome_ra[genome]=genome_readnumbers[genome]/norm_reads_total
    print('   Relative abundances calculated\n')
    return(genome_ra)

#############################################################################

# innanzitutto assegno a ciascuna ko un index in un dictionary, o ciascun module un tuple in un dictionary: (indice, lista di ko che lo formano)
def index_function(annot_file,which): # which sarebbe quale funzione considerare (per ora ko o moduli). 
    function_index={}
    if which=='ko' or which=='KO':    
        f=open(annot_file)
        index=0
        for line in f:
            if line[0]!='#':
                    fields=line.strip().split()
                    for field in fields:
                        if field[:3]=='ko:':
                            for ko in field.split(','): 
                                try:                                     
                                    function_index[ko]
                                except KeyError:
                                    function_index[ko]=index
                                    index+=1
        f.close()
    elif which=='module' or which=='Module':
        f=open(annot_file)
        index=0
        for line in f:
            if line[0]!='#':
                    fields=line.strip().split()
                    for field in fields:
                        if field[:2]=='M0': # cambiare modo per riconoscere il field dei modules??
                            for module in field.split(','): 
                                try:                                     
                                    function_index[module]
                                except KeyError:
                                    # qui inserire la cosa per verificare quali ko formano un certo modulo
                                    function_index[module]=index # e aggiungere il tuple 
                                    index+=1

    return(function_index)

# quindi un indice a tutti i genomi sempre in un dictionary (servirà per la jaccard distance)
def index_genome(genome_contig):
    genome_index={}
    i=0
    for genome in genome_contig:
        genome_index[genome]=i
        i+=1
    return(genome_index)

# quindi faccio un dictionary dove per ogni genoma (key) metto una lista con 1 se ha la funzione a quell'index e zero se non ce l'ha
def gen_listof_function(function_index,annot_file,which,gene_ko_1or0): # l'ultimo parametro serve per verificare la completezza del modulo, mettere None se which è ko
    gen_function_1or0={}
    number_of_functions=len(function_index)
    if which=='ko' or which=='KO':
        f=open(annot_file)
        for line in f:
            if line[0]!='#':
                    fields=line.strip().split()
                    genome=fields[0].split('.')[0] # capire come definire generalizzando il genoma 
                    for field in fields:
                        if field[:3]=='ko:':
                            for ko in field.split(','): 
                                try:                                     
                                    gen_function_1or0[genome][function_index[ko]]=1
                                except KeyError:
                                    gen_function_1or0[genome]=[0]*number_of_functions
                                    gen_function_1or0[genome][function_index[ko]]=1
        f.close()
    elif which=='module' or which=='Module':
        f=open(annot_file)
        for line in f:
            if line[0]!='#':
                    fields=line.strip().split()
                    genome=fields[0].split('.')[0] # capire come definire generalizzando il genoma
                    for field in fields:
                        if field[:2]=='M0': # cambiare modo per riconoscere il field dei modules??
                            for module in field.split(','):
                                # qui inserire la cosa per verificare la completezza del modulo nel genoma
                                try:                                     
                                    gen_function_1or0[genome][function_index[module]]=1
                                except KeyError:
                                    gen_function_1or0[genome]=[0]*number_of_functions
                                    gen_function_1or0[genome][function_index[module]]=1
    return(gen_function_1or0)

# la funzione per calcolare la jaccard distance fra tutti i genomi e tenerla in un dictionary (ogni genoma key ha una lista di tutti i genomi con jd nell'index dato da genome_index)
def jaccard_distance(gen_function_1or0,genome_index):
    genome_jds={}
    for gen1 in gen_function_1or0:
        genome_jds[gen1]=[0]*len(genome_index)
        function_gen1=gen_function_1or0[gen1]
        for gen2 in gen_function_1or0:
            if gen1==gen2:
                genome_jds[gen1][genome_index[gen2]]=0
                continue
            function_gen2=gen_function_1or0[gen2]
            intersection=0
            not_intersection=0                     # la not_intersection è tutto ciò che esiste e non è intersection (ko presenti in una OR l'altra); quindi: union = intersection + not_intersection
            for i in range(len(function_gen1)):
                if function_gen1[i]==1 and function_gen2[i]==1:
                    intersection+=1
                elif function_gen1[i]!=function_gen2[i]:
                    not_intersection+=1
            genome_jds[gen1][genome_index[gen2]]=not_intersection/(not_intersection+intersection) # union - intersection = not_intersection + intersection - intersection = not_intersection
    return(genome_jds)

# la funzione per calcolare tutte le metriche associate a una singola funzione
def single_functions(function_index,genome_index,gen_function_1or0,gen_ra,genome_jd,total_jd):
    topr='' # ritorna una stringa con tutte le metriche. ordine: ko, numero, proporzione, jd media di chi ha la funzia, proporzione corretta per la jd, proporzione corretta ma con le rel ab, shannon index, simpson index, prodotto abbondanze, stdev abbondanze(per controllo)

    # conta e proporzione
    for function in function_index:
        topr+=function
        gen_with_function=[] # ricordo quali genomi hanno la funzione 
        abundances=[]  # e mi segno parallelamente l'abundance (per la trait contribution)
        index=function_index[function]
        for gen in gen_function_1or0:
            gen_funct_profile=gen_function_1or0[gen]
            if gen_funct_profile[index]==1:
                gen_with_function.append(gen)
                abundances.append(gen_ra[gen])
        topr+='\t'+str(len(gen_with_function))+'\t'+str(len(gen_with_function)/len(gen_function_1or0)) # la conta e la proporzione sono segnate
        
    # segue il calcolo della proporzione corretta per la jd
        i=0 # l'indice per sapere dove sono arrivato nella lista di genomi col ko
        intra_jd=0
        while i<len(gen_with_function)-1:
            gen1=gen_with_function[i]
            for j in range(i+1,len(gen_with_function)):
                gen2=gen_with_function[j]
                intra_jd+=genome_jd[gen1][genome_index[gen2]]
            i+=1
        if len(gen_with_function)!=1:
            intra_jd=intra_jd/((len(gen_with_function)*(len(gen_with_function)-1))/2)
        topr+='\t'+str(intra_jd)+'\t'+str((len(gen_with_function)/len(gen_function_1or0))*(intra_jd/total_jd)) # proporzione corretta segnata

    # segue il calcolo della proporzione corretta ma usando le abbondanze
        i=0 # l'indice per sapere dove sono arrivato nella lista di genomi col ko
        intra_jd=0
        total_ra=0
        while i<len(gen_with_function)-1:
            gen1=gen_with_function[i]
            total_ra+=gen_ra[gen1]
            for j in range(i+1,len(gen_with_function)):
                gen2=gen_with_function[j]
                intra_jd+=genome_jd[gen1][genome_index[gen2]]
            i+=1
        if len(gen_with_function)!=1:
            intra_jd=intra_jd/((len(gen_with_function)*(len(gen_with_function)-1))/2)
        total_ra+=gen_ra[gen_with_function[-1]]
        topr+='\t'+str(total_ra*(intra_jd/total_jd)) # proporzione corretta e usando le rel ab segnata
        
    
    # segue il calcolo degli indici di TCE
        total_abundance=sum(abundances)
        num_shan=0
        num_simp=0
        for ab in abundances: # shannon index
            num_shan+= -1*ab/total_abundance*math.log(ab/total_abundance)
            num_simp+=(ab/total_abundance)**2
        if len(abundances)==1:
            topr+='\t'+'0'+'\t'+'0'
        else:
            topr+='\t'+str(num_shan/math.log(len(gen_function_1or0)))+'\t'+str(1-num_simp) # shannon index e simpson index segnati

    # segue il prodotto delle abbondanze pairwise (in ambi i sensi)
        metric=0
        for i in range(len(gen_with_function)):
            for j in range(len(gen_with_function)):
                if i!=j:
                    metric+=abundances[i]*abundances[j]
        topr+='\t'+str(metric)
        
    # segue la std delle abbondanze per controllo
        if len(abundances)==1:
            topr+='\t'+'0'+'\n'
        else:
            topr+='\t'+str(statistics.stdev(abundances))+'\n'
    return(topr) # magari cambiare il return

# la funzione per calcolare la functional redundancy a livello di comunità (tutte le metriche)
def all_functions(gen_function_1or0,genome_jds,gen_ra,genome_index,function_index):
    topr='#' # per ora, ritorna una stringa con tutte le community-level metriche. ordine: Tian\tRicotta\tEsponente\tNumeromedio 

    # tian e ricotta
    value_tian=0
    num_ricotta=0
    den_ricotta=0  
    for gen_i in genome_jds:
        gen_i_jds=genome_jds[gen_i]
        gen_i_ra=gen_ra[gen_i]
        den_ricotta+=gen_i_ra*(1-gen_i_ra)
        for gen_j in genome_jds:
            num_ricotta+=gen_i_ra*gen_ra[gen_j]*gen_i_jds[genome_index[gen_j]]
            if gen_i==gen_j:
                continue
            value_tian+=(1-gen_i_jds[genome_index[gen_j]])*gen_i_ra*gen_ra[gen_j]
    topr+=str(value_tian)+'\t'+str(1-(num_ricotta/den_ricotta))

    # esponente
    average_noffunctions=0
    for gen in gen_function_1or0:
        average_noffunctions+=sum(gen_function_1or0[gen])
    average_noffunctions=average_noffunctions/len(gen_function_1or0)
    fr=math.log(len(function_index),average_noffunctions*len(gen_function_1or0))
    topr+='\t'+str(fr)

    # average number of functions per genome
    topr+='\t'+str(len(gen_function_1or0)/len(function_index))

    return(topr) # magari cambiare il return

###############################################################################################

f=open(sys.argv[-1]) # generation of genome_listofcontig dictionary from input.txt
genome_contigs={}
for line in f:
    a=line.strip().split('\t')
    try:
        genome_contigs[a[1]].append(a[0])
    except KeyError:
        genome_contigs[a[1]]=[a[0]]
f.close()

gen_ra=rel_ab(sys.argv[-3],sys.argv[-2],genome_contigs) # relative abundances

genome_index=index_genome(genome_contigs) # indexing genomes

if default['ko']:
    print('Calculating functional redundancy for ko')

    ko_index=index_function(sys.argv[-4],'ko') # indexing ko
 
    gen_ko_1or0=gen_listof_function(ko_index,sys.argv[-4],'ko',None) # data structure for all genomes and ko

    genome_jds=jaccard_distance(gen_ko_1or0,genome_index) # pairwise jaccard distance for each pair of genome

    total_jd=0 # la media della jd su tutti i genomi è la somma di tutte le liste del dictionary diviso n*(n-1) perchè le conto tutte due volte (e quando è uguale non si conta)
    for gen in genome_jds:
        total_jd+=sum(genome_jds[gen])
    total_jd=total_jd/(len(genome_contigs)*(len(genome_contigs)-1))

    new=open(default['o']+'_ko.tsv','w')

    print('#Tian\tRicotta\tEsponente\tNumeromedio',file=new)
    print(all_functions(gen_ko_1or0,genome_jds,gen_ra,genome_index,ko_index),file=new)

    print('ko\tnumero\tproporzione\tjd_media\tproporzione _corretta_per_la_jd\tprop_corretta_e_rel_ab\tshannon _index\tsimpson_index\tprodotto_abbondanze\tstdev_abbondanze',file=new)
    print(single_functions(ko_index,genome_index,gen_ko_1or0,gen_ra,genome_jds,total_jd),file=new)

    print('\tdone')
    new.close()

if default['module']:
    print('Calculating functional redundancy for module')

    if not default['ko']: # se non è già stato fatto, faccio quello che serve ai moduli a livello di ko (cioè per vedere la module completeness, ancora da fare) 
        ko_index=index_function(sys.argv[-4],'ko') 
        gen_ko_1or0=gen_listof_function(ko_index,sys.argv[-4],'ko',None)

    module_index=index_function(sys.argv[-4],'module')

    gen_module_1or0=gen_listof_function(module_index,sys.argv[-4],'module',None)

    genome_jds=jaccard_distance(gen_module_1or0,genome_index)

    total_jd=0 
    for gen in genome_jds:
        total_jd+=sum(genome_jds[gen])
    total_jd=total_jd/(len(genome_contigs)*(len(genome_contigs)-1))

    new=open(default['o']+'_module.tsv','w')

    print('#Tian\tRicotta\tEsponente\tNumeromedio',file=new)
    print(all_functions(gen_module_1or0,genome_jds,gen_ra,genome_index,module_index),file=new)

    print('ko\tnumero\tproporzione\tjd_media\tproporzione _corretta_per_la_jd\tprop_corretta_e_rel_ab\tshannon _index\tsimpson_index\tprodotto_abbondanze\tstdev_abbondanze',file=new)
    print(single_functions(module_index,genome_index,gen_module_1or0,gen_ra,genome_jds,total_jd),file=new)

    print('\tdone')
    new.close()












