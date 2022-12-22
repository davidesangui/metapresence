import sys

if len(sys.argv)==1:
    print('usage: bowtie2 alignment launch ! [options] -s coverm_input.txt all_sequences.fasta')
    print('use the exclamation mark to separate the bowtie2 command line and the other part')
    print('if bowtie2 is not needed, just put the indexed sorted bam file before the exclamation mark')
    print('for help:',sys.argv[0],'-h')
    exit()
if len(sys.argv)==2 and sys.argv[1]=='-h':
    print('usage: bowtie2 alignment launch ! [options] -s coverm_input.txt all_sequences.fasta')
    print('use the exclamation mark to separate the bowtie2 command line and the other part')
    print('if bowtie2 is not needed, just put the indexed sorted bam file before the exclamation mark')
    print('')
    print('mandatory:')
    print('')
    print(' before the exclamation mark:')
    print('     either the bowtie2 command line or an indexed sorted bam file, with bam and index (.bai) in the same directory')
    print('')
    print(' after the exclamation mark:')
    print('     -s input file for coverm. A file with each line listing a scaffold and a bin name, tab-seperated')
    print('     all_sequences.fasta: the fasta file on which the reads are aligned')        
    print('')    
    print('options:')
    print(' -p [int]','number of processes to be executed in parallel. Default: 1')
    print(' -l [int] minimum length of a contig to be included in the distance-metrics calculation. Default: 0')
    print(' -v [float] value of the distance metrics to be given to a scaffold with less than three reads mapping on it. Default: 0.5')
    print(' -o the prefix of the coverm output and of the .tsv output files. Default: metrics')
    print(' -h','help')
    exit()

import pandas as pd
import os
import random
import numpy as np                              
import pysam
import multiprocessing    

# Options and parameters setting
default={'p':1,'o':'metrics','l':0,'v':0.5,'s':''}
command=' '.join(sys.argv[1:])
bowtie2_launch=command.split('!')[0]
other_launch=command.split('!')[1]
options=other_launch.split('-')
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
    try:
        default[a[0]]
        default[a[0]]=a[1]
    except KeyError:
        print('-'+a[0], 'is not an available option')
        print('help:',sys.argv[0],'-h')
processes=int(default['p'])
output=default['o']
min_scaffold_length=int(default['l'])
penalty=float(default['v'])
coverm_input=default['s']

# The function to calculate the distances between all the possible paris of reads.
# It is used when few reads map on a certain contig, so that it is less convenient do a random samplings multiple times.
def all_distances(array): 
    which=-1
    sum_infunct=0    
    for i in range(len(array)):
        which+=1
        for j in range(which+1,len(array)):            
            sum_infunct+=abs(array[i]-array[j])
    return(sum_infunct/((len(array)*(len(array)-1))/2))


# The function to divide the lists of genomes in n sublists, with n as the number of processes. 
def chunks(l, processi):
    lists=[]
    sizes=len(l)//processi+1
    for i in range(0,len(l),sizes):
        lists.append(l[i:i+sizes])
    return(lists)


# The function that calculates all the metrics on a sublist of genomes starting from the coverm output and the bam file.
# Each process handles one sublist of genome and creates a tmp file. All the tmp files will be orderly joined together (thanks to the second parameter 
# of the function) by the last process to end.
def main(list_of_genomes,order,contig_dict,bam):
    if len(bam.split())>1:
        sam=pysam.AlignmentFile(output+"_sorted.bam","rb")
    else:
        sam=pysam.AlignmentFile(bam.split()[0],"rb")    
    z=open(str(order)+'tmp_multiall_metrics_genomes.tsv','w')
    zz=open(str(order)+'tmp_multiall_metrics_scaffolds.tsv','w')
    if order==0: # the process handling the first sublist prints the header in the final output
        print('genome\tlength\tcoverage\tbreadth\tbr_realexp_ratio\tnext_realexp_ratio\trandom_realexp_ratio',file=z)
        print('genome\tlength\tcoverage\tbreadth\tnext_realexp_ratio\trandom_realexp_ratio',file=zz)
    for seq in list_of_genomes:            
        to_print=seq       
        f=open(output+'_coverm/coverm_raw.tsv') #some values of coverm are taken
        contig_length={}
        contig_cov={}
        contig_breadth={}
        contig_lengths=[]       
        for line in f:
            a=line.split('\t')
            if a[0]=='Sample':
                continue
            try:
                if contig_dict[a[1]]==seq:                       
                    contig_length[a[1]]=int(a[4])
                    contig_breadth[a[1]]=str(int(a[3])/int(a[4]))
                    contig_cov[a[1]]=a[2]
            except KeyError:
                continue            
        seq_length=df.query("genome==@seq")['length'].item() # from coverm the genome length is taken               
        to_print+='\t'+str(seq_length)+'\t'+str(df.query("genome==@seq")['coverage'].item())+'\t'+str(df.query("genome==@seq")['breadth'].item())

        # creation of the list with all the map positions of the first mate of each fragment.   
        next_read_dist=[]
        ratio_50000=[]    
        for contig in contig_length:
            if contig_length[contig]<min_scaffold_length: # Exclusion of too short contig (how much depends on the option)
                seq_length-=contig_length[contig]
                continue    
            pair_taken={}
            positions=[]                        
            for line in sam.fetch(contig=contig):    
                a=str(line).split('\t')            
                try:
                    pair_taken[a[0]]                    
                    continue
                except KeyError:                    
                    pair_taken[a[0]]=1                                               
                positions.append(int(a[3]))
            contig_lengths.append(contig_length[contig])            
            if len(positions)<3: # if less than three fragments map on the contig, it will count in the wheighted average according to the penalty (option)
                next_read_dist.append(penalty)
                ratio_50000.append(penalty)
                continue
                 
            
            # Average distance between consecutive fragments
            sum_1=0
            for i in range(len(positions)-1):
                sum_1+=positions[i+1]-positions[i]
            next_read_ratio=(sum_1/(len(positions)-1))/(contig_length[contig]/len(positions))        
            next_read_dist.append(next_read_ratio) # each contig has its real/expected ratio

            # average distance between reads of randomly sampled pairs (fifty-thousands pairs)
            if len(positions)>316: # if it is convenience (316*315/2 is around fifty-thousand), all the distances are calculated.
                sum=0
                for i in range(50000):
                    a=random.choice(positions)
                    while 1: # to avoid the sampling of the same read
                        b=random.choice(positions)
                        if b!=a: 
                            break
                    sum+=abs(a-b)            
                random_ratio=(sum/50000)/(contig_length[contig]*(1/3)) 
            else:
                tmp=all_distances(positions)
                random_ratio=tmp/(contig_length[contig]*(1/3))

            beb=str(float(contig_breadth[contig])/(1-np.exp(-0.883*float(contig_cov[contig])))) 
            # the distance metrics and those calculated by coverm are written in the contig output, and the distace metrics are added to the list for the weighted average
            print(contig+'\t'+str(contig_length[contig])+'\t'+contig_cov[contig]+'\t'+contig_breadth[contig]+'\t'+str(next_read_ratio)+'\t'+str(random_ratio),file=zz)        
            ratio_50000.append(random_ratio)

           
        for i in range(len(contig_lengths)):  # Each contig value is weighted for the ratio between the contig length and the total one
            next_read_dist[i]=next_read_dist[i]*(contig_lengths[i]/seq_length)
            ratio_50000[i]=ratio_50000[i]*(contig_lengths[i]/seq_length)
        next_read_dist.append(penalty*((seq_length-np.sum(contig_lengths))/seq_length)) # contig without reads will count in the weighted average, the penalty is the value
        ratio_50000.append(penalty*((seq_length-np.sum(contig_lengths))/seq_length)) 

        beb=str(df.query("genome==@seq")['breadth'].item()/(1-np.exp(-0.883*df.query("genome==@seq")['coverage'].item()))) # breadth and expected breadth ratio calculation
        nextratio=str(np.sum(next_read_dist)) #final calculations of the distance metrics
        randomratio=str(np.sum(ratio_50000))
        to_print+='\t'+beb+'\t'+nextratio+'\t'+randomratio
        print(to_print,file=z)    
    queue.put(0,block=False) # when a process ends, it puts something in a queue which has max size equal to the number of processes
    if queue.full(): # the process that finds a full queue is the last to end, it deletes join tmp files into the final output and deletes them.  
        z.close()
        zz.close()        
        tmp_scaf=''
        tmp_gen=''
        for i in range(processes):
            tmp_scaf+=str(i)+'tmp_multiall_metrics_scaffolds.tsv '
            tmp_gen+=str(i)+'tmp_multiall_metrics_genomes.tsv '        
        os.system("cat "+tmp_scaf+">"+output+"_scaffolds.tsv")
        os.system("cat "+tmp_gen+">"+output+"_genomes.tsv")
        os.system("rm -f "+tmp_scaf)
        os.system("rm -f "+tmp_gen)
        os.system("rm -f tmp_qp.txt")
    

        

# The function to give to each process its sublist 
def multi_main(array,processori):
    blocchi=chunks(array,processori) 
    jobs=[]
    for i in range(len(blocchi)):        
        j=multiprocessing.Process(target=main,args=(blocchi[i],i,all_contigs_and_bin,bowtie2_launch))
        jobs.append(j)
    for j in jobs:        
        j.start()   


# generation of a all-contigs dictionary from the coverm input
f=open(coverm_input)
all_contigs_and_bin={}
for line in f:
    a=line.strip().split('\t')
    all_contigs_and_bin[a[0]]=a[1]

# bowtie2 and coverM via inStrain
if len(bowtie2_launch.split())>1: # if it is requested, launch also bowtie2, otherwise just launch coverm with the given bam file 
    os.system(bowtie2_launch+" | samtools sort -O bam > "+output+"_sorted.bam")
    os.system("samtools index -b "+output+"_sorted.bam")
    os.system("inStrain quick_profile -p "+str(default['p'])+" -s"+coverm_input+" -o "+output+"_coverm "+output+"_sorted.bam"+" "+sys.argv[len(sys.argv)-1])
else:
    os.system("inStrain quick_profile -p "+str(default['p'])+" -s"+coverm_input+" -o "+output+"_coverm "+bowtie2_launch.split()[0]+" "+sys.argv[len(sys.argv)-1])





# Main list of genomes that will be spit in sublists for each process
df=pd.read_csv(output+'_coverm/genomeCoverage.csv')
genomes=[] 
for i in df['genome']:
    genomes.append(i)


if __name__ == "__main__": 
    queue = multiprocessing.Queue(maxsize=processes) # the queue is needed for understanding when an ending process is the last to end              
    multi_main(genomes,processes) 
    
    
















        



