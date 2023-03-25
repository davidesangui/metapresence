import sys

if len(sys.argv)==1:
    print('usage: [options] scaffolds_to_bin.txt indexed_sorted.bam all_sequences.fasta')
    print('for help:',sys.argv[0],'-h')
    exit()
if len(sys.argv)==2 and sys.argv[1]=='-h':
    print('usage: [options] scaffolds_to_bin.txt indexed_sorted.bam all_sequences.fasta')
    print('')
    print('mandatory:')
    print('  in this order:')
    print('  scaffold_to_bin.txt: file with each line listing a scaffold and a bin name, tab-seperated')
    print('  indexed_sorted.bam: indexed sorted bam file, with bam and index (.bai) in the same directory')
    print('  all_sequences.fasta: the fasta file on which the reads are aligned')        
    print('')    
    print('options:')
    print(' -p [int]','number of processes to be executed in parallel. Default: 1')
    print(' -l [int] minimum length of a contig to be included in the distance-metrics calculation. Default: 0')
    print(' -v [float] value of the distance metrics to be given to a scaffold with less than three reads mapping on it. Default: 0.5')
    print(' -o the prefix of the .tsv output files. Default: metrics')
    print(' -h','help')
    exit()


import os
import random
import numpy as np                              
import pysam
import multiprocessing    

# Options and parameters setting
default={'p':1,'o':'metrics','l':0,'v':0.5}
command=' '.join(sys.argv[1:])
options=' '.join(sys.argv[1:-3]).split('-')
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
scaffold_input=sys.argv[-3]
sorted_bam=sys.argv[-2]

# The function to calculate the distances between all the possible paris of reads.
# It is used when few reads map on a certain contig, so that it is less convenient doing a random samplings multiple times.
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

# The function to give to each process its sublist 
def multi_main(array,processori):
    blocchi=chunks(array,processori) 
    jobs=[]
    for i in range(len(blocchi)):        
        j=multiprocessing.Process(target=main,args=(blocchi[i],i,sorted_bam,contigs_lengths_d,genome_contigs))
        jobs.append(j)
    for j in jobs:        
        j.start()


#contig lengths calculation into a dictionary
print('reading fasta file')
f=open(sys.argv[-1])
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
print('finished')
# generation of genome_listofcontig dictionary from input.txt
f=open(scaffold_input)
genome_contigs={}
for line in f:
    a=line.strip().split('\t')
    try:
        genome_contigs[a[1]].append(a[0])
    except KeyError:
        genome_contigs[a[1]]=[a[0]]
f.close()

# Main list of genomes that will be spit in sublists for each process
genomes=[] 
for i in genome_contigs:
    genomes.append(i)

def main(list_of_genomes,order,bam,contig_lengths_d,genome_contigs):
    if len(bam.split())>1:
        sam=pysam.AlignmentFile(output+"_sorted.bam","rb")
    else:
        sam=pysam.AlignmentFile(bam.split()[0],"rb")    
    z=open(str(order)+'tmp_multiall_metrics_genomes.tsv','w')
    zz=open(str(order)+'tmp_multiall_metrics_scaffolds.tsv','w')
    if order==0: # the process handling the first sublist prints the header in the final output
        print('genome\tlength\tcoverage\tbreadth\tbr_realexp_ratio\tnext_realexp_ratio\trandom_realexp_ratio\tread_count',file=z)
        print('scaffold\tlength\tcoverage\tbreadth\tbr_realexp_ratio\tnext_realexp_ratio\trandom_realexp_ratio\tread_count',file=zz)
    for seq in list_of_genomes:
        genome_covbases_lengths={}
        next_read_dist=[]
        ratio_50000=[]
        contig_lengths=[]
        nocount_in_wa=0
        noread_contig_length=0
        genome_readcount=0         
        for contig in genome_contigs[seq]:
            genome_change=[0]*contigs_lengths_d[contig]
            rn=0 
            pair_taken={}
            positions=[] 
            for line in sam.fetch(contig=contig):
                rn+=1
                a=str(line).split('\t')                 
                try:
                    pair_taken[a[0]]                                              
                except KeyError:                    
                    pair_taken[a[0]]=1                                               
                    positions.append(int(a[3]))
                clip_lengths={'L':0,'R':0}
                cigar=a[5].split('S')                
                if len(cigar)>1:      # cigar: change way to manage it and include other operations      
                    for i in range(len(cigar)):
                        try:
                            clip_length=int(cigar[i])                                
                            clip_lengths['L']=clip_length                                        
                        except ValueError:
                            if i==0:
                                try:
                                    clip_lengths['L']=int(cigar[0].split('H')[-1])
                                except ValueError:
                                    clip_lengths['R']=int(cigar[0].split('M')[-1])
                            elif i>1:
                                clip_lengths['R']=int(cigar[1].split('M')[-1])                            
                start_pos=int(a[3])-1+clip_lengths['L'] 
                end_pos=start_pos+len(a[9])-clip_lengths['R']-clip_lengths['L']
                try:
                    genome_change[start_pos]+=1
                except IndexError:
                    continue
                try:                 
                    genome_change[end_pos]-=1
                except IndexError:
                    continue
            if rn==0:
                noread_contig_length+=contigs_lengths_d[contig]
                try:
                    genome_covbases_lengths[seq][0]+=contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][1]+=contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][2]+=0
                    continue
                except KeyError:
                    genome_covbases_lengths[seq]=[contigs_lengths_d[contig],contigs_lengths_d[contig],0]
                    continue
            current_coverage=0
            total_depth=0
            zerodepth=0
            for position in range(contigs_lengths_d[contig]): # iterate over the genome change list to calculate: 
                current_coverage=current_coverage+genome_change[position]
                if current_coverage==0:
                    zerodepth+=1                              # the positions with zero depth
                    continue
                total_depth+=current_coverage                 # the sum of all the per-base depths to calculate the average coverage
            
            try:
                genome_covbases_lengths[seq][0]+=zerodepth
                genome_covbases_lengths[seq][1]+=contigs_lengths_d[contig]
                genome_covbases_lengths[seq][2]+=total_depth
            except KeyError:
                genome_covbases_lengths[seq]=[zerodepth,contigs_lengths_d[contig],total_depth]
            genome_readcount+=rn

            if contig_lengths_d[contig]<min_scaffold_length: # ignoring contigs that are shorter than minimum-length value (option) for the distance metrics calculation                
                nocount_in_wa+=contig_lengths_d[contig]
                continue

            contig_lengths.append(contig_lengths_d[contig])

            if len(positions)<3: # if less than three fragments map on the contig, it will count in the wheighted average according to the penalty (option)
                next_read_dist.append(penalty)
                ratio_50000.append(penalty)
                print(contig+'\t'+str(contig_lengths_d[contig])+'\t'+str(total_depth/contig_lengths_d[contig])+'\t'+str((contig_lengths_d[contig]-zerodepth)/contig_lengths_d[contig])+'\t'+str(((contig_lengths_d[contig]-zerodepth)/contig_lengths_d[contig])/(1-np.exp(-0.883*(total_depth/contig_lengths_d[contig]))))+'\t'+'0.5'+'\t'+'0.5'+'\t'+str(rn),file=zz)        
                continue 
                                
            # Average distance between consecutive fragments            
            sum_1=0
            for i in range(len(positions)-1):
                sum_1+=positions[i+1]-positions[i]
            next_read_ratio=(sum_1/(len(positions)-1))/(contig_lengths_d[contig]/len(positions))        
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
                random_ratio=(sum/50000)/(contigs_lengths_d[contig]*(1/3)) 
            else:
                tmp=all_distances(positions)
                random_ratio=tmp/(contig_lengths_d[contig]*(1/3))   
            ratio_50000.append(random_ratio)
            print(contig+'\t'+str(contig_lengths_d[contig])+'\t'+str(total_depth/contig_lengths_d[contig])+'\t'+str((contig_lengths_d[contig]-zerodepth)/contig_lengths_d[contig])+'\t'+str(((contig_lengths_d[contig]-zerodepth)/contig_lengths_d[contig])/(1-np.exp(-0.883*(total_depth/contig_lengths_d[contig]))))+'\t'+str(next_read_ratio)+'\t'+str(random_ratio)+'\t'+str(rn),file=zz)

        if genome_readcount==0:
            continue
        if len(contig_lengths)==0: # this may happen if the contig-length minimum is higher than the length of any contig in a genome
            continue

        for i in range(len(contig_lengths)):  # Each contig value is weighted for the ratio between the contig length and the total one
            next_read_dist[i]=next_read_dist[i]*(contig_lengths[i]/(genome_covbases_lengths[seq][1]-nocount_in_wa))
            ratio_50000[i]=ratio_50000[i]*(contig_lengths[i]/(genome_covbases_lengths[seq][1]-nocount_in_wa))
        next_read_dist.append(penalty*(noread_contig_length/(genome_covbases_lengths[seq][1]-nocount_in_wa))) # contig without reads will count in the weighted average, the penalty is the value
        ratio_50000.append(penalty*(noread_contig_length/(genome_covbases_lengths[seq][1]-nocount_in_wa)))

        # calculation and printing for all the values
        coverage=genome_covbases_lengths[seq][2]/genome_covbases_lengths[seq][1]
        breadth=(genome_covbases_lengths[seq][1]-genome_covbases_lengths[seq][0])/genome_covbases_lengths[seq][1]
        beb_ratio=breadth/(1-np.exp(-0.883*coverage))
        nextratio=str(np.sum(next_read_dist)) 
        randomratio=str(np.sum(ratio_50000))
        print(seq+'\t'+str(genome_covbases_lengths[seq][1])+'\t'+str(coverage)+'\t'+str(breadth)+'\t'+str(beb_ratio)+'\t'+str(nextratio)+'\t'+str(randomratio)+'\t'+str(genome_readcount),file=z)
    z.close()
    zz.close()
    queue.put(0,block=False) # when a process ends, it puts something in a queue which has max size equal to the number of processes
    if queue.full(): # the process that finds a full queue is the last to end, it joins the tmp files into the final output and deletes them.          
        tmp_scaf=''
        tmp_gen=''        
        for i in range(processes):
            tmp_scaf+=str(i)+'tmp_multiall_metrics_scaffolds.tsv '
            tmp_gen+=str(i)+'tmp_multiall_metrics_genomes.tsv '        
        os.system("cat "+tmp_scaf+">"+output+"_scaffolds.tsv")
        os.system("cat "+tmp_gen+">"+output+"_genomes.tsv")
        os.system("rm -f "+tmp_scaf)
        os.system("rm -f "+tmp_gen)

print('reading bam file')                
if __name__ == "__main__": 
    queue = multiprocessing.Queue(maxsize=processes) # the queue is needed for understanding when an ending process is the last to end              
    multi_main(genomes,processes)
