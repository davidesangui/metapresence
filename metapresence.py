import argparse


parser=argparse.ArgumentParser()
parser.add_argument("contigs_genomes",help="scaffolds_genomes: either a text file with each line listing a contig and a bin name, tab-seperated, or the path to a folder containing (only) the genomes to analyze")
parser.add_argument("indexed_sorted_bam",help="indexed sorted bam file, with bam and index (.bai) in the same directory")
parser.add_argument("all_sequences_fasta",help="the fasta file on which the reads are aligned")
parser.add_argument("-p",type=int,metavar='processes',default=1, help="number of processes to be executed in parallel. Default: 1")
parser.add_argument("-o",metavar='output_prefix',default='metapout',help="the prefix of the output files. Default: metapout")
parser.add_argument("--unpaired",help='set this flag if the aligned reads are not paired. Default: FALSE',action="store_true")
parser.add_argument("--plot_metrics",help='set this flag if a scatterplot of the metric values has to be generated. Requires Matplotlib. Default: FALSE',action="store_true")

args=parser.parse_args()


import os
import random
import numpy as np                              
import pysam
import multiprocessing
    

# Options and parameters setting
processes=args.p
output=args.o
scaffold_input=args.contigs_genomes
sorted_bam=args.indexed_sorted_bam
fastafile=args.all_sequences_fasta

# The function to divide the lists of genomes in n sublists, with n as the number of processes. 
def chunks(l, processi):
    lists=[[] for x in range(processi)]
    ind=0
    for ob in l:
        lists[ind].append(ob)
        ind+=1
        if ind==processi:
            ind=0
    end=len(lists)
    for i in range(len(lists)):
        if len(lists[i])==0:
            end=i
            break
    return(lists[:end])


#contig lengths calculation into a dictionary
contigs_lengths_d={}
if not os.path.isdir(scaffold_input):
    print('reading fasta file')
    f=open(fastafile)
    contig_name=''
    for line in f:
        if line[0]=='>':
            header=line[1:].strip().split()
            contig_name=header[0]
            contigs_lengths_d[contig_name]=0
        else:
            contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()

# generation of genome_listofcontig dictionary from input.txt or folder (and contig lengths calculation into dict)
genome_contigs={}
if os.path.isdir(scaffold_input):
    print('reading fasta files')
    folder=scaffold_input[:-1] if scaffold_input[-1]=='/' else scaffold_input
    for fa in os.listdir(scaffold_input):
        f=open(folder+'/'+fa)
        for line in f:
            if line[0]=='>':
                contig_name=line.strip().split()[0][1:]
                contigs_lengths_d[contig_name]=0
                if fa in genome_contigs:
                    genome_contigs[fa].append(contig_name)
                else:
                    genome_contigs[fa]=[contig_name]
            else:
                contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()
else:
    if os.path.isfile(scaffold_input):
        f=open(scaffold_input)
        for line in f:
            a=line.strip().split('\t')
            try:
                genome_contigs[a[1]].append(a[0])
            except KeyError:
                genome_contigs[a[1]]=[a[0]]
        f.close()
    else:
        print('Error: cannot open',scaffold_input)
        print('Aborted')
        exit()

# Main list of genomes that will be spit in sublists for each process
genomes=[x for x in genome_contigs]
print('\tdone\n') 
print('reading bam file')
try:
    save = pysam.set_verbosity(0) # to avoid missing index error coming up
    sam=pysam.AlignmentFile(sorted_bam.split()[0],"rb")
    pysam.set_verbosity(save) # to make it normal back again
    totlen=0
    c=0
    for line in sam.fetch():
        a=str(line).split('\t')
        totlen+=len(a[9])
        c+=1
        if c==1000:
            break
    avreadlen=round(totlen/c)
    del sam
except ValueError:
    print('Error: cannot find index for '+sorted_bam+'. Index and bam file should be placed in the same directory.\nAborted')
    exit()


def main(list_of_genomes,order,bam,contig_lengths_d,genome_contigs,arl):
    sam=pysam.AlignmentFile(bam.split()[0],"rb")    
    z=open(str(order)+'tmp_multiall_metrics_genomes.tsv','w')
    for seq in list_of_genomes:
        genome_covbases_lengths={}
        noread_contig_length=0
        genome_readcount=0
        genome_pointer=0
        allpositions=[0]
        allpositions2=[0]        
        for contig in genome_contigs[seq]:
            genome_change=[0]*contigs_lengths_d[contig]
            rn=0 
            pair_taken={}
            try:
                for read in sam.fetch(contig=contig):
                    rn+=1
                    start_pos=read.reference_start
                    end_pos=start_pos+read.query_length               
                    if read.query_name not in pair_taken:
                        allpositions.append(start_pos+genome_pointer)
                        pair_taken[read.query_name]=1                                              
                    else:                                                                   
                        allpositions2.append(start_pos+genome_pointer)                            
                    if start_pos<len(genome_change):
                        genome_change[start_pos]+=1
                    else:
                        continue
                    if end_pos<len(genome_change):
                        genome_change[end_pos]-=1
                    else:
                        continue
            except ValueError:
                pass
            genome_pointer+=contig_lengths_d[contig] #####################
            
            ###################

            if rn==0:
                noread_contig_length+=contigs_lengths_d[contig]
                if seq in genome_covbases_lengths:
                    genome_covbases_lengths[seq][0]+=contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][1]+=contigs_lengths_d[contig]
                    genome_covbases_lengths[seq][2]+=0
                    continue
                else:
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
            
            if seq in genome_covbases_lengths:
                genome_covbases_lengths[seq][0]+=zerodepth
                genome_covbases_lengths[seq][1]+=contigs_lengths_d[contig]
                genome_covbases_lengths[seq][2]+=total_depth
            else:
                genome_covbases_lengths[seq]=[zerodepth,contigs_lengths_d[contig],total_depth]
            genome_readcount+=rn           
            
        if genome_readcount==0:
            continue

        allpositions.append(genome_covbases_lengths[seq][1]-arl) # average read length
        allpositions2.append(genome_covbases_lengths[seq][1]-arl) 

        

        # calculation and printing for all the values
        window_size=round(genome_covbases_lengths[seq][1]/(len(allpositions)))
        if window_size>0:
            distances=np.array(allpositions[1:])-np.array(allpositions[:-1])
            hist,edges=np.histogram(distances,bins=range(window_size,np.max(distances)+1))
            frequency=hist/(len(allpositions)-1)
            gennext_read_ratio=(window_size-(np.sum(frequency*np.arange(len(frequency)))))/window_size
        else:
            gennext_read_ratio='nan'
        if not args.unpaired:
            window_size=round(genome_covbases_lengths[seq][1]/(len(allpositions2)))
            if window_size>0:
                distances=np.array(allpositions2[1:])-np.array(allpositions2[:-1])
                hist,edges=np.histogram(distances,bins=range(window_size,np.max(distances)+1))
                frequency=hist/(len(allpositions2)-1)
                gennext_read_ratio2=(window_size-(np.sum(frequency*np.arange(len(frequency)))))/window_size
            else:
                gennext_read_ratio2='nan'
        else:
            gennext_read_ratio2='unpaired'

        coverage=genome_covbases_lengths[seq][2]/genome_covbases_lengths[seq][1]
        breadth=(genome_covbases_lengths[seq][1]-genome_covbases_lengths[seq][0])/genome_covbases_lengths[seq][1]
        beb_ratio=breadth/(1-np.exp(-0.883*coverage))
        
        print(seq+'\t'+str(genome_covbases_lengths[seq][1])+'\t'+str(coverage)+'\t'+str(breadth)+'\t'+str(beb_ratio)+'\t'+str(gennext_read_ratio)+'\t'+str(gennext_read_ratio2)+'\t'+str(genome_readcount),file=z)
    
    z.close()

                
if __name__ == "__main__":              
    blocchi=chunks(genomes,processes) 
    jobs=[]
    for i in range(len(blocchi)):        
        j=multiprocessing.Process(target=main,args=(blocchi[i],i,sorted_bam,contigs_lengths_d,genome_contigs,avreadlen))
        jobs.append(j)
    for j in jobs:        
        j.start()
    for j in jobs:
        j.join()

    print('\tdone\n')
    print('preparing output files')
    # printing in final tsv file
    tord={}    
    for i in range(len(blocchi)):
        f=open(str(i)+'tmp_multiall_metrics_genomes.tsv')
        for line in f:
            tord[line.split('\t')[0]]=line.strip()
        f.close()
        os.remove(str(i)+'tmp_multiall_metrics_genomes.tsv')
    tord={x:tord[x] for x in sorted(tord.keys())}
    new=open(output+"_metrics.tsv","w")
    print('genome\tlength\tcoverage\tbreadth\tBER\tFUG1\tFUG2\tread_count',file=new)
    for i in tord:
        print(tord[i],file=new)
    new.close()
    if args.plot_metrics:
        # plotting
        import matplotlib.pyplot as plt
        fug=[]
        ber=[]
        f=open(output+"_metrics.tsv")
        for line in f:
            a=line.strip().split('\t')
            if a[0]!='genome':
                if float(a[-1])<80:
                    continue
                ber.append(float(a[4]))
                if not args.unpaired:
                    fug.append((float(a[5])+float(a[6]))/2)
                else:
                    fug.append(float(a[5]))
        f.close()
        plt.scatter(ber,fug)
        plt.xlabel('BER')
        if args.unpaired:
            plt.ylabel('FUG')
        else:
            plt.ylabel('mean FUG')
        plt.savefig(output+"_scatter.png")
    

    
# da fare: error management, lunghezza contig se folder