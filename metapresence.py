import argparse
import os
try:
    import numpy as np
except ImportError:
    print('Error: numpy is not installed.')
    exit(1)
try:                              
    import pysam
except ImportError:
    print('Error: pysam is not installed.')
    exit(1)
import multiprocessing

parser=argparse.ArgumentParser()
parser.add_argument("input_fasta",help="Either a folder of fasta files containing the sequences to analyze, or a (multi)fasta file containing the sequences to analyze. For large reference sizes, an input folder and more than one process are recommended for speeding-up.")
parser.add_argument("indexed_sorted_bam",help="Indexed sorted bam file, with bam and index (.bai) in the same directory")
parser.add_argument("-input_bins",metavar='input_bins.txt',help="A text file with each line listing a contig and a bin name, tab-seperated. This file is needed only if <input_fasta> is a fasta file and not a folder of fasta files.")
parser.add_argument("-p",type=int,metavar='processes',default=1, help="number of processes to be executed in parallel. Default: 1")
parser.add_argument("-o",metavar='output_prefix',default='metapout',help="the prefix of the output files. Default: metapout")
parser.add_argument("-ber_threshold",metavar='[float]',type=float,default=0.8,help="Breadth-Expected breadth Ratio threshold. All genomes with BER value below the threshold are considered absent. Default: 0.8")
parser.add_argument("-fug_threshold",metavar='[float]',type=float,default=0.5,help="Fraction of Unexpected Gaps threshold. All genomes with coverage lower than <max_for_fug> and with FUG value - for both mates - below the threshold are considered absent. Default: 0.5")
parser.add_argument("-min_reads",type=int,metavar='[int]',default=80,help="Number of mapped reads on a given genome below which it is considered absent. Default: 80")
parser.add_argument("-max_for_fug",type=int,metavar='[float]',default=0.1,help="Coverage value above which only BER metric is used. Default=0.1")
parser.add_argument("-fug_criterion",metavar='["all","any","mean"]',default="all",help="Write < all > if a present species must have the FUG values for both the group of mates above the threshold, < any > if only one FUG value, < mean > if the mean FUG value. Irrelevant if --unpaired is set. Default=all")
parser.add_argument("--unpaired",help='--unpaired: set this flag if the aligned reads are not paired. Default: FALSE',action="store_true")
parser.add_argument("--all_contigs",help='--all_contigs: set this flag if each contiguous sequence in < input_fasta > should be evaluated independently and be reported in the output files. If this flag is set < -input_bins > is irrelevant. Default: FALSE',action="store_true")
parser.add_argument("--plot_metrics",help='--plot_metrics: set this flag if a scatterplot of the metric values has to be generated. Requires Matplotlib. Default: FALSE',action="store_true")
parser.add_argument("--quiet",help='--quiet: set this flag to silence warning messages. default: FALSE',action="store_true")
parser.add_argument('--version', '-v', '-version', action='version', version='metapresence 1.0')

args=parser.parse_args()

if not os.path.isdir(args.input_fasta) and not args.input_bins and not args.all_contigs:
    print('\nError: if < input_fasta > is not a folder of fasta files, < -input_bins > must be specified or < --all_contigs > must be set. \n')
    exit(1)
    
# Options and parameters setting
processes=args.p
output=args.o
scaffold_input=args.input_fasta
sorted_bam=args.indexed_sorted_bam
inputbin=args.input_bins 
ber_threshold=args.ber_threshold
fug_threshold=args.fug_threshold
min_reads=args.min_reads
max_for_fug=args.max_for_fug
fug_criterion=args.fug_criterion
if fug_criterion not in ['all','any','mean']: print('Error. "{}" is not a valid argument for -fug_criterion. Write either "all", "any" or "mean"'.format(fug_criterion)); exit(1)

# The function to divide the lists of genomes in n sublists, with n as the number of processes. 
def chunks(l,nproc):
    right_div=len(l)//nproc
    nmore=len(l)%nproc
    return [l[i*(right_div+1):i*(right_div+1)+right_div+1] for i in range(nmore)]+[l[nmore*(right_div+1)+i*right_div:nmore*(right_div+1)+i*right_div+right_div] for i in range(nproc-nmore) if nmore<len(l)]

def parse_metrics(cov,ber,fug1,fug2,nreads,criterion):
    if nreads<min_reads:
        return
    if cov>max_for_fug:
        if ber>=ber_threshold:
            return True
    else:
        if criterion=='all':
            if ber>=ber_threshold and fug1>=fug_threshold and (fug2>=fug_threshold if not args.unpaired else True):
                return True
        elif criterion=='any':
            if ber>=ber_threshold and fug1>=fug_threshold or (fug2>=fug_threshold if not args.unpaired else False):
                return True
        elif criterion=='mean':
            if ber>=ber_threshold and ((fug1+fug2)/2>=fug_threshold if not args.unpaired else fug1>=fug_threshold):
                return True

def parse_SingleFasta(fastafile,scaffold_input,all_contigs):
    genome_contigs={}
    if not all_contigs:
        f=open(scaffold_input)
        for line in f:
            a=line.strip().split('\t')
            if a[1] in genome_contigs:
                genome_contigs[a[1]].append(a[0])
            else:
                genome_contigs[a[1]]=[a[0]]
        f.close()
    contigs_lengths_d={}
    f=open(fastafile)
    contig_name=''
    for line in f:
        if line[0]=='>':
            header=line[1:].strip().split()
            contig_name=header[0]
            if all_contigs: genome_contigs[contig_name]=[contig_name]
            contigs_lengths_d[contig_name]=0
        else:
            contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()
    return contigs_lengths_d,genome_contigs

def pool_parse_MultiFasta(args):
    ff,all_contigs=args
    contigs_lengths_d,genome_contigs={},{}
    for fa in ff:
        f=open(fa)
        fa=fa.split('/')[-1]
        for line in f:
            if line[0]=='>':
                contig_name=line.strip().split()[0][1:]
                contigs_lengths_d[contig_name]=0
                if not all_contigs:
                    if fa in genome_contigs:
                        genome_contigs[fa].append(contig_name)
                    else:
                        genome_contigs[fa]=[contig_name]
                else: genome_contigs[contig_name]=[contig_name]
            else:
                contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()
    return contigs_lengths_d,genome_contigs

# parse bam file, calculate metrics for each genome
def parseBam(list_of_genomes):
    printing_dict={}
    #list_of_genomes,bam,contigs_lengths_d,genome_contigs,arl=these_arguments
    sam=pysam.AlignmentFile(sorted_bam.split()[0],"rb")    
    for seq in list_of_genomes:
        genome_covbases_lengths={}
        noread_contig_length=0
        genome_readcount=0
        genome_pointer=0
        allpositions=[0]
        allpositions2=[0]
        # retrieve mapping positions from bam file        
        for contig in genome_contigs[seq]:
            if contig not in contigs_lengths_d:
                print('Warning: contig {} of sequence {} was not present in the fasta file(s). Skipping sequence...'.format(contig,seq,seq)) if not args.quiet else None
                genome_readcount=0
                break
            genome_change=[0]*contigs_lengths_d[contig] 
            rn=0 
            pair_taken=set()
            try:
                for read in sam.fetch(contig=contig):
                    rn+=1
                    start_pos=read.reference_start
                    end_pos=start_pos+read.query_length               
                    if read.query_name not in pair_taken:
                        allpositions.append(start_pos+genome_pointer)
                        pair_taken.add(read.query_name)                                              
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
            genome_pointer+=contigs_lengths_d[contig] 
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
            # calculate contig depth
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
        # skip genomes with zero mapped reads    
        if genome_readcount==0:
            continue
        # add fake read at the end of the genome
        allpositions.append(genome_covbases_lengths[seq][1]-avreadlen) # average read length
        allpositions2.append(genome_covbases_lengths[seq][1]-avreadlen) 
        # calculation and printing for all the values
        window_size=round(genome_covbases_lengths[seq][1]/(len(allpositions)))
        if window_size>0:
            distances=np.array(allpositions[1:])-np.array(allpositions[:-1])
            hist=np.histogram(distances,bins=range(window_size,np.max(distances)+1))[0]
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
        
        printing_dict[seq]=str(genome_covbases_lengths[seq][1])+'\t'+str(coverage)+'\t'+str(breadth)+'\t'+str(beb_ratio)+'\t'+str(gennext_read_ratio)+'\t'+str(gennext_read_ratio2)+'\t'+str(genome_readcount)
    
    return printing_dict
     
if __name__ == "__main__":  

    if os.path.isdir(scaffold_input):
        print('reading fasta files')
        fastas=[scaffold_input+'/'+x for x in os.listdir(scaffold_input)]
        chunked_fastas=chunks(fastas,processes)
        parseFasta_args=[(x,args.all_contigs) for x in chunked_fastas]
        with multiprocessing.Pool(processes=processes) as pool:
            all_parsing=pool.map(pool_parse_MultiFasta,parseFasta_args)
        contigs_lengths_d,genome_contigs={x:result[0][x] for result in all_parsing for x in result[0]},{x:result[1][x] for result in all_parsing for x in result[1]}
        if args.all_contigs: genome_contigs={y:[y] for x in genome_contigs for y in genome_contigs[x]}
    elif os.path.isfile(scaffold_input):
        print('reading fasta file')
        contigs_lengths_d,genome_contigs=parse_SingleFasta(scaffold_input,inputbin,args.all_contigs)
    else:
        print('Error: cannot open',scaffold_input,'\naborted')
        exit(1)
    print('\tdone\n') 
    print('reading bam file')
    try:
        save=pysam.set_verbosity(0) # to avoid missing index error coming up
        sam=pysam.AlignmentFile(sorted_bam.split()[0],"rb")
        pysam.set_verbosity(save)
        totlen=0
        c=0
        for read in sam.fetch():
            totlen+=read.query_length
            c+=1
            if c==1000:
                break
        avreadlen=round(totlen/c)
        del sam
    except ValueError:
        print('Error: cannot find index for '+sorted_bam+'. Index and bam file should be placed in the same directory.\nAborted')
        exit(1)
    
    genomes=[x for x in genome_contigs]
    arguments=chunks(genomes,processes) 
    with multiprocessing.Pool(processes=processes) as pool:
        all_metrics=pool.map(parseBam,arguments)

    print('\tdone\n')
    print('preparing output files\n')
    # printing in tsv files
    tord={x:result[x] for result in all_metrics for x in result}
    tord={x:tord[x] for x in sorted(tord.keys())}
    new=open(output+"_metrics.tsv","w")
    present_coverage={}
    print('genome\tlength\tcoverage\tbreadth\tBER\tFUG1\tFUG2\tread_count',file=new)
    for i in tord:
        a=tord[i].strip().split('\t')
        cov,ber,fug1,fug2,nreads=float(a[1]),float(a[3]),float(a[4]),(float(a[5]) if not args.unpaired else None),int(a[6])
        if parse_metrics(cov,ber,fug1,fug2,nreads,fug_criterion):
            present_coverage[i]=cov
        print(i+'\t'+tord[i],file=new)
    new.close()
    new=open(output+"_abundances.tsv","w")
    print('genome\trelative_abundance_%',file=new)
    totcov=sum(present_coverage.values())
    for i in present_coverage: 
        print(i+'\t'+str(present_coverage[i]/totcov*100),file=new)
    new.close()
    # plotting
    if args.plot_metrics:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('Warning: not able to generate the plot. Matplotlib is not installed.\n')
            print('All done!\n')
            exit()
        fug=[]
        ber=[]
        f=open(output+"_metrics.tsv")
        for line in f:
            a=line.strip().split('\t')
            if a[0]!='genome':
                if float(a[-1])<min_reads:
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
        plt.savefig(output+"_scatterplot.png")
print('All done!\n')