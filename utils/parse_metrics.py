import argparse

parser=argparse.ArgumentParser(description='Parse metrics file to generate abundance file. Useful to test different metrics thresholds without parsing the bam file.')
parser.add_argument('-metrics',help='Output of metapresence.py containing the metric values (*_metrics.tsv).',required=True)
parser.add_argument('-output_abundances',help='Output file to store the abundance values.',required=True)
parser.add_argument("-ber_threshold",metavar='[float]',type=float,default=0.8,help="Breadth-Expected breadth Ratio threshold. All genomes with BER value below the threshold are considered absent. Default: 0.8")
parser.add_argument("-fug_threshold",metavar='[float]',type=float,default=0.5,help="Fraction of Unexpected Gaps threshold. All genomes with coverage lower than <max_for_fug> and with FUG value - for both mates - below the threshold are considered absent. Default: 0.5")
parser.add_argument("-min_reads",type=int,metavar='[int]',default=80,help="Number of mapped reads on a given genome below which it is considered absent. Default: 80")
parser.add_argument("-max_for_fug",type=int,metavar='[float]',default=0.1,help="Coverage value above which only BER metric is used. Default=0.1")
parser.add_argument("-fug_criterion",metavar='["all","any","mean"]',default="all",help="Write < all > if a present species must have the FUG values for both the group of mates above the threshold, < any > if only one FUG value, < mean > if the mean FUG value. Irrelevant if --unpaired is set. Default=all")
parser.add_argument("--unpaired",help='--unpaired: set this flag if the aligned reads are not paired. Default: FALSE',action="store_true")
parser.add_argument("--plot_metrics",help='--plot_metrics: set this flag if a scatterplot of the metric values has to be generated. Requires Matplotlib. Default: FALSE',action="store_true")
args=parser.parse_args()
    
# Options and parameters setting
ber_threshold=args.ber_threshold
fug_threshold=args.fug_threshold
min_reads=args.min_reads
max_for_fug=args.max_for_fug
fug_criterion=args.fug_criterion
if fug_criterion not in ['all','any','mean']: print('Error. "{}" is not a valid argument for -fug_criterion. Write either "all", "any" or "mean"'.format(fug_criterion)); exit(1)

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

f=open(args.metrics); f.readline()
present_coverage={}
for i in f:
    a=i.strip().split('\t')
    cov,ber,fug1,fug2,nreads=float(a[2]),float(a[4]),(float(a[5]) if not args.unpaired else None),float(a[6]),int(a[7])
    if parse_metrics(cov,ber,fug1,fug2,nreads,fug_criterion):
        present_coverage[a[0]]=cov

new=open(args.output_abundances,"w")
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
        exit()
    fug=[]
    ber=[]
    f=open(args.metrics)
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
    plt.savefig(args.output_abundances+".scatterplot.png")
