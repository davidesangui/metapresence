import argparse
import os

parser=argparse.ArgumentParser(description='Merge together multiple abundance outputs from metapresence.py into a single table')
parser.add_argument("-abundance_folder",help="Folder containing abundance outputs from metapresence.py",required=True,metavar='')
parser.add_argument("-output_table",help="Output file to store merged abundances",required=True,metavar='')

args=parser.parse_args()

files=[x for x in sorted(os.listdir(args.abundance_folder))]
printing_dict={}
for n,file in enumerate(files):
    f=open(args.abundance_folder+'/'+file); f.readline()
    for line in f:
        a=line.strip().split('\t')
        if a[0] in printing_dict:
            printing_dict[a[0]].append(a[1])
        else:
            printing_dict[a[0]]=['0.0']*n
            printing_dict[a[0]].append(a[1])
    f.close()
    for gen in printing_dict: 
        if len(printing_dict[gen])<n+1:
            printing_dict[gen].append('0.0')
printing_dict={x:printing_dict[x] for x in sorted(printing_dict.keys())}
new=open(args.output_table,'w')
print('fasta\t'+'\t'.join(files),file=new)
for gen in printing_dict:
    print(gen+'\t'+'\t'.join(printing_dict[gen]),file=new)
new.close()

