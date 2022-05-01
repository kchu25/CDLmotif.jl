import pandas as pd, seqlogo as sq, numpy as np, os, sys

target_folder = sys.argv[1]

def get_complement(matrix): # matrix should be in pd.DataFrame format
    rev = matrix.copy()
    rev = rev.iloc[::-1]
    a_cp = rev['A'].copy()
    c_cp = rev['C'].copy()
    rev['C'] = rev['G']
    rev['G'] = c_cp
    rev['A'] = rev['T']
    rev['T'] = a_cp
    rev = rev.reset_index(drop=True)
    return rev

filenames = os.listdir(target_folder)
gts = [i.split(".")[0] for i in filenames if 'pfm_gt' in i and 'csv' in i and "_c" not in i]
discovered = [i.split(".")[0] for i in filenames if 'd_pfm' in i and 'csv' in i and "_c" not in i]

for g in gts:
	read = pd.read_csv(target_folder+"/"+g+".csv", header=None)
	ppm = sq.Ppm(read.values)
	sq.seqlogo(ppm, ic_scale = True, \
		format = 'png', size = 'medium', \
		filename=(target_folder+"/"+g.split('.')[0]+'.png'))
	sq.seqlogo(sq.Ppm(get_complement(ppm.ppm)), ic_scale = True, \
		format = 'png', size = 'medium', \
		filename=(target_folder+"/"+g.split('.')[0]+'_c.png'))

for d in discovered:
	read = pd.read_csv(target_folder+"/"+d+".csv", header=None)
	ppm = sq.Ppm(read.values)
	sq.seqlogo(ppm, ic_scale = True, \
		format = 'png', size = 'medium', \
		filename=(target_folder+"/"+d.split('.')[0]+'.png'))
	sq.seqlogo(sq.Ppm(get_complement(ppm.ppm)), ic_scale = True, \
		format = 'png', size = 'medium', \
		filename=(target_folder+"/"+d.split('.')[0]+'_c.png'))

   