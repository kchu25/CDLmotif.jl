import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
perf_coeff = sys.argv[2]
max_lr_score = sys.argv[3]
prob_perseq = sys.argv[4]
template_path = sys.argv[5]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path))
template = env.get_template("render_sim.html")

# embed gt images in html img tags
filenames = os.listdir(target_folder)
# split_size = 3;
gts = sorted([i.split(".")[0] for i in filenames if 'pfm_gt' in i and 'csv' in i and "_c" not in i])
gts_indices = [i+1 for i in range(len(gts))]

def get_num_fn(fn):
    return int(fn.split("pfm")[1])

split_size_d = 5;
discovered = sorted([i.split(".")[0] for i in filenames if 'd_pfm' in i and 'csv' in i and "_c" not in i],key=get_num_fn)
d_indices = [i+1 for i in range(len(discovered))]
discovered = [[discovered[i] for i in range(len(discovered)) if i%split_size_d==j] for j in range(split_size_d)]
d_indices = [[d_indices[i] for i in range(len(d_indices)) if i%split_size_d==j] for j in range(split_size_d)]
# discovered = [discovered[split_size_d*i:split_size_d*i+split_size_d] for i in range(len(discovered)//split_size_d + 1)]
# d_indices = [d_indices[split_size_d*i:split_size_d*i+split_size_d] for i in range(len(d_indices)//split_size_d + 1)]

# read motif_type info
df = pd.read_csv(target_folder+"/motif_type.csv", sep='::', header=None, engine='python')
modes = df.shape[0]
gap_strs = [i for i in df.iloc[:,1].values]
mode_strs = [i for i in df.iloc[:,0].values]
gap_strs_split = [str(gap_strs[i]).split(',') for i in range(modes)]
mode_strs_split = [mode_strs[i].split(',') for i in range(modes)]
# mode_str_pairs = [ [ [i[j],i[j+1]] for j in range(len(i)-1)] for i in mode_strs_split if len(i) > 1]
mode_str_pairs = [ [ [i[j],i[j+1]] for j in range(len(i)-1)] if len(i) > 1 else [] for i in mode_strs_split]
# gap_str_pairs = [ [ int(float(i[j])) for j in range(len(i)) ] for i in gap_strs_split if len(i) > 1  or i[0].replace('.','').isnumeric()]
gap_str_pairs = [ [ int(float(i[j])) for j in range(len(i)) ]  if len(i) > 1  or i[0].replace('.','').isnumeric() else [] for i in gap_strs_split]
mixture_weights = [1]
if modes > 1:
    mixture_weights = [i for i in df.iloc[:,2].values]


f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(gts = gts,\
                        gts_indices = gts_indices, \
                        prob_perseq = prob_perseq, \
						discovered = discovered, \
                        d_indices = d_indices, \
						target_folder=target_folder, \
                        modes=modes,\
                        mode_strs=mode_strs,\
                        mode_str_pairs=mode_str_pairs,\
                        gap_str_pairs=gap_str_pairs,\
                        mixture_weights=mixture_weights,\
                        perf_coeff=perf_coeff,\
                        max_lr_score=max_lr_score
						))
f.close()
