import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
max_lr_score = sys.argv[2]
template_path = sys.argv[3]
num_seq = sys.argv[4]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path))
template = env.get_template("render_fasta.html")

# embed gt images in html img tags
filenames = os.listdir(target_folder)
gts = sorted([i.split(".")[0] for i in filenames if 'pfm_gt' in i and 'csv' in i and "_c" not in i])
gts_indices = [i+1 for i in range(len(gts))]

def get_num_fn(fn):
    return int(fn.split("pfm")[1])

split_size_d = 100;
discovered = sorted([i.split(".")[0] for i in filenames if 'd_pfm' in i and 'csv' in i and "_c" not in i],key=get_num_fn)
d_indices = [i+1 for i in range(len(discovered))]
discovered = [[discovered[i] for i in range(len(discovered)) if i%split_size_d==j] for j in range(split_size_d)]
d_indices = [[d_indices[i] for i in range(len(d_indices)) if i%split_size_d==j] for j in range(split_size_d)]

# read msa_count and ssc
q=pd.read_csv(target_folder+"/info.csv")
msa_count = [i for i in q.msa_count]
ssc = [i for i in q.ssc]
evalues = [i for i in q.evalues]
msa_count = [[msa_count[i] for i in range(len(msa_count)) if i%split_size_d==j] for j in range(split_size_d)]
ssc = [[ssc[i] for i in range(len(ssc)) if i%split_size_d==j] for j in range(split_size_d)]
evalues = [[evalues[i] for i in range(len(evalues)) if i%split_size_d==j] for j in range(split_size_d)]

f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(gts = gts,\
                        gts_indices = gts_indices, \
						discovered = discovered, \
                        d_indices = d_indices, \
						target_folder=target_folder, \
                        max_lr_score=max_lr_score, \
                        msa_count=msa_count, \
                        ssc=ssc, \
                        evalues=evalues, \
                        num_seq=num_seq
						))
f.close()
