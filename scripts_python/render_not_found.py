import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
prob_perseq = sys.argv[2]
template_path = sys.argv[3]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path))
template = env.get_template("render_sim_not_found.html")

# embed gt images in html img tags
filenames = os.listdir(target_folder)
# split_size = 3;
gts = sorted([i.split(".")[0] for i in filenames if 'pfm_gt' in i and 'csv' in i and "_c" not in i])
gts_indices = [i+1 for i in range(len(gts))]

def get_num_fn(fn):
    return int(fn.split("pfm")[1])

# read motif_type info

f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(gts = gts,\
                        gts_indices = gts_indices, \
                        prob_perseq = prob_perseq, \
						target_folder=target_folder
						))
f.close()
