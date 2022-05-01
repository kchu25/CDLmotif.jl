import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
matrix_id = sys.argv[2]
gt_cover = sys.argv[3]
false_cover = sys.argv[4]
gt_n_cover = sys.argv[5]
a_cover = sys.argv[6]
max_lr_score = sys.argv[7]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='./'))
template = env.get_template('jinja_templates/render_jaspar2.html')

# embed gt images in html img tags
filenames = os.listdir(target_folder)
gts = sorted([i.split(".")[0] for i in filenames if 'pfm_gt' in i and 'csv' in i and "_c" not in i])
gts_indices = [i+1 for i in range(len(gts))]

def get_num_fn(fn):
    return int(fn.split("pfm")[1])

split_size_d = 5;
discovered = sorted([i.split(".")[0] for i in filenames if 'd_pfm' in i and 'csv' in i and "_c" not in i],key=get_num_fn)
d_indices = [i+1 for i in range(len(discovered))]
discovered = [[discovered[i] for i in range(len(discovered)) if i%split_size_d==j] for j in range(split_size_d)]
d_indices = [[d_indices[i] for i in range(len(d_indices)) if i%split_size_d==j] for j in range(split_size_d)]

# read motif_type info

f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(gts = gts,\
                        gts_indices = gts_indices, \
						discovered = discovered, \
                        d_indices = d_indices, \
						target_folder=target_folder, \
                        jaspars=matrix_id,\
                        gt_cover=gt_cover, \
                        false_cover=false_cover, \
                        gt_n_cover=gt_n_cover, \
                        a_cover=a_cover, \
                        max_lr_score=max_lr_score
						))
f.close()
