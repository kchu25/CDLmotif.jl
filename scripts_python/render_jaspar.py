import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
template_path = sys.argv[2]

tot_num_seq_data = sys.argv[3]
max_occ_num = sys.argv[4]
max_rows = sys.argv[5]
used_sequences = sys.argv[6:]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path))
template = env.get_template('render_jaspar.html')


def get_num_fn(fn):
    return int(fn.split("pfm")[1])

filenames = os.listdir(target_folder)

split_size_d = 4;
discovered = sorted([i.split(".")[0] for i in filenames if 'd_pfm' in i and 'csv' in i and "_c" not in i], key=get_num_fn)
d_indices = [i+1 for i in range(len(discovered))]
discovered = [discovered[split_size_d*i:split_size_d*i+split_size_d] for i in range(len(discovered)//split_size_d + 1)]
d_indices = [d_indices[split_size_d*i:split_size_d*i+split_size_d] for i in range(len(d_indices)//split_size_d + 1)]

# read motif_type info

f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(discovered = discovered, \
                        d_indices = d_indices, \
						target_folder=target_folder, \
                        used_sequences=used_sequences, \
                        tot_num_seq_data=tot_num_seq_data,\
                        max_occ_num=max_occ_num, \
                        max_rows=max_rows
						))
f.close()
