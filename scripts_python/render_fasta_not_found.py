import jinja2, os, sys
import pandas as pd
target_folder = sys.argv[1]
template_path = sys.argv[2]

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path))
template = env.get_template("render_fasta_not_found.html")

# embed gt images in html img tags
filenames = os.listdir(target_folder)



f = open(target_folder+"/"+"summary.html",'w')
f.write(template.render(                    
						target_folder=target_folder, \
						))
f.close()
