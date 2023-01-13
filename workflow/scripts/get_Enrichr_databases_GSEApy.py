#!/bin/env python
# import gseapy as gp
# import pickle
# import os
    
# # parameters from snakemake
# results_path = snakemake.output['result_file']

# # get enrichr database names
# databases = snakemake.params["databases"] #gp.get_library_name()

# # download all databases to local dict
# db_dict = dict()
# for db in databases:
#     db_dict[db]=gp.parser.gsea_gmt_parser(db, min_size=0, max_size=100000,)

# # save as pickle 
# with open(results_path, 'wb') as f:
#         pickle.dump(db_dict, f, pickle.HIGHEST_PROTOCOL)
        

#!/bin/env python       
import json
import gseapy as gp
import os
    
# config

# output
results_path = snakemake.output['result_file']

# parameters
db = snakemake.params["database"] #gp.get_library_name() #-> gets all available names
organism = "Human" if (snakemake.config["genome"]=='hg19' or snakemake.config["genome"]=='hg38' ) else "Mouse"

# download the database to local dict
# db_dict = gp.parser.gsea_gmt_parser(db, min_size=0, max_size=100000,)
db_dict = gp.parser.download_library(name=db, organism=organism)

# save as json
with open(results_path, 'w') as fp:
    json.dump(db_dict, fp,  indent=4)