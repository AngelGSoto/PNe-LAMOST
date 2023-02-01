import os
import shutil
import glob
from astropy.table import Table

file_list = glob.glob("better-fitModel/*.ecsv")

dir_model =  "Model-PNNew/models/N2242_He091O382Ar624_itere_R166R1705_final"
template = "{:s}"
model_name_ = []
for file_name in file_list:
    tab = Table.read(file_name)
    model_name = tab["Name model"]
    file_name_model = str(model_name).split("\n")[-1]
    all_files = os.path.join(dir_model, file_name_model) + '*'
    all_files_glob  = glob.glob(all_files)
    for i in all_files_glob:
        shutil.copy2(i, "better-fitModel")
    



