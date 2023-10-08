'''
Reddenered the spectra from CLOUDY 
using dust_extinction package
https://dust-extinction.readthedocs.io/en/latest/
'''
from __future__ import print_function
import numpy as np
import glob
import json
import argparse
import matplotlib.pyplot as plt
from astropy.table import Table, QTable
import seaborn as sns
import sys
import astropy.units as u
import pyCloudy as pc
from dust_extinction.parameter_averages import F99


parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("source", type=str,
                    default="model_100000_36.58",
                    help="Name of input model ")


cmd_args = parser.parse_args()
file_ = cmd_args.source 


model_name = file_
# Reading the Cloudy outputs in the Mod CloudyModel object
Mod = pc.CloudyModel(model_name)

#Getting wl and flux
wl = Mod.get_cont_x(unit='Ang')
flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')

#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))
data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Wl"] < 9000)
data_mask = data[mask]
wav = data_mask["Wl"]*u.AA
# define the model
ext = F99(Rv=3.1)

# extinguish (redden) the spectrum
#Iterating over the desired E(B-V) range and saving the reddened spectra
E = [0.0, 0.01, 0.02, 0.03, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2]

for i in range(len(E)):
    spectrum_ext = data_mask["Flux"]*ext.extinguish(wav, Ebv=E[i])
    #Save the new files
    asciifile = "redden/" + model_name + "_{}.dat".format(str(E[i])) 
    file=open(asciifile,'w') #create file  
    for x,y in zip(data_mask["Wl"], spectrum_ext):
        file.write('%f %s\n'%(x,y))     #assume you separate columns by tabs  
    file.close()     #close file 
      
