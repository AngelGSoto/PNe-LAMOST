'''
Reddenered the spectra from CLOUDY 
'''
from __future__ import print_function
import numpy as np
import glob
import json
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pyCloudy as pc

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

curve="../F04_CURVE_3.1.dat"
R=3.1

#Getting wl and flux
wl = Mod.get_cont_x(unit='Ang')
flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')

#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))

#Readding the extiction curve of 
ext=np.loadtxt(curve)

wl_ext  = 10000/ext[:,0] # in angstroms
ext_val= 1+ext[:,1]/R  # in A_l/A_V
ord=np.argsort(wl_ext)
wl_ext = wl_ext.take(ord)
ext_val = ext_val.take(ord)

#Interpolating the extintion curve to the lambda of the spectrum
ext_val_inter_lamspectrum = np.interp(wll, wl_ext, ext_val)
#Iterating over the desired E(B-V) range and saving the reddened spectra
E = [0.0, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

for i in range(len(E)):
    av=R*E[i]
    cor=ext_val_inter_lamspectrum*av
    #Fluxes are being multiplied by "0.12" to correct for the errors in the convertion
    specflu_red=np.array(flux)*0.12*10**(-0.4*cor)
    
    #Save the new files
    asciifile = "redden/" + model_name + "_{}.dat".format(str(E[i])) 
    file=open(asciifile,'w') #create file  
    for x,y in zip(wll, specflu_red):
        file.write('%f %s\n'%(x,y))     #assume you separate columns by tabs  
    file.close()     #close file   



