__author__ = 'mpopovic'

def smoothPlot(x, y, *args, **kwargs):
  kernel = np.ones(Nsmooth)/Nsmooth
  plt.plot(np.convolve(x,kernel,'valid'), np.convolve(y,kernel,'valid'), *args, **kwargs)

import os
import sqlite3 as lite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.io.parsers as prs
import pandas.io.sql as psql

name = '111102'
#name = 'severed'
#name = 'MT_cold'
#name = 'MT_hot'
#name = 'HT_hot'
for name in ['111102', 'severed', 'MT_cold', 'MT_hot', 'HT_hot']:
    if name == '111102':
        inPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Matthias_analysis/Deformation data/movieCoordinates/wild type/111102/blade/'
    if name == 'severed':
        inPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Matthias_analysis/Deformation data/movieCoordinates/severedHinge/130107/blade/'
    if name == 'MT_cold':
        inPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Matthias_analysis/Deformation data/cdc2/homozygous/25C/01/blade/'
    if name == 'MT_hot':
        inPath = '/Network/home/PupalWingDeformation/DEcadherinGFP_wings/Deformation data/cdc2/homozygous/25C_16hAPF_30C/01/blade/'
    if name == 'HT_hot':
        inPath = '/Network/home/PupalWingDeformation/DEcadherinGFP_wings/Deformation data/cdc2/heterozygous/25C_16hAPF_30C/03/blade/'
    inFile_iso = 'triangleState.dat'
    iso_df = prs.read_csv(inPath+inFile_iso, sep = '\t')
    time = np.array(iso_df['# time;'])
    dT = time[1:]-time[:-1]
    inFile_cell_nr = 'cellNumber.dat'
    cell_nr_df = prs.read_csv(inPath+inFile_cell_nr, sep='\t')
    nr = np.array(cell_nr_df['total cell number;'])
    av_area = np.array(iso_df['sum area;'])/nr
    dnr = (nr[1:] - nr[:-1])/(time[1:]-time[:-1])
    dav_area = (av_area[1:] - av_area[:-1])/(time[1:]-time[:-1])
    area_change = dav_area/av_area[:-1]
    number_change = dnr/nr[:-1]
    inFile_nr_change = 'cellNumberChange.dat'
    nr_change_df=prs.read_csv(inPath+inFile_nr_change,sep='\t')
    change_time = nr_change_df['# time (mid interval);']
    division_change = nr_change_df['cell extrusions;']
    apoptosis_change = nr_change_df['cell appearances due to segmentation errors;']
    relative_division_change = 1.*division_change/nr[:-1]
    relative_division_change_rate = relative_division_change/dT
    relative_apoptosis_change = 1.*apoptosis_change/nr[:-1]
    relative_apoptosis_change_rate = relative_apoptosis_change/dT
    inFile_shear = 'dualMarginDeformation.dat'
    shear_df = prs.read_csv(inPath+inFile_shear, sep='\t')
    time_shear = np.array(shear_df['# time (mid interval);'])
    deltaT = time_shear[1:]-time_shear[:-1]
    shear1 = np.array(0.5*(shear_df['Delta u_{xx};']-shear_df['Delta u_{yy};']))
    shear0 = np.array((shear_df['Delta u_{xx};']+shear_df['Delta u_{yy};']))
    Nsmooth = 10
    f=plt.figure(figsize=(10/2.54,6/2.54))
    plt.subplots_adjust(top=0.95, bottom=0.21, hspace=0.05, right=0.92, left=0.17)
    smoothPlot(time[:-1], area_change+number_change, label = r'istropic velocity gradient', linewidth=2.0)
    smoothPlot(time[:-1], area_change, label = 'cell area', linewidth=2.0)
    smoothPlot(time[:-1], relative_division_change_rate, label = 'divisions', linewidth=2.0)
    smoothPlot(time[:-1], -relative_apoptosis_change_rate, label = 'apoptosis', linewidth=2.0)
    smoothPlot(time_shear[:-1], shear0[:-1]/deltaT, "y--", label='check', linewidth=2.0)
    for lbl in plt.gca().xaxis.get_majorticklabels():
        lbl.set_fontsize(10)
    for lbl in plt.gca().yaxis.get_majorticklabels():
        lbl.set_fontsize(10)
    plt.legend(bbox_to_anchor=(.55,0.97), loc=2, prop= {'size':6})
    plt.xlabel(r'time [hAPF]', fontsize = 12)
    plt.ylabel(r'rate[$\mathrm{h}^{-1}$]', fontsize = 12)
    plt.ylim(-0.25,0.25)
    plt.grid()
    f.savefig('figures/'+name+'_isotropic_deformation_contributions.svg')
    plt.show()
    f=plt.figure(figsize=(10/2.54,6/2.54))
    plt.subplots_adjust(top=0.95, bottom=0.21, hspace=0.05, right=0.92, left=0.17)
    plt.plot(time[:-1], np.cumsum((area_change+number_change)*dT), label = r'cumulative isotropic vel. gradient', linewidth=2.0)
    plt.plot(time[:-1], np.cumsum(area_change*dT), label = r'cell area', linewidth=2.0)
    # plt.plot(time[:-1], np.cumsum(number_change*dT), label = r'cell number', linewidth=2.0)
    smoothPlot(time[:-1], np.cumsum(relative_division_change_rate*dT), label = 'divisions', linewidth=2.0)
    smoothPlot(time[:-1], -1.*np.cumsum(relative_apoptosis_change_rate*dT), label = 'apoptosis', linewidth=2.0)
    plt.xlabel(r'time [hAPF]', fontsize = 12)
    plt.grid()
    for lbl in plt.gca().xaxis.get_majorticklabels():
        lbl.set_fontsize(10)
    for lbl in plt.gca().yaxis.get_majorticklabels():
        lbl.set_fontsize(10)
    plt.legend(bbox_to_anchor=(.45,0.95), loc=2, prop= {'size':6})
    plt.savefig('figures/'+name+'_cumulative_isotropic_deformation_contributions.svg')
    plt.show()
    plt.figure()
    smoothPlot(time_shear[:-1],shear_df['Delta u_{xx};'][:-1]/deltaT, label = 'xx', linewidth=2.0)
    smoothPlot(time_shear[:-1],shear_df['Delta u_{yy};'][:-1]/deltaT, label = 'yy', linewidth=2.0)
    plt.legend(loc='best')
    plt.xlabel('time[h]', fontsize=10)
    plt.ylabel('istoropic velocity gradient component', fontsize=12)
    plt.ylim(-0.07,0.07)
    plt.grid()
    plt.savefig('figures/'+name+'_deformations.png')
    plt.figure()
    plt.plot(time_shear[:-1],np.cumsum(shear_df['Delta u_{xx};'][:-1]), label = 'xx', linewidth=2.0)
    plt.plot(time_shear[:-1],np.cumsum(shear_df['Delta u_{yy};'][:-1]), label = 'yy', linewidth=2.0)
    plt.legend(loc='best')
    plt.xlabel('time[h]', fontsize=20)
    plt.ylabel('cumulative velocity gradient', fontsize=20)
    plt.grid()
    plt.savefig('figures/'+name+'_cumulative_deformations.png')


# inPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/WT_25deg_111102/'
# inFile = 'WT_25deg_111102.sqlite'
# con = lite.connect(inPath+inFile)
# df = psql.frame_query('SELECT * FROM cellinfo WHERE cell_id>10000', con)
#
# df_divided = df[df['lost_by']=='Division']
# df_apoptotic = df[df['lost_by']=='Apoptosis']
# df_divided.head()
# df_apoptotic.head()