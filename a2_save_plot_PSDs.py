# Plot power-spectral densities

%load_ext autoreload
%autoreload
from setup_parameters import *
import matplotlib.pyplot as plt
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as np
import os
from pathlib import Path
import pandas as pd

from noisemodels import accNLNM,accNHNM

# %% codecell

plotunit = 0 # 0: displacement, 1: velocity, 2: acceleration

path2bin = './bin/' # Set path to your environment executables

webservice = "IRIS"

xtype = 'period' # period or frequency
sw_width = 0.25 # 1/4 Smoothing window width in octave
sw_shift = 0.125 # 1/8 Smoothing window shift in fraction of octave

# ## CEDAR CITY, UT
# # Define station parameters
# network = 'TA'
# sta = 'S14A'
# loc = '--'
# compstr = "BHZ"
# tstart = "2008-08-14T12:00:00"
# tend = "2008-09-14T12:00:00"


# %% codecell
# LOAD CLIENT
client = Client(webservice)
print(client)
# LOAD STATIONS
t1 = UTCDateTime(tstart)
t2 = UTCDateTime(tend)
inventory = client.get_stations(network=network, station=sta, channel=compstr, starttime=t1, endtime=t2)
print(inventory)
stadir = './stations/'
if not os.path.exists(stadir):
    os.makedirs(stadir)
file = open(stadir+network+'.'+sta+'.'+loc+'.txt', 'w')
for ista in range(0,len(inventory[0])) :
    file.write("%5s %5s %12f %12f %12f\n" % (
                                        inventory[0].code,
                                        inventory[0].stations[ista]._code, 
                                        inventory[0].stations[ista]._latitude, 
                                        inventory[0].stations[ista]._longitude, 
                                        inventory[0].stations[ista]._elevation))
file.close()

# %% codecell
# Run PSD computation

# !{'python3 '+path2bin+'ntk_computePSD.py'} net={network} sta={sta} loc={loc} chan={compstr} start={tstart} end={tend} xtype={xtype} sw_width={sw_width} sw_shift={sw_shift} plot=0

# %% codecell
# Gather PSDs

path2psdDb = './data/psdDb/'+network+'.'+sta+'.'+loc+'/'+compstr+'/'

pathlist = sorted(Path(path2psdDb).glob('*/*/*.'+xtype+'.txt'))

list_of_dfs = []
for path in pathlist:
    df = pd.read_table(path, delim_whitespace=True,skiprows=[0],names=[xtype,'psd'])
    # sta = data.sta[0]
    df = df.set_index(xtype)

    freqs = 1/df.index.values

    df = 10**(df/10) # to non-log units of m**2/s**4/Hz
    df = 10*np.log10(df['psd'] / (2*np.pi*freqs)**(2-plotunit)) # convert to dB rel 1 m**2/s**4/Hz
    df = 10**(df/20) # convert dB rel 1 m**2/s**4/Hz --> m/Hz**0.5 (not log)

    list_of_dfs.append(df)

# Concatinate psd dataframe
dfs = pd.concat(list_of_dfs, ignore_index=True,axis=1)

# Drop nan values
dfs = dfs.dropna()

# %%
# Convert everything to m/sqrt(Hz)

dfs.index = 1/dfs.index
freqs = dfs.index.values

nlnm = 10*np.log10(accNLNM(freqs) / (2*np.pi*freqs)**(2-plotunit))
nlnm_m_hz = 10**(nlnm/20) # convert dB rel 1 m**2/s**4/Hz --> m/Hz**0.5 (not log)

nhnm = 10*np.log10(accNHNM(freqs) / (2*np.pi*freqs)**(2-plotunit))
nhnm_m_hz = 10**(nhnm/20) # convert dB rel 1 m**2/s**4/Hz --> m/Hz**0.5 (not log)

# %%
# Plot PSDs

plt.figure(figsize=(10,6))
plt.loglog(dfs, color='gray', alpha=0.5)
# Plot median
plt.loglog(dfs.median(axis=1), color='red', alpha=1, linewidth=2, label='Median')
# Plot 97.5% quantile
plt.loglog(dfs.quantile(q=0.975,axis=1), linestyle='--', color='red', alpha=1, linewidth=2, label='97.5')
# Plot 2.5% quantile
plt.loglog(dfs.quantile(q=0.025,axis=1), linestyle='--', color='red', alpha=1, linewidth=2, label='2.5')
# Plot noise models
plt.loglog(freqs, nlnm_m_hz,
              c='blue', label='NLNM')
plt.loglog(freqs, nhnm_m_hz,
            c='blue', label='NHNM')
plt.xlim([0.005,20])
plt.ylim([1e-11,1e-3])
plt.xlabel('Frequency (Hz)')
plt.ylabel(r'$m/\sqrt{Hz}$')
plt.title(network+'.'+sta+'.'+loc+' '+compstr+'\n'+tstart+' to '+tend)
# plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
figout = './figs/'
if not os.path.exists(figout):
    os.makedirs(figout)
plt.savefig(figout+network+'.'+sta+'.'+loc+'.'+compstr+'.pdf', dpi=300)
plt.show()



# %%
# Save median and quantiles to csv file
df_out = pd.DataFrame(
    {'freq': dfs.index.values,
     'psd_med': dfs.median(axis=1) ,
     'psd_2.5': dfs.quantile(q=0.025,axis=1),
     'psd_97.5': dfs.quantile(q=0.975,axis=1),
    })
df_out = df_out.reset_index()
df_out = df_out.drop('period', axis=1)

folder_to_create = './PSD_med/'
if not os.path.exists(folder_to_create):
    os.makedirs(folder_to_create)

df_out.to_csv(folder_to_create+network+'.'+sta+'.'+loc+'.'+compstr+'.csv', index=False) 

# # Save NLNM and NHNM to csv file
# df_nnm = pd.DataFrame(
#     {'freq': freqs,
#      'nlnm': nlnm_m_hz,
#      'nhnm': nhnm_m_hz,
#     })
# df_nnm = df_nnm.reset_index(drop=True)
# df_nnm.to_csv(folder_to_create+'NLNM_NHNM.csv', index=False)



# %%
