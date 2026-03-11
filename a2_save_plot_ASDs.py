# Plot amplitude-spectral densities and save quantiles defined between tstart and tend

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
from datetime import datetime, timedelta

from noisemodels import accNLNM,accNHNM

# %% codecell


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

plotunit = 0 # 0: displacement (m/sqrt(Hz)), 1: velocity, 2: acceleration

# Build daily datetime vector
datevec = np.arange(pd.to_datetime(tstart),pd.to_datetime(tend),timedelta(days=1))

list_of_dfs = []
for day in datevec:
    doy = pd.to_datetime(day).dayofyear
    pathlist = sorted(Path(path2psdDb).glob('*/'+str(doy)+'/*.'+xtype+'.txt'))

    for path in pathlist:
        df = pd.read_table(path, delim_whitespace=True,skiprows=[0],names=[xtype,'psd'])
        # sta = data.sta[0]
        df = df.set_index(xtype)

        freqs = 1/df.index.values

        df = 10**(df['psd']/20)/(2*np.pi*freqs)**(2-plotunit) # convert dB rel m**2/s**4/Hz --> m/Hz**0.5 (not log)

        list_of_dfs.append(df)

# Concatinate psd dataframe
dfs = pd.concat(list_of_dfs, ignore_index=True,axis=1)

# Drop nan values
dfs = dfs.dropna()

# %%
# Convert everything to m/sqrt(Hz)

dfs.index = 1/dfs.index
freqs = dfs.index.values

nlnm_db = 10*np.log10(accNLNM(freqs))
nlnm_m_hz = 10**(nlnm_db/20)/(2*np.pi*freqs)**(2-plotunit) # convert dB rel m**2/s**4/Hz --> m/Hz**0.5 (not log)

nhnm_db = 10*np.log10(accNHNM(freqs))
nhnm_m_hz = 10**(nhnm_db/20)/(2*np.pi*freqs)**(2-plotunit) # convert dB rel m**2/s**4/Hz --> m/Hz**0.5 (not log)

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
plt.ylim([1e-12,1e-4])
plt.xlabel('Frequency (Hz)')
plt.ylabel(r'$m/\sqrt{Hz}$')
plt.title(network+'.'+sta+'.'+loc+' '+compstr+'\n'+tstart+' to '+tend)
# plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
figout = './figs/'
if not os.path.exists(figout):
    os.makedirs(figout)
plt.savefig(figout+network+'.'+sta+'.'+loc+'.'+compstr+'.'+datevec[0].item().strftime('%Y-%m-%d')+'.'+datevec[-1].item().strftime('%Y-%m-%d')+'_ASD.pdf', dpi=300)
plt.show()



# %%
# Save median and quantiles to csv file
df_out = pd.DataFrame(
    {'freq': dfs.index.values,
     'asd_1': dfs.quantile(q=0.01,axis=1),
     'asd_2.5': dfs.quantile(q=0.025,axis=1),
     'asd_5': dfs.quantile(q=0.05,axis=1),
     'asd_10': dfs.quantile(q=0.1,axis=1),
     'asd_25': dfs.quantile(q=0.25,axis=1),
     'asd_50': dfs.quantile(q=0.5,axis=1),
     'asd_75': dfs.quantile(q=0.75,axis=1),
     'asd_90': dfs.quantile(q=0.90,axis=1),
     'asd_95': dfs.quantile(q=0.95,axis=1),
     'asd_97.5': dfs.quantile(q=0.975,axis=1),
     'asd_99': dfs.quantile(q=0.99,axis=1),
    })
df_out = df_out.reset_index()
df_out = df_out.drop('period', axis=1)

folder_to_create = './quantiles_ASD/'
if not os.path.exists(folder_to_create):
    os.makedirs(folder_to_create)

df_out.to_csv(folder_to_create+network+'.'+sta+'.'+loc+'.'+compstr+'.'+datevec[0].item().strftime('%Y-%m-%d')+'.'+datevec[-1].item().strftime('%Y-%m-%d')+'.csv', index=False) 

# # Save NLNM and NHNM to csv file
# df_nnm = pd.DataFrame(
#     {'freq': freqs,
#      'nlnm': nlnm_m_hz,
#      'nhnm': nhnm_m_hz,
#     })
# df_nnm = df_nnm.reset_index(drop=True)
# df_nnm.to_csv(folder_to_create+'NLNM_NHNM.csv', index=False)



# %%
