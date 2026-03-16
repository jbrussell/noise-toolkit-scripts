# Plot power-spectral densities

# %load_ext autoreload
# %autoreload
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
import subprocess

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

cmd = [
    "python3", path2bin + "ntk_computePSD.py",
    f"net={network}",
    f"sta={sta}",
    f"loc={loc}",
    f"chan={compstr}",
    f"start={tstart}",
    f"end={tend}",
    f"xtype={xtype}",
    f"sw_width={sw_width}",
    f"sw_shift={sw_shift}",
    "plot=0"
]

subprocess.run(cmd, check=True)

# %% codecell
# Gather PSDs

path2psdDb = './data/psdDb/'+network+'.'+sta+'.'+loc+'/'+compstr+'/'

pathlist = sorted(Path(path2psdDb).glob('*/*/*.'+xtype+'.txt'))

# Build daily datetime vector
datevec = np.arange(pd.to_datetime(tstart),pd.to_datetime(tend),timedelta(days=1))

list_of_dfs = []
freq_ref = None
period_ref = None
for day in datevec:
    doy = pd.to_datetime(day).dayofyear
    year = pd.to_datetime(day).year
    pathlist = sorted(Path(path2psdDb+str(year)+'/').glob('*'+str(doy)+'/*.'+xtype+'.txt'))

    for path in pathlist:
        df = pd.read_table(path, delim_whitespace=True,skiprows=[0],names=[xtype,'psd'])
        # sta = data.sta[0]
        df = df.set_index(xtype)

        periods = df.index.values
        freqs = 1/periods

        # Interpolate to common period axis
        if freq_ref is None:
            freq_ref = freqs
            period_ref = periods
        else:
            interp_vals = np.interp(np.squeeze(period_ref), np.squeeze(df.index.to_numpy()), np.squeeze(df.to_numpy()))
            s = pd.Series(interp_vals, index=period_ref)
            s.index.name = 'period'
            df = s

        list_of_dfs.append(df)

# Concatinate psd dataframe
dfs = pd.concat(list_of_dfs, ignore_index=True,axis=1)

# Drop nan values
dfs = dfs.dropna(how='all')

# %%
# Plot PSDs

plotunit = 2 # 0: displacement, 1: velocity, 2: acceleration

plt.figure(figsize=(10,6))
plt.semilogx(dfs, color='gray', alpha=0.5)
# Plot median
plt.semilogx(dfs.median(axis=1), color='red', alpha=1, linewidth=2, label='Median')
# Plot 97.5% quantile
plt.semilogx(dfs.quantile(q=0.975,axis=1), linestyle='--', color='red', alpha=1, linewidth=2, label='97.5')
# Plot 2.5% quantile
plt.semilogx(dfs.quantile(q=0.025,axis=1), linestyle='--', color='red', alpha=1, linewidth=2, label='2.5')
# Plot noise models
freqs = 1/dfs.index.values
plt.semilogx(1/freqs, 10*np.log10(accNLNM(freqs) / (2*np.pi*freqs)**(2-plotunit)),
              c='blue', label='NLNM')
plt.semilogx(1/freqs, 10*np.log10(accNHNM(freqs) / (2*np.pi*freqs)**(2-plotunit)),
            c='blue', label='NHNM')
plt.xlim([0.05,200])
plt.ylim([-200,0])
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 (m/s^2)^2/Hz)')
plt.title(network+'.'+sta+'.'+loc+' '+compstr)
# plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()


# %%
