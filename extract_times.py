# Dump a particular observing time from a uvh5 file. Write in UVFits format.

# https://pyuvdata.readthedocs.io/en/latest/uvdata_parameters.html
import sys, numpy as np
from pyuvdata import UVData

uvd = UVData()
uvd.read_uvh5(sys.argv[1], read_data=False)
times = np.unique(uvd.time_array)
freqs = uvd.freq_array[0]

print("Num times", len(times), "Num freqs", len(freqs))
print("Start time", times[0], "middle time", times[len(times)//2])
print("Start freq", freqs[0], "middle freq", freqs[len(uvd.freq_array[0])//2])

#uvd = UVData()
#uvd.read_uvh5(sys.argv[1], times=times[0], frequencies=freqs[0])
#uvd.write_uvfits("start.uvfits", force_phase=True, spoof_nonessential=True)

uvd = UVData()
uvd.read_uvh5(sys.argv[1], times=times[len(times)//2], frequencies=freqs[len(freqs)//2])
uvd.write_uvfits("middle.uvfits", force_phase=True, spoof_nonessential=True)


