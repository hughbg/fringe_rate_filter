# https://pyuvdata.readthedocs.io/en/latest/uvdata_parameters.html
import ephem, numpy as np
import sys
from pyuvdata import UVData

PI = 3.14159265359

def pol_str(pa):
  switcher = { 
    1: "pI",
    2: "pQ",
    3: "pU",
    4: "pV",
    -1: "RR",
    -2: "LL",
    -3: "RL",
    -4: "LR",
    -5: "XX",
    -6: "YY",
    -7: "XY",
    -8: "YX"
  } 

  return [ switcher[x] for x in pa ]
    
def julian_to_utc(j):
  DUBLIN_DATE = 2415020.00000
  d = ephem.Date(j-DUBLIN_DATE)
  return d

def print_uvdata(uvdata):
  
  print("Number of antennas with data present:", uvdata.Nants_data)
  print("Number of antennas in the array:", uvdata.Nants_telescope)
  print("Number of baselines:", uvdata.Nbls)
  print("Number of baseline-times:", uvdata.Nblts)
  print("Number of frequency channels:", uvdata.Nfreqs)
  print("Number of polarizations:", uvdata.Npols)
  print("Polarizations:", pol_str(uvdata.polarization_array))
  print("Number of spectral windows:", uvdata.Nspws)
  print("Number of times:", uvdata.Ntimes)
  print("ant_1_array:", uvdata.ant_1_array)
  print("ant_2_array:", uvdata.ant_2_array)
  print("Antenna names:", uvdata.antenna_names)
  print("Antenna numbers:", uvdata.antenna_numbers)
  print("Antenna positions:", uvdata.antenna_positions)
  print("Baseline array:", uvdata.baseline_array)
  print("Channel width:", "{:.2f}".format(uvdata.channel_width), "[Hz]")
  print("Frequencies:")
  for i in range(uvdata.Nspws):
    print("  SPW", str(i)+":", "{:.2f}".format(uvdata.freq_array[i, 0]), "-", "{:.2f}".format(uvdata.freq_array[i, -1]), "[Hz]")
  print("Integration time:", [ "{:.2f}".format(x) for x in np.unique(uvdata.integration_time) ], "[s]")
  print("Nsamples in integration:",uvdata.nsample_array)
  if len(np.unique(uvdata.lst_array)) == 1:
    print("LST:", "{:.2f}".format(12*uvdata.lst_array[0]/PI))
  else: print("LSTS:", "{:.2f}".format(12*uvdata.lst_array[0]/PI), "-", "{:.2f}".format(12*uvdata.lst_array[-1]/PI))
  if len(np.unique(uvdata.time_array)) == 1:
    print("Time:", julian_to_utc(uvdata.time_array[0]), "[UTC]")
  else: print("Times:", julian_to_utc(uvdata.time_array[0]), "-", julian_to_utc(uvdata.time_array[-1]), "[UTC]", 
          uvdata.time_array[0], "-", uvdata.time_array[-1], "[JD]")
  print("Telescope name:", uvdata.telescope_name)
  degree_sign = u"\N{DEGREE SIGN}"
  print("Telescope location:", "Lat", "{:.2f}".format(uvdata.telescope_location_lat_lon_alt_degrees[0])+degree_sign, 
	"Lon", "{:.2f}".format(uvdata.telescope_location_lat_lon_alt_degrees[1])+degree_sign,
	"Alt", "{:.2f}".format(uvdata.telescope_location_lat_lon_alt_degrees[2])+"m")
  print("Vis units:", uvdata.vis_units)
  print(uvdata.history)

uvd = UVData()
uvd.read_uvh5(sys.argv[1], read_data=False)
print_uvdata(uvd)
exit()

print("Start time", times[0], "middle time", times[len(times)//2])

uvd = UVData()
uvd.read_uvh5(sys.argv[1], times=times[0])
print_uvdata(uvd)
uvd.write_uvfits("start.uvfits", force_phase=True, spoof_nonessential=True)

uvd = UVData()
uvd.read_uvh5(sys.argv[1], times=times[len(times)//2])
print_uvdata(uvd)
uvd.write_uvfits("middle.uvfits", force_phase=True, spoof_nonessential=True)


