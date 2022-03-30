import sys
import numpy as np
from pyuvdata import UVData


uvd = UVData()
uvd.read_uvh5(sys.argv[1])
#print(uvd.data_array.shape)
#print(np.mean(np.abs(uvd.data_array.real)), np.mean(np.abs(uvd.data_array.imag)))

if np.any(np.isnan(uvd.data_array)):
    print("Warning: there are NaNs in the data")

ant5 = uvd.get_data(53, 53, "XX")
if np.sum(ant5[:2, :2]) > 0: print("autos present")
else: print("No autos")
