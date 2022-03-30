import sys
import numpy as np
from pyuvdata import UVData

def add_noise_from_autos(uvd_in, input_noise=None, nsamp=1, seed=None, inplace=False):
    """
    Add noise to a simulation, using the autos to estimate the noise rms.
    
    Parameters
    ----------
    uvd_in : UVData
        Input UVData object.
    input_noise : string, optional
        path to noise data file, Default: None (use same visibility file)
    
    nsamp : float, optional
        Rescale the generated noise according to some effective number of 
        samples (e.g. to emulate the effect of LST binning). Default: 1.
    
    seed : int, optional
        Random seed to set before generating random noise. Default: None.
    
    inplace : bool, optional
        Whether to add the noise directly to the input UVData object, or 
        operate on a copy. Default: False (operate on a copy).
    
    Returns
    -------
    uvd : UVData
        Output UVData object, now with noise included.
    """
    # Set random seed
    np.random.seed(seed)
    
    # Make a copy of the UVData object if needed
    if inplace:
        uvd = uvd_in
    else:
        uvd = copy.deepcopy(uvd_in)
    
    # Get channel width and integration time
    dnu = uvd.channel_width # in Hz
    dt = uvd.integration_time[0] # in sec

    #read noise .uvh5 file if exists
    if input_noise is not None:
        uvd_n = UVData()
        uvd_n.read_uvh5(input_noise)

    # Get all autos
    v_auto = {}
    for ant in uvd.antenna_numbers:
        auto_idxs = uvd.antpair2ind(ant, ant)

        #fill autos from noise file if exists
        if (input_noise != None) and (ant in uvd_n.antenna_numbers):
            v_auto[ant] = uvd_n.data_array[auto_idxs]

        #else fill from same file
        else:
            v_auto[ant] = uvd.data_array[auto_idxs]
    
    # Loop over baselines and add noise to each
    for bl in np.unique(uvd.baseline_array):
        
        # Get antenna numbers and indices in data array
        ant1, ant2 = uvd.baseline_to_antnums(bl)
        bl_idxs = uvd.antpair2ind(ant1, ant2)
        
        # Construct rms for each baseline pair, based on autos
        noise_shape = uvd.data_array[bl_idxs].shape
        std_ij = np.sqrt(v_auto[ant1] * v_auto[ant2].conj() / (nsamp * dt * dnu))
        n = 1.0 * np.random.normal(size=std_ij.size).reshape(std_ij.shape) \
          + 1.j * np.random.normal(size=std_ij.size).reshape(std_ij.shape)
        n *= std_ij / np.sqrt(2.) # to handle real and imaginary contributions
        
        # Add noise realisation
        uvd.data_array[bl_idxs,:,:,:] += n
    
    # Rescale nsamples_array by the assumed number of samples
    uvd.nsample_array *= nsamp
    
    return uvd


uvd = UVData()
uvd.read_uvh5(sys.argv[1])
#print(uvd.data_array.shape)
#print(np.mean(np.abs(uvd.data_array.real)), np.mean(np.abs(uvd.data_array.imag)))

if np.any(np.isnan(uvd.data_array)):
    print("Warning: there are NaNs in the data")


"""
noise_scale = 500

noise = np.random.normal(size=uvd.data_array.shape, scale=noise_scale)+np.random.normal(size=uvd.data_array.shape, scale=noise_scale)*1j
uvd.data_array += noise
uvd.history += "Random noise added, scale "+str(noise_scale)+"\n"
"""

outf = sys.argv[1][:-5]+"_noise.uvh5"
add_noise_from_autos(uvd, inplace=True).write_uvh5(outf, clobber=True, data_compression="lzf")

