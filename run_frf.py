"""Script for running a fringe-rate filter on one of many files.

This script was specially-written for the H1C IDR3.2 validation test, in order
to remove artificial fringe-rate structure from simulations of diffuse emission.
Some of the default parameter settings were chosen based on some cursory
investigations into the fringe-rate structure of simulations run in late January
using v0.4.2 of the vis_cpu simulation package. The high fringe-rate structure was
traced back to an issue with sources at the horizon, but the exact issue is still
unresolved at the time of writing this script.

Some parameters for the simulation are as follows:
    17280 times (~5 second integration over a full day)
    1024 frequencies on the full H1C band (100-200 MHz)
    ~20 antennas (minimal layout to get all baselines), compressed by redundancy
    2008 GSM HEALPix map with NSIDE=128 -> ~0.5 degree resolution
    This gives a pixel crossing time of roughly 2 minutes, with a Nyquist frequency
    of about 4.5 mHz.
    The simulations used the pixel beam setting in vis_cpu, which maps the beam
    onto the image plane and interpolates on that grid, rather than directly
    interpolating the beam on the sky.
"""
import argparse
import numpy as np
from astropy import units
from astropy import constants
import uvtools
from pyuvdata import UVData
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
    "data_dir", type=str, help="Directory containing files to fringe-rate filter."
)
parser.add_argument(
    "file_index", type=int, help="Which of the sorted files to perform the FRF on."
)
parser.add_argument("outdir", type=str, help="Where to write the output data.")
parser.add_argument(
    "--max_freq", type=float, default=199.99e6, help="Maximum observed frequency in Hz."
)
parser.add_argument(
    "--frf_buffer", type=float, default=0, help="Extra buffer to FRF in mHz."
)
parser.add_argument(
    "--fr_min", type=float, default=1.0, help="Minimum half-width of FRF in mHz."
)
parser.add_argument(
    "--fr_max", type=float, default=5.85, help="Maximum half-width of FRF in mHz."
)
parser.add_argument(
    "--fft_taper", type=str, default="", help="Taper to use for FFTs."
)
parser.add_argument(
    "--pad_size", type=int, default=0, help="Amount of padding on either end of time series."
)
parser.add_argument(
    "--zeropad",
    default=False,
    action="store_true",
    help="Whether to zero-pad the time series.",
)
parser.add_argument(
    "--datapad",
    default=False,
    action="store_true",
    help="Whether to use the first/last `pad_size` samples to pad time series.",
)
parser.add_argument(
    "--clobber", default=False, action="store_true", help="Overwrite existing files"
)

if __name__=="__main__":
    args = parser.parse_args()
    if args.pad_size and not (args.zeropad or args.datapad):
        raise RuntimeError("Pad size provided, but neither padding option chosen.")
    
    # Filing information setup.
    data_files = sorted(list(Path(args.data_dir).iterdir()))
    infile = data_files[args.file_index]
    outfile = Path(args.outdir) / infile.name

    # Some parameter setup for the FRF.
    omega_e = np.array([0,1,0]) * units.cycle.to("rad") / units.sday.to("s")
    fr_prefac = args.max_freq * units.Hz.to("mHz") * omega_e / constants.c.si.value
    
    # Read in the data and extract relevant metadata.
    uvdata = UVData()
    uvdata.read(infile)
    times = np.unique(uvdata.time_array) * units.day.to("s")
    times -= times[0]
    antpos = dict(zip(*uvdata.get_ENU_antpos()[::-1]))

    # Modify the time axis if we're doing any padding.
    if args.pad_size:
        dt = times[1] - times[0]
        new_times = np.arange(1, args.pad_size + 1) * dt
        times = np.concatenate([-new_times[::-1], times, new_times + times[-1]])
        original_slice = slice(args.pad_size, -args.pad_size)
    else:
        original_slice = slice(None, None)

    # Figure out some FFT-related things.
    fringe_rates = np.fft.fftshift(
        uvtools.utils.fourier_freqs(times) * units.Hz.to("mHz")
    )
    if args.fft_taper:
        edgecut = args.pad_size if args.zeropad else 0
        taper = uvtools.dspec.gen_window(
            args.fft_taper, times.size, edgecut_low=edgecut, edgecut_hi=edgecut
        )
    else:
        taper = np.ones(times.size)

    # Perform the FRF for each baseline.
    for (ant1, ant2) in uvdata.get_antpairs():
        baseline = antpos[ant2] - antpos[ant1]
        max_frate = np.linalg.norm(np.cross(baseline, fr_prefac))
        filter_half_width = np.clip(
            max_frate + args.frf_buffer, args.fr_min, args.fr_max
        )
        frf = np.zeros(times.size, dtype=complex)
        frf[np.abs(fringe_rates) < filter_half_width] = 1
        for pol in uvdata.get_pols():
            antpairpol = (ant1, ant2, pol)
            data = uvdata.get_data(antpairpol).squeeze()
            inds, _, pol_inds = uvdata._key2inds(antpairpol)
            pol_ind = pol_inds[0][0]
            if args.pad_size:
                if args.zeropad:
                    pad = np.zeros(args.pad_size, dtype=complex)
                    pad_left = pad
                    pad_right = pad
                else:
                    pad_left = data[-args.pad_size:]
                    pad_right = data[:args.pad_size]
                data = np.concatenate([pad_left, data, pad_right])
            filtered_data = np.fft.ifft(
                np.fft.fft(data * taper) * frf
            ) / taper
            uvdata.data_array[inds,0,0,pol_ind] = filtered_data[original_slice]

    uvdata.write_uvh5(outfile, clobber=args.clobber)
