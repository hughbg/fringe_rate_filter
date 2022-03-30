from mpi4py import MPI
import matplotlib.pyplot as plt
from hera_sim.visibilities import VisCPU
from hera_sim import io
from hera_sim.beams import PerturbedPolyBeam
from astropy.coordinates import Latitude, Longitude
from astropy.units import Quantity
import numpy as np
from pyuvsim import uvsim
from pyuvdata import UVBeam
from pyuvsim.telescope import BeamList
from pyradiosky import SkyModel
from pyuvsim import simsetup, AnalyticBeam
import time

from vis_cpu.conversions import equatorial_to_eci_coords
from astropy.coordinates import EarthLocation
from astropy.time import Time

plt.rcParams['figure.figsize'] = [12, 6]

def compare_val(puvsim_uvd, hera_sim_uvd, ant1=0, ant2=1, which_freq=0, which_time=0):
    """
    Compare the outputs of two simulations, using a couple of values.
    puvsim_uvd, hera_sim_uvd are UVData objects from the simulations. Extract the 
    autocorrelation and the cross-correlation between antennas
    and print the values. Values are converted to amplitude and phase.
    Print a table with values from pyuvsim in the first column and hera_sim in the second
    column.
    """
   
    sim1_auto = puvsim_uvd.get_data(ant1, ant1, "XX")[which_time][which_freq]
    sim1_cross = puvsim_uvd.get_data(ant1, ant2, "XX")[which_time][which_freq]
    sim2_auto = hera_sim_uvd.get_data(ant1, ant1, "XX")[which_time][which_freq]
    sim2_cross = hera_sim_uvd.get_data(ant1, ant2, "XX")[which_time][which_freq]
    
    print("\n---------------------------------- Accuracy --------------------------------------")
    print("Values from Channel", which_freq, "Time", which_time)
    print("\t\t\t\t pyuvsim \t\t\t hera_sim")
    print("Auto corr ant "+str(ant1), "\t\t", abs(sim1_auto), "\t\t", abs(sim2_auto)) 
    print("Cross corr ant ("+str(ant1)+","+str(ant2)+") Amp", "\t", abs(sim1_cross), "\t\t", abs(sim2_cross))
    print("Cross corr ant ("+str(ant1)+","+str(ant2)+") Phase", "\t", np.angle(sim1_cross), "\t\t", np.angle(sim2_cross))
    print("\n\n")

def compare_baseline(pyuvsim_uvd, hera_sim_uvd, ant1=0, ant2=0, which_time=0):
    """
    Plot amplitude and phase of baseline. 
    """
    print("Plots of baseline("+str(ant1)+","+str(ant2)+")", "Time", which_time)
    
    f, (ax1, ax2) = plt.subplots(1, 2)
    
    ax1.plot(np.abs(pyuvsim_uvd.get_data(ant1, ant2, "XX")[which_time]), label="pyuvsim")
    ax1.plot(np.abs(hera_sim_uvd.get_data(ant1, ant2, "XX")[which_time]), label="hera_sim")
    ax1.set_ylabel("Amplitude")
    ax1.set_xlabel("Channel")
    ax1.legend()
    ax1.set_title("Amplitude for baseline ("+str(ant1)+","+str(ant2)+")")
 
    ax2.plot(np.angle(pyuvsim_uvd.get_data(ant1, ant2, "XX")[which_time]), label="pyuvsim")
    ax2.plot(np.angle(hera_sim_uvd.get_data(ant1, ant2, "XX")[which_time]), label="hera_sim")
    ax2.set_ylabel("Phase")
    ax2.set_xlabel("Channel")
    ax2.legend()
    ax2.set_title("Phase for baseline ("+str(ant1)+","+str(ant2)+")")
    plt.tight_layout()
    return f
 
def compare(pyuvsim_uvd, hera_sim_uvd, pyuvsim_time, hera_sim_time, ant1=0, ant2=0, which_freq=0, which_time=0):
    print("pyuvsim time:", pyuvsim_time)
    print("hera_sim time:", hera_sim_time)
    compare_val(pyuvsim_uvd, hera_sim_uvd, ant1, ant2, which_freq, which_time)
    #compare_baseline(pyuvsim_uvd, hera_sim_uvd, ant1, ant2, which_time)

RA = 40
def telescope_config(which_beam, which_package, nant=2, nfreq=2, ntime=1, nsource=1):
    """
    Setup the configuration parameters for pyuvsim/hera_sim/healvis.
    Different packages require different objects for simulation.
    healvis not used here.
    """
    if which_package not in ["hera_sim", "pyuvsim", "healvis"]:
        raise ValueError("Unknown package: " + which_package)
    if which_beam not in ["AnalyticBeam", "PerturbedPolyBeam"]:
        raise ValueError("Unknown beam: " + which_beam)

    #np.random.seed(10)  # So we always get the same random values

    # Random antenna locations
    x = np.random.random(nant) * 400  # Up to 400 metres
    y = np.random.random(nant) * 400
    z = np.random.random(nant)
    ants = {}
    for i in range(nant):
        ants[i] = (x[i], y[i], z[i])

    # Observing parameters in a UVData object.
    uvdata = io.empty_uvdata(
          Nfreqs=nfreq,  # Need 2 freqs for healvis to work
          start_freq=100e6,
          channel_width=97.3e3,
          start_time=2458902.4,
          integration_time=40.0,
          Ntimes=ntime,
          array_layout=ants,
          polarization_array=np.array(["XX", "YY", "XY", "YX"]),
          Npols=4,
          x_orientation="east"
    )
    
    # Random sources, 1 fixed.
    RA = 135.7
    DEC = -30.72
    sources = [
        [ RA, DEC, 2, 0 ],    
        ]
    if nsource > 1:                  # Add random others
        ra = RA+(np.random.random(nsource-1)-0.5)*100
        dec = DEC+(np.random.random(nsource-1)-0.5)*100
        #flux = np.random.random(nsource-1)*4
        flux = np.full(nsource-1, 2)
        for i in range(nsource-1): sources.append([ ra[i], dec[i], flux[i], 0])

    sources = np.array(sources)


    # Source locations and frequencies.
    ra_dec = np.deg2rad(sources[:, :2])
    freqs = np.unique(uvdata.freq_array)


    if which_package == "hera_sim":
        # calculate source fluxes for hera_sim. pyuvsim does it a different way.
        flux = (freqs[:, np.newaxis] / freqs[0]) ** sources[:, 3].T * sources[:, 2].T
        beam_ids = list(ants.keys())

    # Beam model. PerturbedPolyBeam, which is not symmetrical.
    cfg_beam = dict(ref_freq=1.e8,
                    spectral_index=-0.6975,
                    mainlobe_width=0.3,
                    beam_coeffs=[0.29778665, -0.44821433, 0.27338272,
                                 -0.10030698, -0.01195859, 0.06063853,
                                 -0.04593295, 0.0107879, 0.01390283,
                                 -0.01881641, -0.00177106, 0.01265177,
                                 -0.00568299, -0.00333975, 0.00452368,
                                 0.00151808, -0.00593812, 0.00351559
                                 ])
    if which_beam == "AnalyticBeam":
        beam = [AnalyticBeam("uniform") for i in range(len(ants.keys())) ]

        beam_dict = {}
        for i in range(len(beam)): beam_dict[str(i)] = i
        
    if which_beam == "PerturbedPolyBeam":
        beam = [PerturbedPolyBeam(perturb_coeffs=np.array([-0.20437532, -0.4864951, -0.18577532, -0.38053642, 0.08897764, 0.06367166,
                                        0.29634711, 1.40277112]),
                              mainlobe_scale=1.0, xstretch=0.9, ystretch=0.8, **cfg_beam)
            for i in range(len(ants.keys()))]
        beam_dict = {}
        for i in range(len(beam)): beam_dict[str(i)] = i

    """
    uvb = UVBeam()
    uvb.read_beamfits("/lustre/aoc/projects/hera/H4C/beams/NF_HERA_Vivaldi_efield_beam_healpix.fits")
    uvb.freq_interp_kind = "linear"
    uvb.pixel_coordinate_system = "az_za"
    uvb.select(axis1_inds=np.arange(360), axis2_inds=np.arange(90))
    beam = [uvb for i in range(len(ants.keys())) ]
    beam_dict = {}
    for i in range(len(beam)): beam_dict[str(i)] = i
    """

    # That's enough for hera_sim, but extra objects are required for pyuvsim and healvis.

    if which_package == "pyuvsim":
        # Need a sky model.

        # Stokes for the first frequency only. Stokes for other frequencies
        # are calculated later.
        stokes = np.zeros((4, 1, ra_dec.shape[0]))
        stokes[0, 0] = sources[:, 2]
        reference_frequency = np.full(len(ra_dec), freqs[0])

        # Setup sky model.
        sky_model = SkyModel(name=[str(i) for i in range(len(ra_dec))],
                             ra=Longitude(ra_dec[:, 0], "rad"), dec=Latitude(ra_dec[:, 1], "rad"),
                             spectral_type="spectral_index",
                             spectral_index=sources[:, 3],
                             stokes=stokes,
                             reference_frequency=Quantity(reference_frequency, "Hz")
                             )

        # Calculate stokes at all the frequencies.
        sky_model.at_frequencies(Quantity(freqs, "Hz"), inplace=True)

    if which_package == "healvis":     # NOT USED
        # Need a GSM model and an Observatory.

        baselines = []
        for i in range(len(ants)):
            for j in range(i + 1, len(ants)):
                bl = hv.observatory.Baseline(ants[i], ants[j], i, j)
                baselines.append(bl)

        times = np.unique(uvdata.get_times("XX"))

        obs_latitude = -30.7215277777
        obs_longitude = 21.4283055554
        obs_height = 1073

        # create the observatory
        fov = 360  # Deg
        obs = hv.observatory.Observatory(obs_latitude, obs_longitude, obs_height, array=baselines, freqs=freqs)
        obs.set_pointings(times)
        obs.set_fov(fov)
        obs.set_beam(beam)

        gsm = hv.sky_model.construct_skymodel('gsm', freqs=freqs, Nside=64)

    # Return what is relevant for each package, pyuvsim or hera_sim
    if which_package == "hera_sim":
        return uvdata, beam, beam_dict, freqs, ra_dec, flux
    elif which_package == "pyuvsim":
        return uvdata, beam, beam_dict, sky_model
    elif which_package == "healvis":
        return obs, gsm
 
comm = MPI.COMM_WORLD

uvdata, beam, beam_dict, freqs, ra_dec, flux = telescope_config("AnalyticBeam", "hera_sim", nant=120, nfreq=32, ntime=1, nsource=100)

location = EarthLocation.from_geodetic(lat=-30.7215, lon=21.4283,  height=1073.)         # HERA telescope location
obstime = Time(uvdata.time_array[0], format='jd', scale='utc')

ra, dec = equatorial_to_eci_coords(ra_dec[:, 0], ra_dec[:, 1], obstime, location)
ra_dec[:, 0] = ra
ra_dec[:, 1] = dec
    
simulator = VisCPU(
        uvdata = uvdata,
        beams = beam,
        beam_ids = list(beam_dict.values()),
        sky_freqs = freqs,
        point_source_pos = ra_dec,
        point_source_flux = flux,
        bm_pix = 20,
        use_pixel_beams = False,
        mpi_comm=comm,
        precision = 2,
        polarized = True
    ) 

start = time.time()
simulator.simulate() 
hera_sim_time = time.time()-start

uvdata.write_uvfits("sim.uvfits", force_phase=True, spoof_nonessential=True)
exit()

uvdata, beam, beam_dict, sky_model = telescope_config("AnalyticBeam", "pyuvsim", nant=2, nfreq=1, ntime=1, nsource=1)
start = time.time()
pyuvsim_uvd = uvsim.run_uvdata_uvsim(uvdata, BeamList(beam), beam_dict=beam_dict,
         catalog=simsetup.SkyModelData(sky_model))
pyuvsim_time = time.time()-start
compare(pyuvsim_uvd, simulator.uvdata, pyuvsim_time, hera_sim_time, ant1=0, ant2=1, which_freq=0, which_time=0)

#pyuvsim_uvd.write_uvfits("x.uvfits", force_phase=True, spoof_nonessential=True)


#perturb_coeffs=   and polarization=
