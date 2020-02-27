Return to MHKiT `documentation <https://mhkit-code-hub.github.io/MHKiT/index.html>`_

Wave module
===================
The wave module contains a set of functions to
calculate quantities of interest for wave energy converters (WEC).

File IO
----------
The :class:`~mhkit.wave.io` submodule contains the following function to load National Data Buoy Center (NDBC) 
data file into a pandas DataFrame, including real time and historical data.

===========================================  =========================
Functions                                    Description
===========================================  =========================
read_NDBC_file                               Reads a NDBC wave buoy data file (from https://www.ndbc.noaa.gov) into a structure. 
===========================================  ========================= 


Wave resource
--------------------------------------

The :class:`~mhkit.wave.resource` submodule contains methods to compute wave energy spectra and various metrics from the spectra.

The following options exist to compute wave energy spectra:

===========================================  =========================
Functions                                    Description
===========================================  =========================
create_spectra                               Calculates a spectra of user defined type.
elevation_spectrum                           Calculates wave spectra from wave probe timeseries.
===========================================  ========================= 
   

The following metrics can be computed from the spectra:

===========================================  =========================
Functions                                    Description
===========================================  =========================
average_crest_period                         Calculate the average creat period from spectra. 
average_wave_period                          Calculates the average wave period from spectra
average_zero_crossing_period                 Calculates wave average zero crossing period from spectra
energy_flux                                  Calculates the omnidirectional wave energy flux of the spectra
energy_period                                Calculates the energy period
frequency_moment                             Calculates the Nth frequency moment of the spectrum
peak_period                                  Calculates wave energy period from spectra
significant_wave_height                      Calculates wave height from spectra
spectral_bandwidth                           Calculates bandwidth from spectra
spectral_width                               Calculates wave spectral width from spectra
surface_elevation                            Calculates wave elevation time series from spectrum using a random phase
wave_celerity                                Calculates wave celerity (group velocity)
wave_number                                  Calculates wave number
===========================================  ========================= 
                              


Wave performance
---------------------

The :class:`~mhkit.wave.performance` submodule contains functions to compute capture length, statistics, 
performance matrices, and mean annual energy production.

===========================================  =========================
Functions                                    Description
===========================================  =========================
dc_power                                     Calculates the real power from DC voltage and current. 
ac_power_three_phase                         Calculates the real power from three phase ac voltage and current. 
capture_length                               Calculates the capture length (often called capture width).
capture_length_matrix                        Generates a capture length matrix for a given statistic
mean_annual_energy_production_matrix         Calculates mean annual energy production (MAEP) from matrix data along with data frequency in each bin
mean_annual_energy_production_timeseeries    Calculates mean annual energy production (MAEP) from timeseries
power_matrix                                 Generates a power matrix from a capture length matrix and wave energy flux matrix
wave_energy_flux_matrix                      Generates a wave eneergy flux matrix for a given statistic
===========================================  ========================= 


Graphics
-----------

The :class:`~mhkit.wave.graphics` submodule contains functions to plot wave data and related metrics.  

===========================================  =========================
Functions                                    Description
===========================================  =========================
plot_elevation_timeseries                    Plots wave elevation timeseries 
plot_matrix                                  Plots the matrix with Hm0 and Te on the y and x axis 
plot_spectrum                                Plots wave amplitude spectrum
===========================================  ========================= 
   


.. Note::
    The names of the functions below are of the convention path.path.path.function. Only the function name is used when calling the function in MATLAB. For example, to call on "mhkit.wave.io.read_NDBC_file" simply 
    use read_NDBC_file(). 

.. automodule:: mhkit.wave
    :members:
    :no-undoc-members:
    :show-inheritance:
    
IO
-----------
.. automodule:: mhkit.wave.io
    :members:
    :undoc-members:
    :show-inheritance:
    
Resource
---------
.. automodule:: mhkit.wave.resource
    :members:
    :undoc-members:
    :show-inheritance:

Performance
------------
.. automodule:: mhkit.wave.performance
    :members:
    :undoc-members:
    :show-inheritance:
	
Graphics
----------- 
.. automodule:: mhkit.wave.graphics
    :members:
    :undoc-members:
    :show-inheritance: