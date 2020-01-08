.. _wave:

Wave module
========================================

The :class:`~mhkit.wave` module contains a set of functions to
calculate quantities of interest for wave energy converters (WEC). 
The wave module contains the following submodules:

* :class:`~mhkit.wave.io`: Loads data from standard formats
* :class:`~mhkit.wave.resource`: Computes resource metrics such as spectra and significant wave height
* :class:`~mhkit.wave.performance`: Computes performance metrics such as capture length matrix and mean annual energy production
* :class:`~mhkit.wave.graphics`: Generates graphics

Data formats and IO
--------------------

As with other modules in the MHKiT Matlab package, the MHKit Python modules are used 
within the wave module.  The wave module uses wave elevation and spectra data.

* **Wave elevation data** is time-series data stored as an array which is converted to a pandas DataFrame with a time index.  The column names describe the type of data in each column.

* **Spectra data** can be stored in a structure with frequency, type, and spectrum field names or converted to a pandas DataFrame where the index is a `descriptor` and columns are
  frequency.  The descriptor(type)/frequency can be time (in seconds or a DateTime index), or some other information that
  defines each row of the spectra data.  For example, if spectra data is computed using the JONSWAP method,
  the index or type is a string that includes the significant wave height and peak period used to compute the
  spectra in each row.


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
   

