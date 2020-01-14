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

.. autosummary::

   ~mhkit.wave.io.read_NDBC_file

Wave resource
--------------------------------------

The :class:`~mhkit.wave.resource` submodule contains methods to compute wave energy spectra and various metrics from the spectra.

The following options exist to compute wave energy spectra:

.. autosummary::

   ~elevation_spectrum
   ~create_spectra
   

The following metrics can be computed from the spectra:

.. autosummary::

   ~frequency_moment
   ~significant_wave_height
   ~average_zero_crossing_period
   ~average_crest_period
   ~average_wave_period
   ~peak_period
   ~energy_period
   ~spectral_bandwidth
   ~spectral_width
   ~energy_flux
   ~wave_celerity
   ~wave_number
                              


Wave performance
---------------------

The :class:`~mhkit.wave.performance` submodule contains functions to compute capture length, statistics, 
performance matrices, and mean annual energy production.

.. autosummary::

   ~capture_length
   ~capture_length_matrix
   ~wave_energy_flux_matrix
   ~power_matrix
   ~mean_annual_energy_production_timeseries
   ~mean_annual_energy_production_matrix


Graphics
-----------

The :class:`~mhkit.wave.graphics` submodule contains functions to plot wave data and related metrics.  

.. autosummary::

   ~plot_elevation_timeseries
   ~plot_spectrum
   ~plot_matrix
   

