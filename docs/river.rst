.. _river:

River module
========================================

The :class:`~mhkit/river` module contains a set of functions to 
calculate quantities of interest for river energy converters (REC). 
The river module contains the following submodules:

* :class:`~mhkit/river/io`: Loads data from standard formats
* :class:`~mhkit/river/resource`: Computes resource metrics such as exceedance probability, velocity, and power
* :class:`~mhkit/river/device`: Computes device metrics such as equivalent diameter and capture area
* :class:`~mhkit/river/graphics`: Generates graphics

Data format
--------------

The river module uses discharge data.
**Note, IEC/TS 62600-100 recommends that river resource calculations use 10 years of daily discharge data.**
      
* **Discharge time series data** is stored as a strucuture  
  Time can be specified in datetime or in seconds.  The column names describe the type of data in each column. 

File IO
--------------
The :class:`~mhkit/river/io` submodule contains the following functions to 
load USGS Discharge data into structures.  

.. autosummary::

   ~mhkit.river.io.read_usgs_file
   ~mhkit.river.io.request_usgs_data

River resource
---------------
The :class:`~mhkit.river.resource` submodule uses discharge data to compute 
exeedance probability, velocity, and power.  The module also contains functions
to compute the Froude number and to fit a polynomial to a series of points.
The polynomial is used to estimate the relationship between discharge and velocity 
or velocity and power at an individual turbine.

.. autosummary::

   ~mhkit.river.resource.Froude_number
   ~mhkit.river.resource.polynomial_fit
   ~mhkit.river.resource.exceedance_probability
   ~mhkit.river.resource.discharge_to_velocity
   ~mhkit.river.resource.velocity_to_power
 
.. TODO Add Annual energy produced (AEP)

River device 
---------------------
The :class:`~mhkit.river.device` submodule contains functions to compute equivalent diameter 
and capture area for circular, ducted, rectangular, adn multiple circular devices. 
A circular device is a vertical axis water turbine (VAWT). A 
rectangular device is a horizontal axis water turbine. A ducted device
is an enclosed VAWT. A multiple-circular devices is a device with
multiple VAWTs per device.

.. autosummary::

   ~mhkit.river.device.circular
   ~mhkit.river.device.ducted
   ~mhkit.river.device.rectangular
   ~mhkit.river.device.multiple_circular
   
Graphics
-------------
The :class:`~mhkit.river.graphics` submodule contains functions to plot river data and related metrics.  
The functions are designed to work in parallel with the :class:`~mhkit.river.resource` submodule.

.. autosummary::

   ~mhkit.river.graphics.plot_flow_duration_curve
   ~mhkit.river.graphics.plot_velocity_duration_curve
   ~mhkit.river.graphics.plot_power_duration_curve
   ~mhkit.river.graphics.plot_discharge_timeseries
   ~mhkit.river.graphics.plot_discharge_vs_velocity
   ~mhkit.river.graphics.plot_velocity_vs_power

