.. _river:

River module
--------------

The :class:`~mhkit.river` module contains a set of functions to 
calculate quantities of interest for river energy converters (REC). 
The river module contains the following submodules:

* :class:`~mhkit.river.io`: Loads data from standard formats
* :class:`~mhkit.river.resource`: Computes resource metrics such as exceedance probability, velocity, and power
* :class:`~mhkit.river.device`: Computes device metrics such as equivalent diameter and capture area
* :class:`~mhkit.river.graphics`: Generates graphics

Data format
^^^^^^^^^^^^

The river module uses discharge data.
**Note, IEC/TS 62600-100 recommends that river resource calculations use 10 years of daily discharge data.**
      
* **Discharge time series data** is stored as a strucuture  
  Time can be specified in datetime or in seconds.  The column names describe the type of data in each column. 

IO
^^^^^^^^^^^^
The :class:`~mhkit.river.io` submodule contains the following functions to 
load USGS Discharge data into structures.  

===========================================  =========================
Functions                                    Description
===========================================  =========================
read_usgs_file                               Reads a USGS JSON data file (from https://waterdata.usgs.gov/nwis) into a structure 
request_usgs_file                            Loads USGS data directly from https://waterdata.usgs.gov/nwis using a GET request into a structure
===========================================  ========================= 

Resource
^^^^^^^^^^
The :class:`~mhkit.river.resource` submodule uses discharge data to compute 
exeedance probability, velocity, and power.  The module also contains functions
to compute the Froude number and to fit a polynomial to a series of points.
The polynomial is used to estimate the relationship between discharge and velocity 
or velocity and power at an individual turbine.


===========================================  =========================
Functions                                    Description
===========================================  =========================
Froude_number                                Calculate the Froude Number of the river, channel or duct flow, to check subcritical flow assumption (if Fr <1).
polynomial_fit                               Returns a polynomial fit for y given x of order n with an R-squared score of the fit
exceedance_probability                       Calculates the exceedance probability
discharge_to_velocity                        Calculates velocity given discharge data and the relationship between discharge and velocity at an individual turbine
velocity_to_power                            Calculates power given velocity data and the relationship between velocity and power from an individual turbine
energy_produced                              Returns the energy produced for a given time period provided exceedence probability and power.
===========================================  ========================= 

Device 
^^^^^^
The :class:`~mhkit.river.device` submodule contains functions to compute equivalent diameter 
and capture area for circular, ducted, rectangular, adn multiple circular devices. 
A circular device is a vertical axis water turbine (VAWT). A 
rectangular device is a horizontal axis water turbine. A ducted device
is an enclosed VAWT. A multiple-circular devices is a device with
multiple VAWTs per device.

===========================================  =========================
Functions                                    Description
===========================================  =========================
circular                                     Calculates the equivalent diameter and projected capture area of a circular turbine
ducted                                       Calculates the equivalent diameter and projected capture area of a ducted turbine
multiple_circular                            Calculates the equivalent diameter and projected capture area of a multiple circular turbine
rectangular                                  Calculates the equivalent diameter and projected capture area of a retangular turbine
===========================================  ========================= 

   
Graphics
^^^^^^^^^^
The :class:`~mhkit.river.graphics` submodule contains functions to plot river data and related metrics.  
The functions are designed to work in parallel with the :class:`~mhkit.river.resource` submodule.

===========================================  =========================
Functions                                    Description
===========================================  =========================
plot_discharge_timeseries                    Plots discharge vs time
plot_discharge_vs_velocity                   Plots discharge vs velocity
plot_flow_duration_curve                     Plots discharge vs exceedance probability as a Flow Duration Curve (FDC)
plot_power_duration_curve                    Plots power vs exceedance probability as a Flow Duration Curve (FDC)
plot_velocity_duration_curve                 Plots velocity vs exceedance probability as a Flow Duration Curve (FDC)
plot_velocity_vs_power                       Plots velocity vs power along with a polynomial fit 
===========================================  ========================= 

