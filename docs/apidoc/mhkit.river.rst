Return to MHKiT `documentation <https://mhkit-code-hub.github.io/MHKiT/index.html>`_

River module
===================
The river module contains a set of functions to 
calculate quantities of interest for river energy converters (REC).

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

.. Note::
    The names of the functions below are of the convention path.path.path.function. Only the function name is used when calling the function in MATLAB. For example, to call on "mhkit.river.io.request_usgs_data" simply 
    use request_usgs_data(). 

.. automodule:: mhkit.river
    :members:
    :no-undoc-members:
    
IO
-----------
.. automodule:: mhkit.river.io
    :members:
    :undoc-members:
    
Resource
---------
.. automodule:: mhkit.river.resource
    :members:
    :undoc-members:

Device
---------
.. automodule:: mhkit.river.device
    :members:
    :undoc-members:
	
Graphics
----------- 
.. automodule:: mhkit.river.graphics
    :members:
    :undoc-members: