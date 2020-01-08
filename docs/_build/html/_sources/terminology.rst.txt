.. _terminology:

Terminology
===========

The following terminology is used in MHKiT-Python

======================= ===================================================
Term       	 	Definition [unit]
======================= ===================================================
:math:`A_P`		Projected capture area [m^2]
BS                 	Bretschneider spectrum
:math:`D_E` 		Equivalent diameter [m]
:math:`E`		Energy [J] 
:math:`\eta` 		Incident wave [m]
:math:`f` 		Frequency [Hz]
:math:`F` 		Exceedance probability [%]
Fr			Froude Number 
:math:`g` 		Gravity [m/s/s]
:math:`h` 		Water depth from bottom to water surface (e.g. SWL) [m]
:math:`H` 		Wave height [m]
:math:`H_{s}`		Significant wave height, mean wave height of the tallest third of waves [m]
:math:`H_{m0}`		Spectrally derived significant wave height [m]
:math:`J` 		Wave energy flux [W/m]
JS                 	JONSWAP spectrum
:math:`L` 		Capture length [m]
:math:`k` 		Wave number, :math:`k = \frac{2\pi}{\lambda}` [rad/m]
:math:`m` 		Mass [kg]
:math:`m_k`		Spectral moment of k, for k = 0,1,2,...
:math:`\omega` 		Wave frequency, :math:`\omega = \frac{2\pi}{T}` [rad/s]
:math:`P` 		Power [W]
PM                 	Pierson-Moskowitz specturm
:math:`Q`		Discharge [m^3/s] 
:math:`\rho` 		Density [kg/m3]
:math:`S` 		Spectral density [m^2/Hz]
SWL	 		Still Water Line
:math:`T_{e}` 		Energy period [s]
:math:`T_{m}` 		Mean wave period [s] 
:math:`T_{p}` 		Peak period [s]
:math:`T_{z}` 		Zero-crossing period [s]
:math:`v`		Velocity [m/s] 
:math:`V`		Velocity calculated for river and tidal modules [m/s] 
======================= ===================================================

Units
---------
The methods in MHKiT use the MKS (meters-kilograms-seconds) system, and  
assume data is stored in SI units, for example:

* Acceleration = :math:`m/s^2`
* Distance = :math:`m`
* Energy = :math:`J`
* Frequency = :math:`Hz` 
* Mass = :math:`kg`
* Power = :math:`W`
* Pressure = :math:`Pa`
* Time = :math:`s`
* Velocity = :math:`m/s`
* Voltage = :math:`V`
* Volume = :math:`m^3`

.. Note:: 
	How do we want to handle angles? Radians? Degrees?