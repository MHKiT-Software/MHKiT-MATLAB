function Cp=power_coefficient(power, inflow_speed, capture_area, rho)

%%%%%%%%%%%%%%%%%%%%
%     Function that calculates the power coefficient of MEC device
%     
%     
% Parameters
% ------------
%     power : vector
%         Power output signal of device after losses [W]
%         
%     inflow_speed : vector
%         Velocity of inflow condition [m/s]
%
%     capture_area : double or int
%         Projected area of rotor normal to inflow [m^2]
%
%     rho : double or int
%         Density of environment [kg/m^3]
%
% Returns
% ---------
%     Cp: vector 
%         Power coefficient of device [-]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py_tsr = py.mhkit.river.performance.tip_speed_ratio(py.numpy.array(rotor_speed),double(rotor_diameter),py.numpy.array(inflow_speed));

TSR = double(py.array.array('d',py.numpy.nditer(py_tsr,pyargs("flags",{"refs_ok"}))));