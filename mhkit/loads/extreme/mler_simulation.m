function sim = mler_simulation()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Define the simulation parameters that are used in various MLER
%     functionalities. If no input is given, then default values are returned.  
%     
%     Parameters
%     ----------
%         sim : struct (optional)
%             Simulation parameters
%             Keys:
%             -----
%             'startTime': starting time [s]
%             'endTime': ending time [s]
%             'dT': time-step size [s]
%             'T0': time of maximum event [s]
%             'startx': start of simulation space [m]
%             'endX': end of simulation space [m]
%             'dX': horizontal spacing [m]
%             'X': position of maximum event [m]
%         
%      Returns
%      -------
%         sim : struct 
%             Simulation parameters including spatial and time calculated
%             arrays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    sim.startTime = -150.0;  % .s Starting time
    sim.endTime = 150.0;  % .s Ending time
    sim.dT = 1.0;  % .s Time-step size
    sim.T0 = 0.0;  % .s Time of maximum event

    sim.startX = -300.0;  % .m Start of simulation space
    sim.endX = 300.0;  % .m End of simulation space
    sim.dX = 1.0;  % .m Horiontal spacing
    sim.X0 = 0.0;  % .m Position of maximum event
end


sim.maxIT = ceil((sim.endTime - sim.startTime)/sim.dT + 1);
sim.T = linspace(sim.startTime, sim.endTime, sim.maxIT);

sim.maxIX = ceil((sim.endX - sim.startX)/sim.dX + 1);
sim.X = linspace(sim.startX, sim.endX, sim.maxIX);


end


