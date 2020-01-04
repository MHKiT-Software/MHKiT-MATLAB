function parameters = initParameters()
%
% parameters = initParameters();
%
%
% Initalizes a single instance of a parameters structure
%
% Input:
%   none;
%
% Output:
%   powerData                     MHKiT structure
%
% Dependencies
%   none
%
% Usage
%   parameters = initParameters();
%
% Version 1.01, 11/25/2018Rick Driscoll, NREL

parameters.structType   = 'Parameters'; % type of structure
parameters.g            = 9.81;         % gravitational constant
parameters.waterDensity = 1025;         % density of water
parameters.waterDepth   = [];           % depth of water at measurement location
         