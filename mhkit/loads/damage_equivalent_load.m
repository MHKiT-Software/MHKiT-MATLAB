function DEL=damage_equivalent_load(var,m,options)

%%%%%%%%%%%%%%%%%%%%
%     Calculates the damage equivalent load 
%     
% Parameters
% ------------
%     var: vector  
%         data of variable/channel being analyzed 
%
%     m : double or int 
%         fatigue slope factor of material  
%
%     bin_num : int (optional) 
%         number of bins for rainflow counting method (minimum=100)
%         to call: get_DELs(data,chan_info,"binNum",binNum)
%
%     data_length : double or int (optional)
%         Used to control DEL frequency. Length of data in sec. Default for 1Hz is 600 seconds for 10min data
%         to call: get_DELs(data,chan_info,"data_length",t)
%     
% Returns
% ---------
%     DEL: Structure 
%         Damage equivalent load of each specified variable
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    var 
    m
    options.bin_num {mustBeNumeric} = 100;
    options.data_length {mustBeNumeric} = 600;
end

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');


DEL = py.mhkit.loads.damage_equivalent_load(var,m,pyargs('bin_num',int32(options.bin_num),'data_length',int32(options.data_length)));


