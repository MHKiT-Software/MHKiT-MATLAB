function new_time=excel2datetime(excel_num)

%%%%%%%%%%%%%%%%%%%%
%     conversion of excel datenum format to datetime 
%     
% Parameters
% ------------
%     excel_num: vector  
%         vector of excel datenums to be converted 
%
% Returns
% ---------
%     time: DateTimeIndex 
%         vector of corresponding python datetime values
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

new_time = double(py.array.array('d',py.numpy.nditer(py.mhkit.utils.excel2datetime(py.numpy.array(excel_num)))))/1e9;
new_time = datetime(new_time,'ConvertFrom','posix');



