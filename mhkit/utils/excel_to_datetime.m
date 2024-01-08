function new_time=excel_to_datetime(excel_num)

%%%%%%%%%%%%%%%%%%%%
%     Convert excel datenum format to datetime
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
% py.importlib.import_module('numpy');


new_time = datetime(excel_num,'ConvertFrom','excel');

