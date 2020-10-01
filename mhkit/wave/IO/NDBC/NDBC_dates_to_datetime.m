function ndbc_data=NDBC_dates_to_datetime(parameter,data,options)

%%%%%%%%%%%%%%%%%%%%
%     Takes a structure and converts the NDBC date fileds 
%     (e.g. "#YY  MM DD hh mm") to datetime. 
%     
%     
% Parameters
% ------------
%     parameter : string
%         'swden'	:	'Raw Spectral Wave Current Year Historical Data'
%         'stdmet':   'Standard Meteorological Current Year Historical Data'
%
%     data : Structure
%         Structure with fields (e.g. ['YY', 'MM', 'DD', 'hh', {'mm'}])  
%         
%
%     return_date_col : Logical (optional)
%         Default False. When true will return list of NDBC date columns
%     
% Returns
% ---------
%     ndbc_data: Structure 
%         Structure of structures broken down by years of data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    parameter
    data
    options.return_date_col = py.False;
end

py.importlib.import_module('mhkit');