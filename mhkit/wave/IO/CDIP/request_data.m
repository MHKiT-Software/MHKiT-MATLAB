url = 'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/433p1/433p1_historic.nc';
filename = 'ncdata.nc';


function cdip_data=CDIP_request_data(url,filename)

%%%%%%%%%%%%%%%%%%%%
%     Requests data by url
%     
%     
% Parameters
% ------------
%     parameter : string
%         'swden'	:	'Raw Spectral Wave Current Year Historical Data'
%         'stdmet':   'Standard Meteorological Current Year Historical Data'
%
%     filename  : string
%         
%     
% Returns
% ---------
%     cdip_data: netCDF4 
%         netCDF4 file
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdip_data = websave(filename, url);

end