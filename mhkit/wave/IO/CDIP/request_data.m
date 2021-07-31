% url = 'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/430p1/430p1_d04.nc';
% url = 'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/430p1/430p1_d05.nc';
url = 'http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/430p1/430p1_d02.nc';
filename = 'ncdata.nc';
path = strcat(pwd , '/IO/CDIP/');

CDIP_request_data(strcat(path,filename), url); 

% Open netcdf4 file: 
ncid = netcdf.open(append(path, filename), 'NC_NOWRITE');

% sstTime: 
sstTime_varID = netcdf.inqVarID(ncid, 'sstTime');
sstTime_dimID = netcdf.inqDimID(ncid, 'sstTime');
[sstTime_dimname, sstTime_dimlength] = netcdf.inqDim(ncid, sstTime_dimID);

% Read the content of the variable:
sstTime_datenum = netcdf.getVar(ncid,sstTime_varID,0,sstTime_dimlength);

% Convert to inverse of POSIX time: 
sstTime = datetime(sstTime_datenum, 'ConvertFrom', 'posixtime');
% only for the data I have currently, which has different lengths of time
% and waveHs: 
% sstTime(31501) = [];

% waveHs
waveHs_varID = netcdf.inqVarID(ncid, 'waveHs');
waveHs_dimID = netcdf.inqDimID(ncid, 'waveTime');
[waveHs_dimname, waveHs_dimlength] = netcdf.inqDim(ncid, waveHs_dimID);

% Read the content of the variable:
waveHs = netcdf.getVar(ncid,waveHs_varID,0,waveHs_dimlength);

% waveTp
waveTp_varID = netcdf.inqVarID(ncid, 'waveTp');
waveTp_dimID = netcdf.inqDimID(ncid, 'waveTime');
[waveTp_dimname, waveTp_dimlength] = netcdf.inqDim(ncid, waveTp_dimID);

% Read the content of the variable:
waveTp = netcdf.getVar(ncid,waveTp_varID,0,waveTp_dimlength);


% waveDp
waveDp_varID = netcdf.inqVarID(ncid, 'waveDp');
waveDp_dimID = netcdf.inqDimID(ncid, 'waveTime');
[waveDp_dimname, waveDp_dimlength] = netcdf.inqDim(ncid, waveDp_dimID);

% Read the content of the variable:
waveDp = netcdf.getVar(ncid,waveDp_varID,0,waveDp_dimlength);

% plotting:
tiledlayout(3,1)
% data_title = ncreadatt(filename, '/', 'title');
title_varID = netcdf.inqVarID(ncid,'metaStationName');
plot_title = convertCharsToStrings(netcdf.getVar(ncid, title_varID));
sgtitle(plot_title)

nexttile
plot(sstTime, waveHs)
ylabel('Hs, m')
title('WaveHs vs. sstTime')

nexttile
plot(sstTime, waveTp)
ylabel('Tp, s')
title('WaveTp vs. sstTime')

nexttile
plot(sstTime, waveDp)
ylabel('Dp, deg')
title('WaveDp vs. sstTime')

figure; 

% Box plots: 
time = yyyymmdd(sstTime);
Hs = [time waveHs];
T = table();
T.year=floor(Hs(:,1)/10000);
T.mm=floor((Hs(:,1)-T.year*10000)/100);
T.dd=Hs(:,1)-T.year*10000-T.mm*100;
T.value=Hs(:,2);
boxplot(T.value, T.mm)
% boxplot(T.value, T.mm, 'Labels',{'Jan','Feb','Mar','Apr','May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
ylabel('Significant Wave Height, Hs (m)')
xlabel('Month #')
title('Significant Wave Height by month')

if exist(strcat(path, 'ncdata.nc'), 'file')==2
  delete(strcat(path, 'ncdata.nc'));
end

function [cdip_data, cdip_data_path]=CDIP_request_data(filename, url)

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

cdip_data_path = websave(filename, url);
cdip_data = filename;

end
