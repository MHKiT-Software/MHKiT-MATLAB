function data=swan_read_table(swan_file)

%%%%%%%%%%%%%%%%%%%%
%     Reads in SWAN table format output
%     
% Parameters
% ------------
%     swan_file : string
%         SWAN file name to import    
%     
% Returns
% ---------
%     data: Structure 
%
%
%         data.Data: named according to header row 
%
%         data.units: units for each data field 
%
%         data.metadata: metadata for the SWAN run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(swan_file);
indata = fgetl(fid);
while ischar(indata)
     if startsWith(indata,"% Run")==1
         meta = indata;
         break
     end
     indata = fgetl(fid);
     
end
fclose(fid);

delimiterIn = ' ';
mystructure = importdata(swan_file,delimiterIn);
vars = string(strsplit(mystructure.textdata{5},' '));
vars = vars(2:end-1);
units = string(strsplit(mystructure.textdata{6},' '));
units = units(2:end-1);

data = table2struct(array2table(mystructure.data,'VariableNames',vars),'ToScalar',true);

for i = 1:max(size(units))
    data.units.(string(vars(i))) = units(i);
end

data.metadata = strip(meta,'left','%');


