function data=swan_read_table(swan_file)

%%%%%%%%%%%%%%%%%%%%
%     Reads in SWAN ASCII table format output and returns a Matlab
%     structure with modeled data and assocaited metadata. 
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
meta = string(strsplit(strip(meta,'left','%'),' '));
meta = meta(2:end);
for j=1:max(size(meta))
    if j ~= 3
      m = string(strsplit(meta(j),':'));
      data.meatdata.(m(1)) = m(2);
    end
end



