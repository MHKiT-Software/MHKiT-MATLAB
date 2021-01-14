function data=swan_read_block(swan_file)

%%%%%%%%%%%%%%%%%%%%
%     Reads in SWAN ASCII block format output
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

datatp = py.mhkit.wave.io.swan.read_block(swan_file)
datac=cell(datatp);
datapd=datac{1}
data_keys = cell(py.list(py.numpy.nditer(keys(datapd),pyargs("flags",{"refs_ok"}))));
keys2 = cell(py.list(py.numpy.nditer(data_keys{1},pyargs("flags",{"refs_ok"}))))
%disp(data_keys(2))
datamat=datac{2};
%matstr=struct(datamat)

data = 1;

% fid = fopen(swan_file);
% indata = fgetl(fid);
% while ischar(indata)
%      if startsWith(indata,"% Run")==1
%          meta = indata
%          break
%      end
%      indata = fgetl(fid);
%      
% end
% fclose(fid);
% 
% delimiterIn = ' ';
% mystructure = importdata(swan_file,delimiterIn)
% vars = string(strsplit(mystructure.textdata{5},' '))
% vars = vars(2:end-1)
% units = string(strsplit(mystructure.textdata{6},' '))
% units = units(2:end-1)
% data = mystructure;
%data = table2struct(array2table(mystructure.data,'VariableNames',vars),'ToScalar',true);

% for i = 1:max(size(units))
%     data.units.(string(vars(i))) = units(i);
% end
% meta = string(strsplit(strip(meta,'left','%'),' '));
% meta = meta(2:end);
% for j=1:max(size(meta))
%     if j ~= 3
%       m = string(strsplit(meta(j),':'));
%       data.meatdata.(m(1)) = m(2);
%     end
% end



