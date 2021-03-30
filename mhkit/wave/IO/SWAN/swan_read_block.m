function data=swan_read_block(swan_file)

%%%%%%%%%%%%%%%%%%%%
%     Reads in SWAN ASCII block format output and returns a data structure
%     containing modeled values and assocuated metadata
%     
% Parameters
% ------------
%     swan_file : string
%         SWAN file name to import    
%     
% Returns
% ---------
%     data: Structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datatp = py.mhkit.wave.io.swan.read_block(swan_file);
datac=cell(datatp);
datapd=datac{1};
datamat=datac{2};

for raw_key = py.list(keys(datapd))
   key = raw_key{1};
   value = datapd{key};
   new_string = regexprep(string(key),' ','_');
   sha=cell(value.values.shape);
   x=int64(sha{1,1});
   y=int64(sha{1,2});
   vals = double(py.array.array('d',py.numpy.nditer(value.values,pyargs("flags",{"refs_ok"}))));
   vals=reshape(vals,[x,y]);
   data.(new_string).values = vals.';
   meta = datamat{key};
   for mat_key = py.list(keys(meta))
       key_meta = mat_key{1};
       str_meta = string(key_meta);
       met_val = meta{key_meta};
       if isa(met_val,'py.str')
         met_val = string(met_val);
       end
       data.(new_string).(str_meta)= met_val;
   end
end
