var sourceData138 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB-2/mhkit/wave/IO/SWAN/swan_read_table.m","RawFileContents":["function data=swan_read_table(swan_file)","","%%%%%%%%%%%%%%%%%%%%","%     Reads in SWAN ASCII table format output and returns a Matlab","%     structure with modeled data and assocaited metadata.","%","% Parameters","% ------------","%     swan_file : string","%         SWAN file name to import","%","% Returns","% ---------","%     data: Structure","%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","","fid = fopen(swan_file);","indata = fgetl(fid);","while ischar(indata)","     if startsWith(indata,\"% Run\")==1","         meta = indata;","         break","     end","     indata = fgetl(fid);","","end","fclose(fid);","","delimiterIn = ' ';","mystructure = importdata(swan_file,delimiterIn);","vars = string(strsplit(mystructure.textdata{5},' '));","vars = vars(2:end-1);","units = string(strsplit(mystructure.textdata{6},' '));","units = units(2:end-1);","","data = table2struct(array2table(mystructure.data,'VariableNames',vars),'ToScalar',true);","","for i = 1:max(size(units))","    data.units.(string(vars(i))) = units(i);","end","meta = string(strsplit(strip(meta,'left','%'),' '));","meta = meta(2:end);","for j=1:max(size(meta))","    if j ~= 3","      m = string(strsplit(meta(j),':'));","      data.metadata.(m(1)) = m(2);","    end","end","",""],"CoverageDisplayDataPerLine":{"Function":{"LineNumber":1,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":40,"ContinuedLine":false},"Statement":[{"LineNumber":18,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":23,"ContinuedLine":false},{"LineNumber":19,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":20,"ContinuedLine":false},{"LineNumber":20,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":20,"ContinuedLine":false},{"LineNumber":21,"Hits":3,"StartColumnNumbers":5,"EndColumnNumbers":37,"ContinuedLine":false},{"LineNumber":22,"Hits":1,"StartColumnNumbers":9,"EndColumnNumbers":23,"ContinuedLine":false},{"LineNumber":23,"Hits":1,"StartColumnNumbers":9,"EndColumnNumbers":15,"ContinuedLine":false},{"LineNumber":25,"Hits":2,"StartColumnNumbers":5,"EndColumnNumbers":25,"ContinuedLine":false},{"LineNumber":28,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":12,"ContinuedLine":false},{"LineNumber":30,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":18,"ContinuedLine":false},{"LineNumber":31,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":48,"ContinuedLine":false},{"LineNumber":32,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":53,"ContinuedLine":false},{"LineNumber":33,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":21,"ContinuedLine":false},{"LineNumber":34,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":54,"ContinuedLine":false},{"LineNumber":35,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":23,"ContinuedLine":false},{"LineNumber":37,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":88,"ContinuedLine":false},{"LineNumber":39,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":26,"ContinuedLine":false},{"LineNumber":40,"Hits":6,"StartColumnNumbers":4,"EndColumnNumbers":44,"ContinuedLine":false},{"LineNumber":42,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":52,"ContinuedLine":false},{"LineNumber":43,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":19,"ContinuedLine":false},{"LineNumber":44,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":23,"ContinuedLine":false},{"LineNumber":45,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":13,"ContinuedLine":false},{"LineNumber":46,"Hits":3,"StartColumnNumbers":6,"EndColumnNumbers":40,"ContinuedLine":false},{"LineNumber":47,"Hits":3,"StartColumnNumbers":6,"EndColumnNumbers":34,"ContinuedLine":false}]}}