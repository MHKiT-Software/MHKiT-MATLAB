diary 'test_log_1626.txt'
fnms = {dir('example_ncfiles/').name};
for ifnm = 3:numel(fnms)
    fnm = fnms{ifnm};
    fprintf("Checking File: %s \n",fnm);
    if strcmp(fnm,'QA4ECV_L2_NO2_OMI_20180301T052400_o72477_fitB_v1.nc')
        finfo = ncinfo(strcat('example_ncfiles/',fnm));
        ginfo = finfo.Groups(1).Groups(1).Groups(2);
        vnms = {ginfo.Variables.Name};
        vnms = strcat('PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/',vnms);
        res = read_nc_file_var(strcat('example_ncfiles/',fnm),vnms);
        sz = size(ginfo.Variables);
        idx = randi([1,sz(2)],1);
        tmplst = split(vnms{idx},'/');
        vname = check_name(tmplst{end});
    else 
        res = read_nc_file_var(strcat('example_ncfiles/',fnm));
        ginfo = ncinfo(strcat('example_ncfiles/',fnm));
        vnms = {ginfo.Variables.Name};
        sz = size(ginfo.Variables);
        idx = randi([1,sz(2)],1);
        vname = check_name(vnms{idx});
    end
    fprintf("Check variable: %s, idx = %d\n",vname,idx);
    %check dimension:
    fprintf("- check dimensions: ");
    if ~isempty(ginfo.Variables(idx).Dimensions)
        in_names = {ginfo.Variables(idx).Dimensions.Name};
        in_vals = {ginfo.Variables(idx).Dimensions.Length};
        out_names = res.(vname).dims;
        if ((~isequal(in_names,out_names)) || ...
                ~isequal( ...
                size(ncread(strcat('example_ncfiles/',fnm),vnms{idx})),...
                size(res.(vname).data)))
            ME = MException('MATLAB:test_read_nc_file',['failed  ' ...
            'test']);
            disp(in_names);
            disp(out_names);
            disp(size(ncread(strcat('example_ncfiles/',fnm),vnms{idx})));
            disp(size(res.(vname).data));
            throw(ME);
        else
            fprintf("Passed.\n");
        end
    else
        disp(res.(vname).dims);
       
    end
    fprintf("- check attributes: \n");
    if ~isempty(ginfo.Variables(idx).Attributes)

        in_names = {ginfo.Variables(idx).Attributes.Name};
        in_vals = {ginfo.Variables(idx).Attributes.Value};
        for iattr = 1:numel(in_names)
            fprintf("attr %s: ",in_names{iattr});
            if strcmp(in_names{iattr},'_FillValue')
                out_val = res.(vname).FillValue;
            else
                out_val = res.(vname).attrs.(in_names{iattr});
            end
            if(~isequal(in_vals{iattr},out_val) && ...
                    ~isnan(in_vals{iattr}) && ...
                    ~isnan(out_val))
                ME = MException('MATLAB:test_read_nc_file',['failed  ' ...
            'test']);
                disp(in_vals{iattr});
                disp(out_val);
                throw(ME);
            else
                fprintf("Passed. \n");
            end
        end

    else
        disp(ginfo.Variables(idx).Attributes);
    end
    disp(res.(vname).attrs);
end
diary off
