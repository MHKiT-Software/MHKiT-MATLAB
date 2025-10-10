function [vec_data, vec_sysdata, awac_profile] = nortek_defs()
% Returns nortek instrument type definitions
%
% Returns
% -------
% vec_data : struct
%     Vector data type definitions
% vec_sysdata : struct
%     Vector system data type definitions
% awac_profile : struct
%     AWAC profile data type definitions

vec_data = struct();
vec_sysdata = struct();
awac_profile = struct();

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                        vec_data
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
vec_data.('AnaIn2LSB').('dims')  = {};
vec_data.('AnaIn2LSB').('dtype') = "uint8";
vec_data.('AnaIn2LSB').('group') = "sys";
vec_data.('AnaIn2LSB').('units') = "1";

vec_data.('Count').('dims')  = {};
vec_data.('Count').('dtype') = "uint8";
vec_data.('Count').('group') = "sys";
vec_data.('Count').('units') = "1";

vec_data.('PressureMSB').('dims')  = {};
vec_data.('PressureMSB').('dtype') = "uint8";
vec_data.('PressureMSB').('group') = "data_vars";
vec_data.('PressureMSB').('units') = "dbar";

vec_data.('AnaIn2MSB').('dims')  = {};
vec_data.('AnaIn2MSB').('dtype') = "uint8";
vec_data.('AnaIn2MSB').('group') = "sys";
vec_data.('AnaIn2MSB').('units') = "1";

vec_data.('PressureLSW').('dims')  = {};
vec_data.('PressureLSW').('dtype') = "uint16";
vec_data.('PressureLSW').('group') = "data_vars";
vec_data.('PressureLSW').('units') = "dbar";

vec_data.('AnaIn1').('dims') = {};
vec_data.('AnaIn1').('dtype') = "uint16";
vec_data.('AnaIn1').('group') = "sys";
vec_data.('AnaIn1').('units') = "1";

vec_data.('vel').('dims') = {3};
vec_data.('vel').('dtype') = "float32";
vec_data.('vel').('group') = "data_vars";
vec_data.('vel').('factor') = 0.001;
vec_data.('vel').('default_val') = nan;
vec_data.('vel').('units') = "m s-1";

vec_data.('amp').('dims')  = {3};
vec_data.('amp').('dtype') = "uint8";
vec_data.('amp').('group') = "data_vars";
vec_data.('amp').('units') = "dB";

vec_data.('corr').('dims')  = {3};
vec_data.('corr').('dtype') = "uint8";
vec_data.('corr').('group') = "data_vars";
vec_data.('corr').('units') = "%";
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                        vec_sysdata
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
vec_sysdata.('time').('dims') = {};
vec_sysdata.('time').('dtype') = "datetime";
vec_sysdata.('time').('group') = "coords";
vec_sysdata.('time').('default_val') = nan;
vec_sysdata.('time').('units') = "1";

vec_sysdata.('batt').('dims') = {};
vec_sysdata.('batt').('dtype') = "float32";
vec_sysdata.('batt').('group') = "data_vars";
vec_sysdata.('batt').('default_val') = nan;
vec_sysdata.('batt').('factor') = 0.1;
vec_sysdata.('batt').('units') = "V";

vec_sysdata.('c_sound').('dims') = {};
vec_sysdata.('c_sound').('dtype') = "float32";
vec_sysdata.('c_sound').('group') = "data_vars";
vec_sysdata.('c_sound').('default_val') = nan;
vec_sysdata.('c_sound').('factor') = 0.1;
vec_sysdata.('c_sound').('units') = "m s-1";

vec_sysdata.('heading').('dims') = {};
vec_sysdata.('heading').('dtype') = "float32";
vec_sysdata.('heading').('group') = "data_vars";
vec_sysdata.('heading').('default_val') = nan;
vec_sysdata.('heading').('factor') = 0.1;
vec_sysdata.('heading').('units') = "degree";

vec_sysdata.('pitch').('dims') = {};
vec_sysdata.('pitch').('dtype') = "float32";
vec_sysdata.('pitch').('group') = "data_vars";
vec_sysdata.('pitch').('default_val') = nan;
vec_sysdata.('pitch').('factor') = 0.1;
vec_sysdata.('pitch').('units') = "degree";

vec_sysdata.('roll').('dims') = {};
vec_sysdata.('roll').('dtype') = "float32";
vec_sysdata.('roll').('group') = "data_vars";
vec_sysdata.('roll').('default_val') = nan;
vec_sysdata.('roll').('factor') = 0.1;
vec_sysdata.('roll').('units') = "degree";

vec_sysdata.('temp').('dims') = {};
vec_sysdata.('temp').('dtype') = "float32";
vec_sysdata.('temp').('group') = "data_vars";
vec_sysdata.('temp').('default_val') = nan;
vec_sysdata.('temp').('factor') = 0.01;
vec_sysdata.('temp').('units') = "degree_C";

vec_sysdata.('error').('dims') = {};
vec_sysdata.('error').('dtype') = "uint8";
vec_sysdata.('error').('group') = "data_vars";
vec_sysdata.('error').('default_val') = nan;
vec_sysdata.('error').('units') = "1";

vec_sysdata.('status').('dims') = {};
vec_sysdata.('status').('dtype') = "uint8";
vec_sysdata.('status').('group') = "data_vars";
vec_sysdata.('status').('default_val') = nan;
vec_sysdata.('status').('units') = "1";

vec_sysdata.('AnaIn').('dims') = {};
vec_sysdata.('AnaIn').('dtype') = "float32";
vec_sysdata.('AnaIn').('group') = "sys";
vec_sysdata.('AnaIn').('default_val') = nan;
vec_sysdata.('AnaIn').('units') = "1";
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                        awac_profile
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
awac_profile.("time").("dims") = {};
awac_profile.("time").("dtype") = "datetime";
awac_profile.("time").("group") = "coords";
awac_profile.("time").("units") = "1";

awac_profile.("error").("dims") = {};
awac_profile.("error").("dtype") = "uint16";
awac_profile.("error").("group") = "data_vars";
awac_profile.("error").("units") = "1";

awac_profile.("AnaIn1").("dims") = {};
awac_profile.("AnaIn1").("dtype") = "float32";
awac_profile.("AnaIn1").("group") = "sys";
awac_profile.("AnaIn1").("units") = "n/a";
awac_profile.("AnaIn1").("default_val") = nan;

awac_profile.("batt").("dims") = {};
awac_profile.("batt").("dtype") = "float32";
awac_profile.("batt").("group") = "data_vars";
awac_profile.("batt").("units") = "V";
awac_profile.("batt").("factor") = 0.1;
awac_profile.("batt").("default_val") = nan;

awac_profile.("c_sound").("dims") = {};
awac_profile.("c_sound").("dtype") = "float32";
awac_profile.("c_sound").("group") = "data_vars";
awac_profile.("c_sound").("units") = "m s-1";
awac_profile.("c_sound").("factor") = 0.1;
awac_profile.("c_sound").("default_val") = nan;

awac_profile.("heading").("dims") = {};
awac_profile.("heading").("dtype") = "float32";
awac_profile.("heading").("group") = "data_vars";
awac_profile.("heading").("units") = "degree";
awac_profile.("heading").("factor") = 0.1;
awac_profile.("heading").("default_val") = nan;

awac_profile.("pitch").("dims") = {};
awac_profile.("pitch").("dtype") = "float32";
awac_profile.("pitch").("group") = "data_vars";
awac_profile.("pitch").("units") = "degree";
awac_profile.("pitch").("factor") = 0.1;
awac_profile.("pitch").("default_val") = nan;

awac_profile.("roll").("dims") = {};
awac_profile.("roll").("dtype") = "float32";
awac_profile.("roll").("group") = "data_vars";
awac_profile.("roll").("units") = "degree";
awac_profile.("roll").("factor") = 0.1;
awac_profile.("roll").("default_val") = nan;

awac_profile.("pressure").("dims") = {};
awac_profile.("pressure").("dtype") = "float32";
awac_profile.("pressure").("group") = "data_vars";
awac_profile.("pressure").("units") = "dbar";
awac_profile.("pressure").("factor") = 0.001;
awac_profile.("pressure").("default_val") = nan;

awac_profile.("status").("dims") = {};
awac_profile.("status").("dtype") = "float32";
awac_profile.("status").("group") = "data_vars";
awac_profile.("status").("units") = "1";
awac_profile.("status").("default_val") = nan;

awac_profile.("temp").("dims") = {};
awac_profile.("temp").("dtype") = "float32";
awac_profile.("temp").("group") = "data_vars";
awac_profile.("temp").("units") = "degree_C";
awac_profile.("temp").("factor") = 0.01;
awac_profile.("temp").("default_val") = nan;

awac_profile.("vel").("dims") = {3, 'nbins', 'n'};
awac_profile.("vel").("dtype") = "float32";
awac_profile.("vel").("group") = "data_vars";
awac_profile.("vel").("units") = "m s-1";
awac_profile.("vel").("factor") = 0.001;
awac_profile.("vel").("default_val") = nan;

awac_profile.("amp").("dims") = {3, 'nbins', 'n'};
awac_profile.("amp").("dtype") = "uint8";
awac_profile.("amp").("group") = "data_vars";
awac_profile.("amp").("units") = "counts";
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

end
