function reference_model = reference_model_power_matrix(device)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load the power matrix of a Reference Model wave energy converter
%
% The Reference Model Project, sponsored by the U.S. Department of
% Energy, developed open-source marine energy point designs as reference
% models to benchmark technology performance and costs, see
% https://openei.org/wiki/PRIMRE/Signature_Projects/Reference_Model.
% The power matrices are from the National Laboratory of the Rockies
% (NLR) System Advisor Model (SAM) wave energy converter library,
% https://github.com/NatLabRockies/SAM. They ship with MHKiT as
% pre-built .mat files, mhkit/wave/reference_models/<device>.mat; see the
% data_source and last_modified fields for provenance.
%
% Parameters
% ------------
%     device : string
%         "RM3" (wave point absorber), "RM5" (oscillating surge flap), or
%         "RM6" (oscillating water column)
%
% Returns
% ---------
%     reference_model : structure
%         reference_model.values : matrix
%             Device power [W], Hs rows by Te columns
%         reference_model.stat : char
%             'power'
%         reference_model.x_bins, x_edges : vector
%             Te bin centers and edges [s]
%         reference_model.y_bins, y_edges : vector
%             Hs bin centers and edges [m]
%         reference_model.device, technology_type, pto_type, description : char
%             Device metadata from the SAM library
%         reference_model.characteristic_diameter : double
%             Characteristic diameter [m]
%         reference_model.mass : double
%             Unballasted structural mass [kg]
%         reference_model.data_source : char
%             Where the matrix came from
%         reference_model.reference : char
%             Report describing the devices, Neary et al. (2014), SAND2014-9040
%         reference_model.last_modified : datetime
%             When the .mat file was generated [UTC]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    device {mustBeTextScalar, mustBeMember(device, ["RM3", "RM5", "RM6"])}
end

mat_file = fullfile(fileparts(mfilename('fullpath')), char(device) + ".mat");
loaded = load(mat_file, 'reference_model');
reference_model = loaded.reference_model;

end
