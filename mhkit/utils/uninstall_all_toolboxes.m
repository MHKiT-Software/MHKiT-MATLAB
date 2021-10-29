function result=uninstall_all_toolboxes()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Uninstall all currently installed toolboxes
%     
% Parameters
% ------------
%     none 
%
% Returns
% ---------
%     result: boolean
%         0 = success, 1 = failure, not all toolboxes uninstalled
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toolboxes = matlab.addons.toolbox.installedToolboxes;

if height(toolboxes) < 1
    result = 0;
else
    % Note: this cannot be vectorized
    result = 0;
    for i = 1:height(toolboxes)
        try
            % output = matlab.addons.toolbox.uninstallToolbox(toolboxes(i));
            output = matlab.addons.disableAddon(toolboxes(i));
        catch
            result = 1;
            return;
        end

        if isempty(output)
            result = 0 | result;
        else
            result = 1;
        end
    end
end
