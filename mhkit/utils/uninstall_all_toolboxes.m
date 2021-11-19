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

%%
installedToolbox = struct;
installedToolbox.Name = "Mixed-Signal Blockset";
installedToolbox.Version = "2.1";
installedToolbox.Guid = "AM";
matlab.addons.toolbox.uninstallToolbox(installedToolbox);
result = 0;
%%

% toolboxes = matlab.addons.toolbox.installedToolboxes;
% addons = matlab.addons.installedAddons;
% 
% if height(addons) < 1
%     result = 0;
% else
%     % Note: this cannot be vectorized
%     result = 0;
%     for i = 1:height(addons)
%         % output = matlab.addons.toolbox.uninstallToolbox(toolboxes(i));
%         if strcmpi('matlab', addons{i,1})
%             continue;
%         else
%             matlab.addons.disableAddon(addons{i,1});
%             result = matlab.addons.isAddonEnabled(addons{i,1});
%             result = 0 | result;
%         end
%     end
% end
