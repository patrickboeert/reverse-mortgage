function [exists,varargout] = classinworkspace(class)
% CLASSINWORKSPACE Is a specific class currently in workspace?
%   Determines whether or not an instance of a specified class currently
%   exists in the workspace, 
%
%   input:  class name [char]
%   output: existence of instance [0/1]
%           name of the existent instance [char]

    % Load data from workspace; search for instance
    workspace = evalin('base','whos');
    found = 0;
    instance = {};
    varargout = {};
    
    % loop over all variables
    for i=min(size(workspace)):length(workspace)
        if strcmpi(workspace(i).class,class) == 1
            instance{i} = workspace(i).name;
            found = found + 1;
        end
    end
   
    % prepare output
    if found >= 1;      % whether or not instance has been found
        exists = 1;
        if nargout == 2
            for i=1:found
                varargout(i) = {instance{i}}; % assign name of instance in workspace if demanded
            end
        end
    else exists = 0;
    end  
    
end % CLASSINWORKSPACE function