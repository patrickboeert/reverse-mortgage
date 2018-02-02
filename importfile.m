function importfile(filetoread)
% IMPORTFILE  Loads data from file into workspace
%   Variables are labelled according to column

    % Import the file
    newdata = importdata(filetoread);

    % Break the data up into a new structure with one field per column.
    colheaders = genvarname(newdata.colheaders);
    for i = 1:length(colheaders)
        dataByColumn1.(colheaders{i}) = newdata.data(:, i);
    end

    % Create new variables in the base workspace from those fields.
    vars = fieldnames(dataByColumn1);
    for i = 1:length(vars)
        assignin('caller', vars{i}, dataByColumn1.(vars{i}));
    end
    
end