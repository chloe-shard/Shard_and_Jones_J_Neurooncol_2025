

%% Use this code to Import Slide_info tab

opts = spreadsheetImportOptions("NumVariables", 2);
% Specify sheet and range
opts.Sheet = "slides_info";
opts.DataRange = "A2:B100";

% Specify column names and types
opts.VariableNames = ["slide_name", "Data_range"];
opts.VariableTypes = ["string", "string"];

% Specify variable properties
opts = setvaropts(opts, ["slide_name", "Data_range"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["slide_name", "Data_range"], "EmptyFieldRule", "auto");

% Import the data
Slide_Info = readmatrix(strcat(Source_dir,"\",Source_file), opts, "UseExcel", false);


%% Clear temporary variables
clear opts


