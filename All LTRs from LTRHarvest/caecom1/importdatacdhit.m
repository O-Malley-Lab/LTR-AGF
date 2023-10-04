function resultsdefault95 = importdatacdhit(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  RESULTSDEFAULT95 = IMPORTFILE3(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  RESULTSDEFAULT95 = IMPORTFILE3(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  resultsdefault95 = importfile3("/Users/nabilalansari/Desktop/cd-hit code/resultsdefault95.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 30-Sep-2022 22:50:53

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ["...", ">"];

% Specify column names and types
opts.VariableNames = ["Var1", "Clusterno_LTRName", "Var3"];
opts.SelectedVariableNames = "Clusterno_LTRName";
opts.VariableTypes = ["string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Clusterno_LTRName", "Var3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Clusterno_LTRName", "Var3"], "EmptyFieldRule", "auto");

% Import the data
resultsdefault95 = readtable(filename, opts);

end