function ltrsequencesgfma = importdata_ltrharvest(filename, dataLines)
%IMPORTFILE Import data from a text file
%  LTRSEQUENCESGFMA = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  LTRSEQUENCESGFMA = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  ltrsequencesgfma = importfile("/Users/nabilalansari/Desktop/LTRs/ltrsequences_gfma.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 31-Mar-2024 19:11:37

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ">";

% Specify column names and types
opts.VariableNames = ["Var1", "scaffold_1_dbseqnr_0_1570523681"];
opts.SelectedVariableNames = "scaffold_1_dbseqnr_0_1570523681";
opts.VariableTypes = ["string", "string"];

% Specify file level properties
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "scaffold_1_dbseqnr_0_1570523681"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "scaffold_1_dbseqnr_0_1570523681"], "EmptyFieldRule", "auto");

% Import the data
ltrsequencesgfma = readtable(filename, opts);

end