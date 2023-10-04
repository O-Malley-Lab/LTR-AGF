function G1substraterawexpectedcounts = importsubstrateraw(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  G1SUBSTRATERAWEXPECTEDCOUNTS = IMPORTFILE3(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the data as a
%  table.
%
%  G1SUBSTRATERAWEXPECTEDCOUNTS = IMPORTFILE3(FILE, DATALINES) reads
%  data for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  G1substraterawexpectedcounts = importfile3("/Users/nabilalansari/Desktop/Research/june 18/Genomes/neosp1/G1_substrate_raw_expectedcounts.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 22-Jun-2023 21:44:02

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 24);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Transcript", "A1", "A2", "A3", "AS1", "AS2", "AS3", "CB1", "CB2", "CB3", "CS1", "CS2", "CS3", "G1", "G2", "G3", "M1", "M2", "M3", "RCG1", "RCG2", "RCG3", "SG1", "SG2"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Transcript", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Transcript", "EmptyFieldRule", "auto");

% Import the data
G1substraterawexpectedcounts = readtable(filename, opts);

end