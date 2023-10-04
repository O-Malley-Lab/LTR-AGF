function Aresults = importdeseqresCB(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  ARESULTS = IMPORTFILE3(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  ARESULTS = IMPORTFILE3(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  Aresults = importfile3("/Users/nabilalansari/Desktop/Research/June 11/Genomes/neosp1/Aresults.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 12-Jun-2023 21:33:56

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Transcriptome", "meanExp_CB", "log2FC_CB", "Var4", "Var5", "Var6", "FDR_CB"];
opts.SelectedVariableNames = ["Transcriptome","meanExp_CB", "log2FC_CB", "FDR_CB"];
opts.VariableTypes = ["string", "double", "double", "string", "string", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Transcriptome", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Transcriptome", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");

% Import the data
Aresults = readtable(filename, opts);

end