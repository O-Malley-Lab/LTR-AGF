clc
close all
clear all
% Set plot formatting
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultAxesFontWeight','bold')
set(groot,'DefaultAxesFontSize',12)
set(groot,'DefaultTextFontName','Arial')
set(groot,'DefaultTextFontWeight','bold')
set(groot,'DefaultTextFontSize',12)

genomes = {'anasp1', 'uh31', 'pirfi3','neolan1','caecom1','gfma','neosp1'};

classified = [];
for j = 1:length(genomes)
    current_genome = genomes{j};
    LTRlengths_data_tesorter.(current_genome) = [];
    
    ltrdata.(current_genome) = importltrlength(['ltrdata_', current_genome, '.txt']);
    ltrname.(current_genome) = importltrname(['ltrsequences_', current_genome, '.txt']);
    ltrclassification.(current_genome) = importtesorter(['tesorterresults_', current_genome, '.txt']);
    classified(1,j) = [size(ltrclassification.(current_genome),1)];
    unclassified(1,j) = size(ltrname.(current_genome),1)-size(ltrclassification.(current_genome),1);
end

x = categorical({'Anaeromyces robustus', 'Piromyces sp. UH3-1' ,'Piromyces finnis' ,'Neocallimastix lanati' ,'Caecomyces churrovis', 'Neocallimastix giraffae', 'Neocallimastix californiae'});

%% classified vs unclassified
stacked = [size(ltrclassification,1), (size(ltrname,1)-size(ltrclassification,1))]; 
figure;
% Plot the stacked bar chart
bar(x, [classified; unclassified], 'stacked');
ylabel('Number of LTR Retrotransposons')
xlabel('Isolate')
legend('Classified by TESorter','Unclassified by TESorter')
saveas(gcf, 'classvsunclass.png')
ax = gca;

% Set the tick labels to italic
ax.XTickLabel = cellfun(@(x) ['\it' x], ax.XTickLabel, 'UniformOutput', false);
%% Classification Plot, does not include unclassified by tesorter


for j = 1:length(genomes)
    current_genome = genomes{j};
    labels.(current_genome) = table2array(unique(ltrclassification.gfma(:,2)));
    for i = 1:length(labels.(current_genome))
        label = labels.(current_genome)(i);
        counts.(current_genome).(label) = sum(strcmp(ltrclassification.(current_genome){:, 2}, labels.(current_genome)(i)));
    end
end

for j = 1:length(genomes)
    current_genome = genomes{j};
    gypsy(1,j) = counts.(current_genome).Gypsy;
    copia(1,j) = counts.(current_genome).Copia;
    mixture(1,j) = counts.(current_genome).mixture;
    Caulimoviridae(1,j) = counts.(current_genome).Caulimoviridae;
    unknown(1,j) = counts.(current_genome).unknown;

end

classifications = [gypsy; copia; mixture; unknown; Caulimoviridae];

figure;
bar(x, classifications, 'stacked');
ylabel('Number of LTR Retrotransposons')
xlabel('Genome')
legend('Gypsy', 'Copia', 'Mixture', 'Unknown', 'Caulimoviridae')
saveas(gcf, 'ltrnumbersclass.png')
ax = gca;

% Set the tick labels to italic
ax.XTickLabel = cellfun(@(x) ['\it' x], ax.XTickLabel, 'UniformOutput', false);
%% percent classified

for j = 1:length(genomes)
    current_genome = genomes{j};
    ltrsize(1,j) = size(ltrclassification.(current_genome),1);
    
end

percentclass = 100.*classifications./ltrsize;


figure;
bar(x, percentclass, 'stacked');
ylabel('Percentage of LTR Retrotransposons')
xlabel('Isolate')
legend('Gypsy', 'Copia', 'Mixture', 'Unknown', 'Caulimoviridae')
saveas(gcf, 'ltrpercentclass.png')
ax = gca;

% Set the tick labels to italic
ax.XTickLabel = cellfun(@(x) ['\it' x], ax.XTickLabel, 'UniformOutput', false);
%% classified vs unclassified - includes classifications
stacked = [size(ltrclassification,1), (size(ltrname,1)-size(ltrclassification,1))]; 
figure;
% Plot the stacked bar chart
bar(x, [classifications; unclassified], 'stacked');
ylabel('Number of LTR Retrotransposons')
xlabel('Isolate')
legend('Gypsy', 'Copia', 'Mixture', 'Unknown', 'Caulimoviridae','Unclassified by TESorter')
saveas(gcf, 'classvsunclasswclassifications.png')
ax = gca;

% Set the tick labels to italic
ax.XTickLabel = cellfun(@(x) ['\it' x], ax.XTickLabel, 'UniformOutput', false);