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


%% Find tesorter lengths and create histogram
genomes = {'neosp1', 'anasp1', 'pirfi3','caecom1','gfma','neolan1','uh31'};
names = {'N. californiae', 'A. robustus', 'P. finnis', 'C. churrovis', 'N. sp. Gf-Ma', 'N. lanati', 'P. sp. UH3-1'};

for j = 1:length(genomes)
    current_genome = genomes{j};%set current genome
    LTRlengths_data_tesorter.(current_genome) = []; %initialize matrix
    
    ltrdata.(current_genome) = importltrlength(['ltrdata_', current_genome, '.txt']); %import ltr lengths
    ltrname.(current_genome) = importltrname(['ltrsequences_', current_genome, '.txt']);%import ltr names
    ltrclassification.(current_genome) = importtesorter(['tesorterresults_', current_genome, '.txt']); %import tesorter results
    
    for i = 1:size(ltrclassification.(current_genome), 1)
        idx = find(strcmp([ltrname.(current_genome){:, 1}], ltrclassification.(current_genome){i, 1})); %find index where ltrname = classification
        LTRlengths_data_tesorter.(current_genome)(i, 1) = ltrdata.(current_genome)(idx, 1); %save ret length to new matrix
    end

    %create figure
    figure;
    histogram(LTRlengths_data_tesorter.(current_genome),'normalization','probability','facecolor',[0.5, 0.5, 0.5]);
    title(names{j})
    xlabel('Sequence Length','fontweight','bold')
    ylabel("Frequency",'fontweight','bold')
    
    %save figure
    saveas(gcf, [current_genome '_retlengths.png'])
end

