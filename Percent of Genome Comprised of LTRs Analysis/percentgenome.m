clc
close all
clear all
set(groot,'defaultFigurePaperPositionMode','auto')
assemblylength.neosp1 = 193032486; %from mycocosm
assemblylength.caecom1 = 165495782;
assemblylength.pirfi3 = 56455805;
assemblylength.anasp1 = 71685009;
assemblylength.gfma = 209503801;
assemblylength.neolan1 = 200974851;
assemblylength.uh31 =  84096456;

%% Find tesorter lengths
genomes = {'anasp1', 'uh31', 'pirfi3','neolan1','caecom1','gfma','neosp1'};

for j = 1:length(genomes)
    current_genome = genomes{j};
    LTRlengths_data_tesorter.(current_genome) = [];
    
    ltrdata.(current_genome) = importltrlength(['ltrdata_', current_genome, '.txt']);
    ltrname.(current_genome) = importltrname(['ltrsequences_', current_genome, '.txt']);
    ltrclassification.(current_genome) = importtesorter(['tesorterresults_', current_genome, '.txt']);
    
    LTRlength_total.(current_genome) = sum(ltrdata.(current_genome));
    
    for i = 1:size(ltrclassification.(current_genome), 1)
        idx = find(strcmp([ltrname.(current_genome){:, 1}], ltrclassification.(current_genome){i, 1}));
        LTRlengths_data_tesorter.(current_genome)(i, 1) = ltrdata.(current_genome)(idx, 1);
    end
    
    LTRpercent_tesorter.(current_genome) = 100 * sum(LTRlengths_data_tesorter.(current_genome)) /assemblylength.(current_genome);
    LTRpercent_total.(current_genome) = 100 * LTRlength_total.(current_genome) / assemblylength.(current_genome);
end


%% Plot graph - includes LTRs not classified by TEsorter
% x = categorical({'Anaeromyces robustus', 'Piromyces sp. UH3-1' ,'Piromyces finnis' ,'Neocallimastix lanati' ,'Caecomyces churrovis', 'Neocallimastix sp. GF-Ma3-1', 'Neocallimastix californiae G1'});
x = categorical({'A. robustus', 'P. sp. UH3-1' ,'P. finnis' ,'N. lanati' ,'C. churrovis', 'N. sp. Gf-Ma', 'N. californiae'});
y_totalpercent = [LTRpercent_total.anasp1 LTRpercent_total.uh31 LTRpercent_total.pirfi3 LTRpercent_total.neolan1 LTRpercent_total.caecom1 LTRpercent_total.gfma LTRpercent_total.neosp1];

figure;
bar(x,y_totalpercent,'facecolor',[0.5, 0.5, 0.5]);
set(gca,'fontweight','bold')
set(gcf,'color','w');
ylabel('Percentage of genome composed of LTR Retrotransposons', 'FontWeight','bold')
xlabel('Isolate')
set(gcf,'color','w');
saveas(gcf, 'ltrharvestbar.png')
saveas(gcf, 'ltrharvest.svg')


%% Plot graph - ONLY includes LTRs classified by TEsorter
y_tesorterpercent = [LTRpercent_tesorter.anasp1 LTRpercent_tesorter.uh31 LTRpercent_tesorter.pirfi3 LTRpercent_tesorter.neolan1 LTRpercent_tesorter.caecom1 LTRpercent_tesorter.gfma LTRpercent_tesorter.neosp1];

figure;
bar(x,y_tesorterpercent,'facecolor',[0.5, 0.5, 0.5]);
set(gca,'fontweight','bold')
set(gcf,'color','w');
ylabel('Percentage of genome composed of LTR Retrotransposons', 'FontWeight','bold')
xlabel('Isolate')
set(gcf,'color','w');
saveas(gcf, 'tesorterbar.png')
saveas(gcf, 'tesorterbar.svg')


%% Stacked bar graph with classifications


% create array with classification ltr lengths for each genome

for j = 1:length(genomes)
    current_genome = genomes{j};
    labels.(current_genome) = table2array(unique(ltrclassification.gfma(:,2)));
    
    for i = 1:length(labels.(current_genome))
        label = labels.(current_genome)(i);
        ltrclassnames.(current_genome).(label) = ltrclassification.(current_genome){find(strcmp(ltrclassification.(current_genome){:, 2}, labels.(current_genome)(i))),1};
        
        ltrclasslengths.(current_genome).(label) = [];% initialize ltr class lengths array

        for k = 1:size(ltrclassnames.(current_genome).(label), 1)
            idx = find(strcmp([ltrname.(current_genome){:, 1}], ltrclassnames.(current_genome).(label){k, 1}));
            
            ltrclasslengths.(current_genome).(label)(k, 1) = ltrdata.(current_genome)(idx, 1);
         
        end
        ltrclasspercent.(current_genome).(label) = 100*sum(ltrclasslengths.(current_genome).(label))./assemblylength.(current_genome);
    end

end

for j = 1:length(genomes)
    current_genome = genomes{j};
    gypsy(1,j) = ltrclasspercent.(current_genome).Gypsy;
    copia(1,j) = ltrclasspercent.(current_genome).Copia;
    mixture(1,j) = ltrclasspercent.(current_genome).mixture;
    Caulimoviridae(1,j) = ltrclasspercent.(current_genome).Caulimoviridae;
    unknown(1,j) = ltrclasspercent.(current_genome).unknown;

end

classifications = [gypsy; copia; mixture; unknown; Caulimoviridae];


figure;
bar(x, [classifications; y_totalpercent-y_tesorterpercent], 'stacked');
ylabel('Percentage of genome composed of LTR Retrotransposons', 'FontWeight','bold')
xlabel('Isolate')
legend('Gypsy', 'Copia', 'Mixture', 'Unknown', 'Caulimoviridae','Unclassified by TESorter')
saveas(gcf, 'stackedperc.png')
saveas(gcf, 'stackedperc.svg')
