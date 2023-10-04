function f = substrate(cdhit) %creates rpkm versus cluster graph


% Set plot formatting
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultAxesFontWeight','bold')
set(groot,'DefaultAxesFontSize',12)
set(groot,'DefaultTextFontName','Arial')
set(groot,'DefaultTextFontWeight','bold')
set(groot,'DefaultTextFontSize',12)
substrateraw= importsubstrateraw('G1_substrate_raw_expectedcounts.txt');

substratedata.A = importdeseqresA('AvsGresults_fromall.csv');
substratedata.AS = importdeseqresAS('ASvsGresults_fromall.csv');
substratedata.CB = importdeseqresCB('CBvsGresults_fromall.csv');
substratedata.CS = importdeseqresCS('CSvsGresults_fromall.csv');
substratedata.SG = importdeseqresSG('SGvsGresults_fromall.csv');
substratedata.RCG = importdeseqresRCG('RCGvsGresults_fromall.csv');
substratedata.M = importdeseqresM('MvsGresults_fromall.csv');

%create new heatshock empty table for each time
conditions = {'A','AS','CB','CS','SG','RCG','M'};

k=1;
%matches clustered ltrs to heatshock values, saves heatshock values to cdhit table
for j=1:length(conditions)
    log2fc_fieldName = ['log2FC_' conditions{j}];
    meanexp_fieldName = ['meanExp_' conditions{j}];
    fdr_fieldName = ['FDR_' conditions{j}];
    
    for i=1:size(cdhit,1)
        if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript 
            idx = find(strcmp([substratedata.(conditions{j}){:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
            
            if isempty(idx)  
            
            else
                cdhit{i,log2fc_fieldName} = substratedata.(conditions{j}){idx,log2fc_fieldName}; %save all that data to the cdhit table
                cdhit{i,fdr_fieldName} = substratedata.(conditions{j}){idx,fdr_fieldName}; %save all that data to the cdhit table
                cdhit{i,meanexp_fieldName} = substratedata.(conditions{j}){idx,meanexp_fieldName}; %save all that data to the cdhit table
            end
        end
    end
end

%going to try to save only data from one each cluster
k=1;
for j=1:length(conditions)
     log2fc_fieldName = ['log2FC_' conditions{j}];
     meanexp_fieldName = ['meanExp_' conditions{j}];
     fdr_fieldName = ['FDR_' conditions{j}];
     k=1;
    for i=1:size(cdhit,1)-1
        if startsWith(cdhit{i,1},"Cluster") && cdhit{i+1,log2fc_fieldName} ~= 0
            log2fcdata.(conditions{j})(k,1) = cdhit{i+1,log2fc_fieldName};
            meanexprdata.(conditions{j})(k,1) = cdhit{i+1,meanexp_fieldName};
            FDRdata.(conditions{j})(k,1) = cdhit{i+1,fdr_fieldName};
            clustercountdata.(conditions{j})(k,1) = cdhit{i+1,'ClusterSize'};
            k=k+1;
        end
    end
end


for j=1:length(conditions)
figure;
   scatter(meanexprdata.(conditions{j})(FDRdata.(conditions{j})>=0.05),log2fcdata.(conditions{j})(FDRdata.(conditions{j})>=0.05),30.*clustercountdata.(conditions{j})(FDRdata.(conditions{j})>=0.05),'MarkerFaceColor', 'black','MarkerEdgeColor', [.7 .7 .7]);
   hold on;
   scatter(meanexprdata.(conditions{j})(FDRdata.(conditions{j})<0.05),log2fcdata.(conditions{j})(FDRdata.(conditions{j})<0.05),30.*clustercountdata.(conditions{j})(FDRdata.(conditions{j})<0.05),'MarkerFaceColor', 'red','MarkerEdgeColor', 'black');
   set(gca, 'XScale', 'log')
   xlabel('Mean Expression')
   ylabel('log2FC')
   title(conditions{j})
   box on
   yline(0)
   saveas(gcf, [conditions{j} 'ltrsonly.png'])
    saveas(gcf, [conditions{j} 'ltrsonly.svg'])
end



%find data for stacked bar chart

for i = 1:length(conditions)
   x.(conditions{i}) = unique(clustercountdata.(conditions{i}));
    for j=2:size(x.(conditions{i}),1)
        sig_up.(conditions{i})(j) = numel(meanexprdata.(conditions{i})(FDRdata.(conditions{i})<0.05 & log2fcdata.(conditions{i})>0 & clustercountdata.(conditions{i}) == x.(conditions{i})(j))); %sig up
        sig_down.(conditions{i})(j) = numel(meanexprdata.(conditions{i})(FDRdata.(conditions{i})<0.05 & log2fcdata.(conditions{i})<0 & clustercountdata.(conditions{i}) == x.(conditions{i})(j))); %significant-down
        notsig_up.(conditions{i})(j) = numel(meanexprdata.(conditions{i})(FDRdata.(conditions{i})>=0.05 & log2fcdata.(conditions{i})>0 & clustercountdata.(conditions{i}) == x.(conditions{i})(j)));
        notsig_down.(conditions{i})(j) = numel(meanexprdata.(conditions{i})(FDRdata.(conditions{i})>=0.05 & log2fcdata.(conditions{i})<0 & clustercountdata.(conditions{i}) == x.(conditions{i})(j)));
    end
end

%define colors for the stacked bar chart
stack1Color = 'red'; % significant upregulated
stack2Color = [1 0.5 0.5]; %  significant-down
stack3Color = 'black'; % notsig-up
stack4Color = [0.7 0.7 0.7]; %notsig-down

%plot stacked bar chart for each time point
for i = 1:length(conditions)
    figure;
    h = bar(x.(conditions{i}), [sig_up.(conditions{i})' sig_down.(conditions{i})'  notsig_up.(conditions{i})' notsig_down.(conditions{i})'], 'stacked');
    %set colors
    xticks('auto')
    h(1).FaceColor = stack1Color;
    h(2).FaceColor = stack2Color;
    h(3).FaceColor = stack3Color;
    h(4).FaceColor = stack4Color;
    xlabel('Cluster Size')
    ylabel('Number of Clusters')
    title(conditions{i})
    %add legend
    legend('Significant & Overexpressed', 'Significant & Underexpressed','Not Significant & Overexpressed','Not Significant & Underexpressed');
    saveas(gcf, [conditions{i} 'stackedbarsubstrate.png'])
    saveas(gcf, [conditions{i} 'stackedbarsubstrate.svg'])
    
end

columnlabels = substrateraw.Properties.VariableNames(2:end);

%matches clustered ltrs to heatshock values, saves heatshock values to cdhit table
cdhitend = size(cdhit,2);
for i=1:size(cdhit,1)
    if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript
        idx = find(strcmp([substrateraw{:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
        newColumns = substrateraw{idx, 2:end}; % Save all that data to the newColumns variable
        cdhit{i, cdhitend+1:cdhitend+size(newColumns, 2)} = newColumns; % Add the newColumns at the end of the cdhit table
    end
end
cdhit.Properties.VariableNames(cdhitend+1:cdhitend+size(newColumns, 2))= columnlabels;

f=cdhit;
end