function f = heatshockanalysis(cdhit) %creates rpkm versus cluster graph

heatshockraw = importheatshockrawdata('G1_heatshock_raw_expectedcountsWithB3.txt');

% Set plot formatting
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultAxesFontWeight','bold')
set(groot,'DefaultAxesFontSize',12)
set(groot,'DefaultTextFontName','Arial')
set(groot,'DefaultTextFontWeight','bold')
set(groot,'DefaultTextFontSize',12)


heatshockdata.t0 = importdeseqrest0('t0results.csv');
heatshockdata.t15 = importdeseqrest15('t15results.csv');
heatshockdata.t30 = importdeseqrest30('t30results.csv');
heatshockdata.t45 = importdeseqrest45('t45results.csv');
heatshockdata.t60 = importdeseqrest60('t60results.csv');



%create new heatshock empty table for each time
times = {'before','t0', 't15', 't30','t45', 't60'};

for i = 1:length(times)
    heatshock.(times{i}) = table;
end

times = {'t0','t15','t30','t45','t60'};

k=1;
%matches clustered ltrs to heatshock values, saves heatshock values to cdhit table
for j=1:length(times)
    log2fc_fieldName = ['log2FC_' times{j}];
    meanexp_fieldName = ['meanExp_' times{j}];
    fdr_fieldName = ['FDR_' times{j}];
    
    for i=1:size(cdhit,1)
        if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript 
            idx = find(strcmp([heatshockdata.(times{j}){:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
            
            if isempty(idx)  
            
            else
                cdhit{i,log2fc_fieldName} = heatshockdata.(times{j}){idx,log2fc_fieldName}; %save all that data to the cdhit table
                cdhit{i,fdr_fieldName} = heatshockdata.(times{j}){idx,fdr_fieldName}; %save all that data to the cdhit table
                cdhit{i,meanexp_fieldName} = heatshockdata.(times{j}){idx,meanexp_fieldName}; %save all that data to the cdhit table
            end
        end
    end
end


%going to try to save only data from one each cluster
k=1;
for j=1:length(times)
     log2fc_fieldName = ['log2FC_' times{j}];
     meanexp_fieldName = ['meanExp_' times{j}];
     fdr_fieldName = ['FDR_' times{j}];
     k=1;
    for i=1:size(cdhit,1)-1
        if startsWith(cdhit{i,1},"Cluster") && cdhit{i+1,log2fc_fieldName} ~= 0
            log2fcdata.(times{j})(k,1) = cdhit{i+1,log2fc_fieldName};
            meanexprdata.(times{j})(k,1) = cdhit{i+1,meanexp_fieldName};
            FDRdata.(times{j})(k,1) = cdhit{i+1,fdr_fieldName};
            clustercountdata.(times{j})(k,1) = cdhit{i+1,'ClusterSize'};
            k=k+1;
        end
    end
end

titles= {'t_0','t_{15}','t_{30}','t_{45}','t_{60}'};
for j=1:length(times)
figure;
   scatter(meanexprdata.(times{j})(FDRdata.(times{j})>=0.05),log2fcdata.(times{j})(FDRdata.(times{j})>=0.05),30.*clustercountdata.(times{j})(FDRdata.(times{j})>=0.05),'MarkerFaceColor', 'black','MarkerEdgeColor', [.7 .7 .7]);
   hold on;
   scatter(meanexprdata.(times{j})(FDRdata.(times{j})<0.05),log2fcdata.(times{j})(FDRdata.(times{j})<0.05),30.*clustercountdata.(times{j})(FDRdata.(times{j})<0.05),'MarkerFaceColor', 'red','MarkerEdgeColor', 'black');
   set(gca, 'XScale', 'log')
   xlabel('Mean Expression')
   ylabel('log2FC')
   title(titles{j})
   box on
   yline(0)
   saveas(gcf, [times{j} 'ltrsonly.png'])
    saveas(gcf, [times{j} 'ltrsonly.svg'])
end





%find data for stacked bar chart
x = unique(clustercountdata.t0);
for i = 1:length(times)
    
    for j=2:size(x,1)
        sig_up.(times{i})(j) = numel(meanexprdata.(times{i})(FDRdata.(times{i})<0.05 & log2fcdata.(times{i})>0 & clustercountdata.(times{i}) == x(j))); %sig up
        sig_down.(times{i})(j) = numel(meanexprdata.(times{i})(FDRdata.(times{i})<0.05 & log2fcdata.(times{i})<0 & clustercountdata.(times{i}) == x(j))); %significant-down
        notsig_up.(times{i})(j) = numel(meanexprdata.(times{i})(FDRdata.(times{i})>=0.05 & log2fcdata.(times{i})>0 & clustercountdata.(times{i}) == x(j)));
        notsig_down.(times{i})(j) = numel(meanexprdata.(times{i})(FDRdata.(times{i})>=0.05 & log2fcdata.(times{i})<0 & clustercountdata.(times{i}) == x(j)));
    end
end

%define colors for the stacked bar chart
stack1Color = 'red'; % significant upregulated
stack2Color = [1 0.5 0.5]; %  significant-down
stack3Color = 'black'; % notsig-up
stack4Color = [0.7 0.7 0.7]; %notsig-down

%plot stacked bar chart for each time point
for i = 1:length(times)
    figure;
    h = bar(x, [sig_up.(times{i})' sig_down.(times{i})'  notsig_up.(times{i})' notsig_down.(times{i})'], 'stacked');
    %set colors
    h(1).FaceColor = stack1Color;
    h(2).FaceColor = stack2Color;
    h(3).FaceColor = stack3Color;
    h(4).FaceColor = stack4Color;
    xlabel('Cluster Size')
    ylabel('Number of Clusters')
    title(titles{i})
    %add legend
    legend('Significant & Overexpressed', 'Significant & Underexpressed','Not Significant & Overexpressed','Not Significant & Underexpressed');
    saveas(gcf, [times{i} 'stackedbarheatshock.png'])
    saveas(gcf, [times{i} 'stackedbarheatshock.svg'])
end

columnlabels = heatshockraw.Properties.VariableNames(2:end);

%matches clustered ltrs to heatshock values, saves heatshock values to cdhit table
cdhitend = size(cdhit,2);
for i=1:size(cdhit,1)
    if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript
        idx = find(strcmp([heatshockraw{:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
        newColumns = heatshockraw{idx, 2:22}; % Save all that data to the newColumns variable
        cdhit{i, cdhitend+1:cdhitend+size(newColumns, 2)} = newColumns; % Add the newColumns at the end of the cdhit table
    end
end
cdhit.Properties.VariableNames(cdhitend+1:cdhitend+size(newColumns, 2))= columnlabels;
f=cdhit;
end