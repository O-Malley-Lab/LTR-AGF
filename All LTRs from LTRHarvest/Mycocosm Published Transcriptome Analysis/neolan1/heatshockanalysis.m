function f = heatshockanalysis(cdhit, heatshockdata)

% Set plot formatting
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultAxesFontWeight','bold')
set(groot,'DefaultAxesFontSize',12)
set(groot,'DefaultTextFontName','Arial')
set(groot,'DefaultTextFontWeight','bold')
set(groot,'DefaultTextFontSize',12)


%intialize new heatshock empty table for each time
times = {'before','t0', 't15', 't30','t45', 't60'};
for i = 1:length(times)
    heatshock.(times{i}) = table;
end

heatshock.before = heatshockdata(:,{'B1', 'B2', 'B3', 'B4'}); %save before heatshock values to new array
heatshock.t0 = heatshockdata(:,{'t0R1','t0R2','t0R3','t0R4'}); %save time 0 after heatshock values to new array
heatshock.t15 = heatshockdata(:,{'t15R1','t15R2','t15R3','t15R4'}); %save time 15 after heatshock values to new array
heatshock.t30 = heatshockdata(:,{'t30R1','t30R2','t30R3'}); %save time 30 after heatshock values to new array
heatshock.t45 = heatshockdata(:,{'t45R2','t45R3','t45R4'}); %save time 45 after heatshock values to new array
heatshock.t60 = heatshockdata(:,{'t60R1','t60R3','t60R4'}); %save time 60 after heatshock values to new array

for i = 1:length(times)
    heatshock.(times{i}) = table2array(heatshock.(times{i})).';
end


times = {'t0', 't15'};
for i = 1:length(times)
    log2fc.(times{i}) = [];
    FDR.(times{i}) = [];
    meanexp.(times{i}) = [];
    [log2fc.(times{i}),FDR.(times{i}),meanexp.(times{i})]  = DESeq2(heatshock.(times{i}),heatshock.before); 
end

times = {'t30', 't45', 't60'};
for i = 1:length(times)
    log2fc.(times{i}) = [];
    FDR.(times{i}) = [];
    meanexp.(times{i}) = [];
    [log2fc.(times{i}),FDR.(times{i}),meanexp.(times{i})]  = DESeq2(heatshock.(times{i}),heatshock.before(1:3,:)); 
end

times = {'t0','t15','t30', 't45', 't60'};

for i = 1:length(times)
    log2fc_fieldName = ['log2FC_' times{i}];
    meanexp_fieldName = ['meanExp_' times{i}];
    fdr_fieldName = ['FDR_' times{i}];

    heatshockdata(:,log2fc_fieldName) = array2table(log2fc.(times{i}));
    heatshockdata(:,fdr_fieldName) = array2table(FDR.(times{i}));
    heatshockdata(:,meanexp_fieldName) = array2table(meanexp.(times{i})');
    
end
% ,'log2FC_t15','FDR_t15','meanExp_t15','log2FC_t30','FDR_t30','meanExp_t30','log2FC_t45','FDR_t45','meanExp_t45','log2FC_t60','FDR_t60','meanExp_t60') = 


%matches clustered ltrs to heatshock values, saves heatshock values to cdhit table
for j=1:length(times)
    for i=1:size(cdhit,1)
       
        log2fc_fieldName = ['log2FC_' times{j}];
        meanexp_fieldName = ['meanExp_' times{j}];
        fdr_fieldName = ['FDR_' times{j}];
    
        if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript
            idx = find(strcmp([heatshockdata{:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
            cdhit{i,log2fc_fieldName} = heatshockdata{idx,log2fc_fieldName}; %save all that data to the cdhit table
            cdhit{i,fdr_fieldName} = heatshockdata{idx,fdr_fieldName}; %save all that data to the cdhit table
            cdhit{i,meanexp_fieldName} = heatshockdata{idx,meanexp_fieldName}; %save all that data to the cdhit table
            
        end
    end
end
% 
% for i=1:size(cdhit,1)
%     if startsWith(cdhit{i,'Transcriptome'},'Locus') %if incriment hits a transcript
%         idx = find(strcmp([heatshockdata{:,1}], cdhit{i,'Transcriptome'})); %find the index on the heatshock data which matches the transcript
%         cdhit{i,'log2fc_t15'} = heatshockdata{idx,'log2FC_t15'}; %save all that data to the cdhit table
%         cdhit{i,'FDR_t15'} = heatshockdata{idx,'FDR_t15'}; %save all that data to the cdhit table
%         cdhit{i,'meanexp_t15'} = heatshockdata{idx,'meanExp_t15'}; %save all that data to the cdhit table
%         
%     end
% end


count= 0; 
j=1;
d=1;
%while loop to count number of ltrs in a cluster and get rpkm data for each
%cluster
while d < size(cdhit,1)
    if startsWith(cdhit{d,1},"Cluster")
        count = 0;
        d=d+1;

        if d == size(cdhit,1)
            break;
        end

        while startsWith(cdhit{d,1},"scaffold")
            count = count + 1;
            d=d+1;
        end

        clustercountdata(j,1) = count;
        
        cdhit{d - count - 1,'ClusterSize'} = count;

        j=j+1;
    end
end

%put deseq output into another array so it is easier to graph later
%this is also used to get an array with just one transcript per cluster,
%this is done to make sure the data array matches the size one
for k = 1:length(times)
    log2fc_fieldName = ['log2FC_' times{k}];
    meanexp_fieldName = ['meanExp_' times{k}];
    fdr_fieldName = ['FDR_' times{k}];
    j=1;
    
    for i=1:size(cdhit,1)-1
        
        if startsWith(cdhit{i,'Clusterno_LTRName'}, 'Cluster') &&  startsWith(cdhit{i+1,'Clusterno_LTRName'}, 'scaffold') % if i = cluster header and i+1 = scafold (i.e. start of cluster) start while loop
            log2fcdata.(times{k})(j,1) = cdhit{i+1,log2fc_fieldName};
            FDRdata.(times{k})(j,1) = cdhit{i+1,fdr_fieldName};
            meanexpdata.(times{k})(j,1) = cdhit{i+1,meanexp_fieldName};
            j=j+1;
        end

    end
end

titles = {'t_0','t_{15}','t_{30}', 't_{45}', 't_{60}'};


%% plot mean expression data at all times
for i = 1:length(times)
    figure;
    scatter(meanexp.(times{i})(FDR.(times{i})>=0.1),log2fc.(times{i})(FDR.(times{i})>=0.1),'MarkerFaceColor', 'black','MarkerEdgeColor', [.7 .7 .7]);
    hold on;
    scatter(meanexp.(times{i})(FDR.(times{i})<0.1),log2fc.(times{i})(FDR.(times{i})<0.1),'MarkerFaceColor', 'red','MarkerEdgeColor', 'black');
    hold on;
    yline(-1)
    hold on;
    yline(1)
    set(gca, 'XScale', 'log')
    box on;
    title(titles{i})
    set(gcf,'color','w');
    xlabel('Mean expression','fontweight','bold')
    ylabel('Log2 fold change','fontweight','bold')
    ylim([-5 11]);
    saveas(gcf, [times{i} 'allgenes_meanexpr.png'])
    saveas(gcf, [times{i} 'allgenes_meanexpr.svg'])
end
%% plot only ltrs
for i = 1:length(times)
    figure;
    scatter(meanexpdata.(times{i})(FDRdata.(times{i})>=0.1),log2fcdata.(times{i})(FDRdata.(times{i})>=0.1),30.*clustercountdata(FDRdata.(times{i})>=0.1),'MarkerFaceColor', 'black','MarkerEdgeColor', [.7 .7 .7]);
    hold on;
    scatter(meanexpdata.(times{i})(FDRdata.(times{i})<0.1),log2fcdata.(times{i})(FDRdata.(times{i})<0.1),30.*clustercountdata(FDRdata.(times{i})<0.1),'MarkerFaceColor', 'red','MarkerEdgeColor', 'black');
    hold on;
    yline(-1)
    hold on;
    yline(1)
    set(gca, 'XScale', 'log')
    box on;
    title(titles{i})
    set(gcf,'color','w');
    xlabel('Mean expression','fontweight','bold')
    ylabel('Log2 fold change','fontweight','bold')
    saveas(gcf, [times{i} 'ltrsonly_meanexprnolim.png'])
    saveas(gcf, [times{i} 'ltrsonly_meanexprnolim.svg'])
end


%find data for stacked bar chart
x = unique(clustercountdata);
for i = 1:length(times)
    for j=2:size(x,1)
        sig_up.(times{i})(j) = numel(meanexpdata.(times{i})(FDRdata.(times{i})<0.1 & log2fcdata.(times{i})>1 & clustercountdata == x(j))); %sig up
        sig_down.(times{i})(j) = numel(meanexpdata.(times{i})(FDRdata.(times{i})<0.1 & log2fcdata.(times{i})<-1 & clustercountdata == x(j))); %significant-down
        notsig_up.(times{i})(j) = numel(meanexpdata.(times{i})(FDRdata.(times{i})>=0.1 & log2fcdata.(times{i})>1 & clustercountdata == x(j)));
        notsig_down.(times{i})(j) = numel(meanexpdata.(times{i})(FDRdata.(times{i})>=0.1 & log2fcdata.(times{i})<-1 & clustercountdata == x(j)));
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

f=cdhit;
end