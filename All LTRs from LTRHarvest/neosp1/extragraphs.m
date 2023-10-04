function f = extragraphs(clusteredltrs_rpkm)

%What does histogram of LTR lengths look like for those that passed LTRdigest/TEsorter?

%What does heat shock data look like just for LTRs that passed TEsorter?

k=1;

for i=1:size(clusteredltrs_rpkm, 1)-1
    if startsWith(clusteredltrs_rpkm{i,'Clusterno_LTRName'}, 'Cluster') && startsWith(clusteredltrs_rpkm{i+1,'Clusterno_LTRName'}, 'scaffold')
        classifcations(k,1) = clusteredltrs_rpkm{i+1,'classification'};
        clustercountdata(k,1) = clusteredltrs_rpkm{i+1,'ClusterSize'};
        k=k+1;
    end
end

x = unique(clustercountdata);

for j=1:size(x,1)
    gypsy(j) = numel(classifcations(classifcations == "Gypsy" & clustercountdata == x(j))); %gypsy
    copia(j) = numel(classifcations(classifcations == "Copia" & clustercountdata == x(j))); %gypsy
    mixture(j) = numel(classifcations(classifcations == "mixture" & clustercountdata == x(j))); %gypsy
    unknown(j) = numel(classifcations(classifcations == "unknown" & clustercountdata == x(j))); %gypsy
    notclassified(j) = numel(classifcations(classifcations == "0" & clustercountdata == x(j))); %gypsy
end

%define colors for the stacked bar chart
stack1Color = 'red'; % significant upregulated
stack2Color = 'blue'; %  significant-down
stack3Color = 'green'; % notsig-up
stack4Color = 'yellow'; % notsig-up
stack5Color = [0.7 0.7 0.7]; %notsig-down

%plot stacked bar chart for each time point
figure;
h = bar(x, [gypsy' copia'  mixture' unknown' notclassified'], 'stacked');
%set colors
h(1).FaceColor = stack1Color;
h(2).FaceColor = stack2Color;
h(3).FaceColor = stack3Color;
h(4).FaceColor = stack4Color;
h(5).FaceColor = stack5Color;
xlabel('Cluster Size')
ylabel('Number of Classifications')
%add legend
legend('Gypsy', 'Copia','Mixture','Unknown','Not Classified by TESorter');
saveas(gcf,  'class.png')
saveas(gcf, 'class.svg')

breakyaxis([1000 2600]);
legend('Gypsy', 'Copia','Mixture','Unknown','Not Classified by TESorter');

total = gypsy + copia +mixture +notclassified +unknown;
figure;
h = bar(x, [(100.*gypsy./total)' (100.*copia./total)'  (100.*mixture./total)' (100.*unknown./total)' (100.*notclassified./total)'], 'stacked');
%set colors
h(1).FaceColor = stack1Color;
h(2).FaceColor = stack2Color;
h(3).FaceColor = stack3Color;
h(4).FaceColor = stack4Color;
h(5).FaceColor = stack5Color;
ylim([0 100])
xlabel('Cluster Size')
ylabel('Percentage of Classifications')
%add legend
legend('Gypsy', 'Copia','Mixture','Unknown','Not Classified by TESorter');
saveas(gcf,  'class.png')
saveas(gcf, 'class.svg')

legend('Gypsy', 'Copia','Mixture','Unknown','Not Classified by TESorter');
f=1;
end
