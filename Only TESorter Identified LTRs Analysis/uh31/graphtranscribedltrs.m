function f = graphtranscribedltrs(clusteredltrs_ltrdata)

count= 0; 
j=1;
d=1;


%enters 0 for rpkm values for cluster rows - needs to be done so the cluster number rows aren't removed
for i=1:size(clusteredltrs_ltrdata,1)
    if startsWith(clusteredltrs_ltrdata{i,1},"Cluster")
        clusteredltrs_ltrdata{i,2} = "0";
    end
end

%removes ltrs that aren't transcribed
emptyStringIndices = strcmp(clusteredltrs_ltrdata{:, 2}, "");

clusteredltrs_ltrdata = clusteredltrs_ltrdata(~emptyStringIndices, :);


%removes the cluster headers for the clusters that weren't transcribed
k=0;
l=0;
while k < size(clusteredltrs_ltrdata,1)
    l=l+1;
    k=l+1;
    if startsWith(clusteredltrs_ltrdata{l,'Clusterno_LTRName'},"Cluster") && startsWith(clusteredltrs_ltrdata{k,'Clusterno_LTRName'},"Cluster")
       clusteredltrs_ltrdata(l,:) = [];
       l=l-1;
    end
end



%while loop to count number of ltrs in a cluster
while d < size(clusteredltrs_ltrdata,1)
    if startsWith(clusteredltrs_ltrdata{d,'Clusterno_LTRName'},"Cluster")
        count = 0;
        d=d+1;
           if d == size(clusteredltrs_ltrdata,1)
                break;
            end
        while startsWith(clusteredltrs_ltrdata{d,'Clusterno_LTRName'},"scaffold") || startsWith(clusteredltrs_ltrdata{d,'Clusterno_LTRName'},"chromosome")
            count = count + 1;
            d=d+1;
             if d > size(clusteredltrs_ltrdata,1)
                ltrsinclustercount(j,1)= count;
                break;
            end
        end
        ltrsinclustercount(j,1)= count;
        j=j+1;
    end
end

%remove clusters that are of size 1
ltrsinclustercount(ltrsinclustercount==1) = [];
% 
% % Define the desired bin width
% binWidth = 2;
% 
% % Calculate the number of bins based on the bin width
% numBins = ceil((max(ltrsinclustercount) - min(ltrsinclustercount)) / binWidth);
% 
% % Calculate the custom bin edges
% binEdges = linspace(min(ltrsinclustercount), max(ltrsinclustercount), numBins + 1);

%graph the data
figure
f= histogram(ltrsinclustercount,'facecolor',[0.5, 0.5, 0.5]);
title('P. sp. UH3-1')
xlabel('Cluster Size','fontweight','bold')
ylabel("Number of Clusters",'fontweight','bold')
set(gcf,'color','w');
end

