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
clusteredltrs_ltrdata = rmmissing(clusteredltrs_ltrdata);

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
        while startsWith(clusteredltrs_ltrdata{d,'Clusterno_LTRName'},"scaffold")
            count = count + 1;
            d=d+1;
        end
        ltrsinclustercount(j,1)= count;
        j=j+1;
    end
end

%remove clusters that are of size 1
ltrsinclustercount(ltrsinclustercount==1) = [];

%graph the data
figure
f= histogram(ltrsinclustercount,'facecolor',[0.5, 0.5, 0.5]);
title('N. californiae')
xlabel('Cluster Size','fontweight','bold')
ylabel("Number of Clusters",'fontweight','bold')
set(gcf,'color','w');
end

