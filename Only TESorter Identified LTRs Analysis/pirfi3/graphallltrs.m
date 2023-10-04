function f = graphallltrs(cdhit, ltrclassification)
%this function graphs all LTRs in clusters, includes non transcribed LTRs

%matches classification to cluster
for i=1:size(cdhit,1)
    idx = find(strcmp([ltrclassification{:,1}], cdhit{i,1}));
    if (isempty(idx)) %this is for the rows that say Cluster #
        cdhit{i,'classification'} = "0";
    else
        cdhit{i,'classification'}= ltrclassification{idx,2};
    end
end



cdhit((strcmp(cdhit.classification,'0') & startsWith(cdhit.Clusterno_LTRName,'scaffold')), :) = [];

%removes the cluster headers for the clusters that weren't transcribed
k=0;
l=0;
while k < size(cdhit,1)
    l=l+1;
    k=l+1;
    if startsWith(cdhit{l,'Clusterno_LTRName'},"Cluster") && startsWith(cdhit{k,'Clusterno_LTRName'},"Cluster")
       cdhit(l,:) = [];
       l=l-1;
    end
end

count= 0; 
j=1;
d=1;

while d < size(cdhit,1)
    if startsWith(cdhit{d,'Clusterno_LTRName'},"Cluster")
        count = 0;
        d=d+1;

           if d == size(cdhit,1)
                break;
           end

        while startsWith(cdhit{d,'Clusterno_LTRName'},"scaffold")
            count = count + 1;
            d=d+1;
        end
        
        data(j,1)= count;
        j=j+1;
    
    end
end


%remove clusters that are of size 1
data(data==1) = [];

%graph the data
figure
f= histogram(data,'facecolor',[0.5, 0.5, 0.5]);
title('P. finnis')
xlabel('Cluster Size')
ylabel("Number of Clusters")
set(gcf,'color','w');

end

