function f = clustercount(clusteredltrs_ltrdata)

%add alignment length column
ClusterSize = zeros(height(clusteredltrs_ltrdata),1);
clusteredltrs_ltrdata = addvars(clusteredltrs_ltrdata,ClusterSize);

%while loop to count number of ltrs in a cluster and get rpkm data for each
%cluster
d=1;
j=1;
while d <= size(clusteredltrs_ltrdata,1)
    if startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
        count = 0;
        d=d+1;

        if d == size(clusteredltrs_ltrdata,1) && ~startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            clusteredltrs_ltrdata{d - 1,'ClusterSize'} = 1;
            break;
        end
        
        if d == size(clusteredltrs_ltrdata,1) && startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            break;
        end

        while ~startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            count = count + 1;
            d=d+1;
            if d == size(clusteredltrs_ltrdata,1)
                clusteredltrs_ltrdata{d - 1,'ClusterSize'} = count;
                break;
            end
        end

        clustercountdata(j,1) = count;
        
        clusteredltrs_ltrdata{d - count - 1,'ClusterSize'} = count;

        j=j+1;


        if d == size(clusteredltrs_ltrdata,1) && ~startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            clusteredltrs_ltrdata{d - 1,'ClusterSize'} = 1;
            break;
        end
        
        if d == size(clusteredltrs_ltrdata,1) && startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            break;
        end

    end
end

count= 0; 
j=1;
d=1;
%save cluster data size throughout the rows in cdhit table
while d < size(clusteredltrs_ltrdata,1)
    if startsWith(clusteredltrs_ltrdata{d,1},"Cluster") 
        count = clusteredltrs_ltrdata{d,'ClusterSize'} ;
        d=d+1;

        if d == size(clusteredltrs_ltrdata,1) && ~startsWith(clusteredltrs_ltrdata{d,1},"Cluster") 
            clusteredltrs_ltrdata{d,'ClusterSize'} = count;
            break;
        end

        if d == size(clusteredltrs_ltrdata,1) && startsWith(clusteredltrs_ltrdata{d,1},"Cluster") 
            break;
        end

        while ~startsWith(clusteredltrs_ltrdata{d,1},"Cluster")
            clusteredltrs_ltrdata{d,'ClusterSize'} = count;
            d=d+1;
        end

    end
end

m = size(clusteredltrs_ltrdata,1);
while startsWith(clusteredltrs_ltrdata{m,'Clusterno_LTRName'}, "Cluster")
    clusteredltrs_ltrdata(size(clusteredltrs_ltrdata,1),:)=[];
    m = size(clusteredltrs_ltrdata,1);
end


for i=1:size(clusteredltrs_ltrdata,1)
    if startsWith(clusteredltrs_ltrdata{i,1},"Cluster") 
        clusteredltrs_ltrdata{i,2:end} =0;
    end
end

f =clusteredltrs_ltrdata;
end