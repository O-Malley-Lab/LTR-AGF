function f = graphallltrs(cdhit)
%this function graphs all LTRs in clusters, includes non transcribed LTRs

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
title('N. californiae')
xlabel('Cluster Size')
ylabel("Number of Clusters")
set(gcf,'color','w');
breakyaxis([80 295]);

end

