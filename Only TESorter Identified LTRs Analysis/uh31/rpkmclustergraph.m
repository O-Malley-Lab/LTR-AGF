function f = rpkmclustergraph(clusteredltrs_rpkm) %creates rpkm versus cluster graph


count= 0; %count for number of LTRs in a cluster
j=1; %increment for data table - the table where data is saved
d=1; % increment for cdhit cluster file

while d < size(clusteredltrs_rpkm,1)
    if startsWith(clusteredltrs_rpkm{d,1},"Cluster")
        count = 0; %at new cluster start count at zero
        d=d+1;
           if d == size(clusteredltrs_rpkm,1)
                break; %termination statement, if code reaches end of table
           end
        while ~startsWith(clusteredltrs_rpkm{d,1},"Cluster")
            count = count + 1;
            d=d+1; %if it hits an LTR it increases count of cluster and continues
        end
        clustercount(j,1)= count; %saves cluster count to data array
        rpkm(j,1)= clusteredltrs_rpkm{d-1,2}; %saves rpkm value right next to it on data array
        j=j+1;
    end
end
%now we have a data table which includes cluster size and rpkm value next
%to it, we will use this to graph on a scatter plot

figure
scatter(clustercount,rpkm, 80, 'x','k'); %plot rpkm vs cluster size
xlabel('Cluster Size [# of LTRs]','fontweight','bold')
ylabel('Expression [rpkm]','fontweight','bold')
set(gcf,'color','w'); %set bg to white
set(gca, 'YScale', 'log')
box on %add surrounding box
ax = gca;
ax.Clipping = 'on';
xlim([1 70])


figure
f = scatter(clustercount(clustercount>1),rpkm(clustercount>1),80, 'x','k'); %plot rpkm vs cluster sizeylabel('Cluster Size','fontweight','bold')
ylabel('Expression [rpkm]','fontweight','bold')
xlabel('Cluster Size [# of LTRs]','fontweight','bold')
set(gcf,'color','w'); %set bg to white
set(gca, 'YScale', 'log')
box on %add surrounding box
ax = gca;
ax.Clipping = 'on';
xlim([2 70])


% Add a tick mark for 2 on the x-axis
xticks([2, get(gca, 'XTick')]);

% Refresh the plot
hold off

end