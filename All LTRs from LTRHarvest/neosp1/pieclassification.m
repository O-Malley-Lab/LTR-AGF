function f= pieclassification(ltrclassification)

%make categoricl array for pie chart
X = categorical(table2array(ltrclassification(:,2)));
%labels
labels = table2array(unique(ltrclassification(:,2)));

%create pie chart and legend
figure
f= pie(X);
set(gcf,'color','w');
colors = [0.8500 0.3250 0.0980;0 0.4470 0.7410;.5 .5 .5;0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560];
legend(labels,'Location','northeastoutside')
set(findobj(f,'type','text'),'fontweight','bold', 'fontsize',14)

patches = findobj(f, 'Type', 'patch');
for i = 1:size(labels,1)
    set(patches(i), 'FaceColor', colors(i,:))
end


%delete slices with small percentages
del = findobj(f,'Type','Text');
isSmall = startsWith({del.String}, 'u');
isSmall2 = startsWith({del.String}, 'm');
delete(del(isSmall));
delete(del(isSmall2));


end
