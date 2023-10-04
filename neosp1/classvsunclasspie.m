
function f = classvsunclasspie(ltrclassification, ltrname)

percentages = [size(ltrclassification,1)/size(ltrname,1), (size(ltrname,1)-size(ltrclassification,1))/size(ltrname,1)]; % Two percentages: 30% and 70%
labels = {'Classified', 'Unclassified'}; % Labels for the percentages

stacked = [size(ltrclassification,1), (size(ltrname,1)-size(ltrclassification,1))]; % Two percentages: 30% and 70%
figure;
% Plot the stacked bar chart
bar(categorical({'N. Californiae'}), stacked, 'stacked');
ylabel('LTRs Identified')
xlabel('Genome')

figure;
% Create the pie chart
f = pie(percentages);

colors = [0.8500 0.3250 0.0980;0 0.4470 0.7410;.5 .5 .5;0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560];

colormap(colors);

legend(labels,'Location','northeastoutside')
set(findobj(f,'type','text'),'fontweight','bold', 'fontsize',14)



end