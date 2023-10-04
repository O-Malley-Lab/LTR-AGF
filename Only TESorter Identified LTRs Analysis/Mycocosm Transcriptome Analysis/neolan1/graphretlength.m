function f=graphretlength(ltrlengths)
%% creates histogram with retrotransposon lengths identified by LTR Harvest%%

%NOTE: this includes non transcribed LTRs 

ltrlengths= table2array(ltrlengths);

figure
f = histogram(ltrlengths(:,1),'facecolor',[0.5, 0.5, 0.5],'normalization','probability');
title('Histogram - LTR Retrotransposon Lengths')
xlabel('LTR Retrotransposon Lengths','fontweight','bold')
ylabel("Frequency",'fontweight','bold')
set(gcf,'color','w');
title('A. robustus')

end