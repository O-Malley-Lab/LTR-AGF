close all
clear all
clc
set(groot, 'DefaultAxesFontSize', 14);

%% importing data

%import LTR name for matching from LTRHarvest
ltrname = importltrname('ltrsequences_neosp1.txt');

%import LTR data from LTRHarvest
ltrlengths = importltrlength('ltrdata_neosp1.txt');

%Import BLAST transcriptome matches
BLASTresults = importdatablast('70cov90ident_dbltrquerytranscriptome_G1.txt');
rpkmresults = importrpkm('70cov90ident_dbltrquerytranscriptome_G1.txt');

%Import clustered LTRs
cdhit = importdatacdhit('cdhitresultsaSaL70_neosp1.txt');

%Import LTR classification from TEsorter
ltrclassification = importtesorter('tesorterresults_neosp1.txt');

%Import BLAST results from promoters

BLASTresults_promoters = importdata_promoters('shortened_querypromoter_dbneosp1.txt');


%% function to match clustered ltr retrotransposons to rpkm values
disp('Matching rpkm to LTRs in clusters')
clusteredltrs_rpkm = matchrpkmtoltr(BLASTresults, cdhit);

%% matches clustered ltr retrotransposons to ltr data from ltrharvest
disp('Matching LTR data to clusters')
clusteredltrs_ltrdata = matchltrtodata(clusteredltrs_rpkm, ltrlengths, ltrname, ltrclassification);

%% function to match clustered ltr retrotransposons to promoters 
clusteredltrs_ltrdata = matchltrtopromoter(clusteredltrs_ltrdata, BLASTresults_promoters);

%% counting clusters
disp('Counting clusters')
clusteredltrs_ltrdata = clustercount(clusteredltrs_ltrdata);

%% extra graphs
extra = extragraphs(clusteredltrs_ltrdata);

%% graph heatshock data
disp('Graphing Heatshock data')
heatshock = heatshockanalysis(clusteredltrs_ltrdata);

%% graph substrate data
substratedata = substrate(clusteredltrs_ltrdata);

%% size of cluster vs number of clusters histogram - transcribed ltrs only
disp('Graphing transcribed LTRs vs Cluster size')
graph_transcribedltrs = graphtranscribedltrs(clusteredltrs_ltrdata);

%% size of cluster vs number of clusters histogram - all ltrs
disp('Graphing all LTRs vs Cluster size')
graph_all_ltrs = graphallltrs(cdhit);

%% histogram for ltr retrotransposon length
disp('Graphing histogram for LTR length')
graphretlength = graphretlength(ltrlengths);

%% graph for rpkm vs cluster
disp('Creating rpkm vs cluster size graph')
graph_rpkm_vs_cluster=rpkmclustergraph(clusteredltrs_rpkm);

%% pie chart ltr classifcations
disp('Creating pie chart for LTR classifications')
pie_chart = pieclassification(ltrclassification);



writetable(heatshock,'ltrdata_neosp1_heatshock.xlsx')
writetable(substratedata,'ltrdata_neosp1_substrate.xlsx')

saveas(graph_transcribedltrs, 'transcribedclusters_neosp1.png')
saveas(graph_all_ltrs, 'all_ltrs_clustergraph_neosp1.png')
saveas(graphretlength, 'retrolength_neosp1.png')
saveas(graph_rpkm_vs_cluster, 'rpkmclustergraph_neosp1.png')
saveas(gcf, 'piechart_neosp1.png')

%save svg
saveas(graph_transcribedltrs, 'transcribedclusters_neosp1.svg')
saveas(graph_all_ltrs, 'all_ltrs_clustergraph_neosp1.svg')
saveas(graphretlength, 'retrolength_neosp1.svg')
saveas(graph_rpkm_vs_cluster, 'rpkmclustergraph_neosp1.svg')
saveas(gcf, 'piechart_neosp1.svg')

%% function to find percent classified LTRs from TESorter
classvsunclass = classvsunclasspie(ltrclassification, ltrname);

saveas(gcf, 'classvsunclass_neosp1.svg')
saveas(gcf, 'classvsunclass_neosp1.png')