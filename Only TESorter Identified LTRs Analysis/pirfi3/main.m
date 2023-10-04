close all
clear all
clc
set(groot, 'DefaultAxesFontSize', 14);

%% importing data

%import LTR name for matching from LTRHarvest
ltrname = importltrname('ltrsequences_pirfi3.txt');

%import LTR data from LTRHarvest
ltrlengths = importltrlength('ltrdata_pirfi3.txt');

%Import BLAST transcriptome matches
BLASTresults = importdatablast('blastresultspirfi3.txt');
rpkmresults = importrpkm('blastresultspirfi3.txt');

%Import clustered LTRs
cdhit = importdatacdhit('cdhitresultsaSaL70_pirfi3.txt');

%Import LTR classification from TEsorter
ltrclassification = importtesorter('tesorterresults_pirfi3.txt');

%Import heatshock data (from Candice's data)
%heatshockdata = importdataheatshock('G1_heatshock_raw_expectedcountsWithB3.txt');

%% function to match clustered ltr retrotransposons to rpkm values
disp('Matching rpkm to LTRs in clusters')
clusteredltrs_rpkm = matchrpkmtoltr(BLASTresults, cdhit);

%% matches clustered ltr retrotransposons to ltr data from ltrharvest
disp('Matching LTR data to clusters')
clusteredltrs_ltrdata = matchltrtodata(clusteredltrs_rpkm, ltrlengths, ltrname, ltrclassification);
%% counting clusters
disp('Counting clusters')
clusteredltrs_ltrdata = clustercount(clusteredltrs_ltrdata);
%% graph heatshock data
%disp('Graphing Heatshock data')
%heatshock = heatshock(clusteredltrs_ltrdata, heatshockdata);

%% size of cluster vs number of clusters histogram - transcribed ltrs only
disp('Graphing transcribed LTRs vs Cluster size')
graph_transcribedltrs = graphtranscribedltrs(clusteredltrs_ltrdata);

%% size of cluster vs number of clusters histogram - all ltrs
disp('Graphing all LTRs vs Cluster size')
graph_all_ltrs = graphallltrs(cdhit,ltrclassification);

%% histogram for ltr retrotransposon length
%disp('Graphing histogram for LTR length')
%graphretlength = graphretlength(ltrlengths);

%% graph for rpkm vs cluster
disp('Creating rpkm vs cluster size graph')
graph_rpkm_vs_cluster=rpkmclustergraph(clusteredltrs_rpkm);

%% pie chart ltr classifcations
disp('Creating pie chart for LTR classifications')
pie_chart = pieclassification(ltrclassification);

writetable(clusteredltrs_ltrdata,'ltrdata_pirfi3.xlsx')

saveas(graph_transcribedltrs, 'transcribedclusters_pirfi3.png')
saveas(graph_all_ltrs, 'all_ltrs_clustergraph_pirfi3.png')
%saveas(graphretlength, 'retrolength_pirfi3.png')
saveas(graph_rpkm_vs_cluster, 'rpkmclustergraph_pirfi3.png')
saveas(gcf, 'piechart_pirfi3.png')


saveas(graph_transcribedltrs, 'transcribedclusters_pirfi3.svg')
saveas(graph_all_ltrs, 'all_ltrs_clustergraph_pirfi3.svg')
%saveas(graphretlength, 'retrolength_pirfi3.svg')
saveas(graph_rpkm_vs_cluster, 'rpkmclustergraph_pirfi3.svg')
saveas(gcf, 'piechart_pirfi3.png')

%% function to find percent classified LTRs from TESorter
classvsunclass = classvsunclasspie(ltrclassification, ltrname);

saveas(gcf, 'classvsunclass_pirfi3.svg')
saveas(gcf, 'classvsunclass_pirfi3.png')