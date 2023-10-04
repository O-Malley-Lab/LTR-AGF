%function to match rpkm values to each individual ltr and cluster
function f = matchrpkmtoltr(BLASTresults, cdhit)

rpkm = zeros(height(BLASTresults),1);
BLASTresults = addvars(BLASTresults,rpkm);


BLASTresults = unique(BLASTresults, 'stable');

BLASTresults = sortrows(BLASTresults , 'AlignmentLength' );

%% add new columns to cdhit table

%add transcript column
Transcript = strings(height(cdhit),1);
cdhit = addvars(cdhit,Transcript);
cdhit = standardizeMissing(cdhit,0);

%add percent identity column
AlignmentPercentIdentity = strings(height(cdhit),1);
cdhit = addvars(cdhit,AlignmentPercentIdentity);

%add alignment length column
AlignmentLength = strings(height(cdhit),1);
cdhit = addvars(cdhit,AlignmentLength);


%matches rpkm to ltr
for i=1:size(BLASTresults,1)
    idx = strcmp([cdhit{:,1}], BLASTresults{i,'LTRName'});
    cdhit{idx,'Transcript'} = BLASTresults{i,'Transcriptome'};
    cdhit{idx,'AlignmentPercentIdentity'} = BLASTresults{i,'AlignmentPercentIdentity'};
    cdhit{idx,'AlignmentLength'} = BLASTresults{i,'AlignmentLength'};
end

%enters 0 for rpkm values for cluster rows - needs to be done so the cluster number rows aren't removed
for i=1:size(cdhit,1)
    if startsWith(cdhit{i,1},"Cluster")
        cdhit{i,2} = "0";
    end
end



f = cdhit;

end