function f = matchltrtopromoter(clusteredltrs_ltrdata, BLASTresults_promoters)

promoters = strings(height(clusteredltrs_ltrdata),1);
clusteredltrs_ltrdata = addvars(clusteredltrs_ltrdata,promoters);
BLASTresults_promoters = sortrows(BLASTresults_promoters,5,'descend');

%matches classification to cluster
for i=1:size(clusteredltrs_ltrdata,1)
    idx = find(strcmp([BLASTresults_promoters{:,2}], clusteredltrs_ltrdata{i,1}));
    if (isempty(idx)) %this is for the rows that say Cluster #
        clusteredltrs_ltrdata{i,'promoters'} = "0";
    else
        clusteredltrs_ltrdata{i,'promoters'}= BLASTresults_promoters{idx(1),1};
    end
end

f= clusteredltrs_ltrdata;
