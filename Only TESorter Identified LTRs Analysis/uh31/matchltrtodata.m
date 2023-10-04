%this function matches ltr data from ltr harvest and classifies each LTR
function f = matchltrtodata(clusteredltrs_rpkm, ltrlengths,ltrname,ltrclassification)

%creates new columns in table for data
retlength = zeros(height(clusteredltrs_rpkm),1);
clusteredltrs_rpkm = addvars(clusteredltrs_rpkm,retlength);

classification = strings(height(clusteredltrs_rpkm),1);
clusteredltrs_rpkm = addvars(clusteredltrs_rpkm,classification);


%matches retrotransposon ltr length to cluster
for i=1:size(clusteredltrs_rpkm,1)
    idx = find(strcmp([ltrname{:,1}], clusteredltrs_rpkm{i,1}));
    if (isempty(idx)) %this is for the rows that say Cluster #
        clusteredltrs_rpkm{i,'retlength'} = 0;
    else
        clusteredltrs_rpkm{i,'retlength'}= ltrlengths{idx,1};
    end
end

%matches classification to cluster
for i=1:size(clusteredltrs_rpkm,1)
    idx = find(strcmp([ltrclassification{:,1}], clusteredltrs_rpkm{i,1}));
    if (isempty(idx)) %this is for the rows that say Cluster #
        clusteredltrs_rpkm{i,'classification'} = "0";
    else
        clusteredltrs_rpkm{i,'classification'}= ltrclassification{idx,2};
    end
end

clusteredltrs_rpkm((strcmp(clusteredltrs_rpkm.classification,'0') & ~startsWith(clusteredltrs_rpkm.Clusterno_LTRName,'Cluster')), :) = [];

%removes the cluster headers for the clusters that weren't transcribed
k=0;
l=0;
while k < size(clusteredltrs_rpkm,1)
    l=l+1;
    k=l+1;
    if startsWith(clusteredltrs_rpkm{l,'Clusterno_LTRName'},"Cluster") && startsWith(clusteredltrs_rpkm{k,'Clusterno_LTRName'},"Cluster")
       clusteredltrs_rpkm(l,:) = [];
       l=l-1;
    end
end


f=clusteredltrs_rpkm;
