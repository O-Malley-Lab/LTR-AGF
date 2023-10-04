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

j=1;
d=1;
e=1;
f=1;
h=1;
for i=1:height(clusteredltrs_rpkm)
    if startsWith(clusteredltrs_rpkm{i,'classification'},"Gypsy")
        classifictiontable.gypsy(j,1) = clusteredltrs_rpkm(i,'rpkm');
        j=j+1;
    end
    if startsWith(clusteredltrs_rpkm{i,'classification'},"Copia")
        classifictiontable.copia(d,1) = clusteredltrs_rpkm(i,'rpkm');
        d=d+1;
    end
    if startsWith(clusteredltrs_rpkm{i,'classification'},"mixture")
        classifictiontable.mixture(e,1) = clusteredltrs_rpkm(i,'rpkm');
        e=e+1;
    end
    if startsWith(clusteredltrs_rpkm{i,'classification'},"unknown")
        classifictiontable.unknown(f,1) = clusteredltrs_rpkm(i,'rpkm');
        f=f+1;
    end
    if startsWith(clusteredltrs_rpkm{i,'classification'}, '0') && startsWith(clusteredltrs_rpkm{i,'Clusterno_LTRName'},'scaffold')
        classifictiontable.nontesorter(h,1) = clusteredltrs_rpkm(i,'rpkm');
        h=h+1;
    end
end

f=clusteredltrs_rpkm;
