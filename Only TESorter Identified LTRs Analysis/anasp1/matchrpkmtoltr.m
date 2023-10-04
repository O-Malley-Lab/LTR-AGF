%function to match rpkm values to each individual ltr and cluster
function f = matchrpkmtoltr(BLASTresults, cdhit)

rpkm = zeros(height(BLASTresults),1);
BLASTresults = addvars(BLASTresults,rpkm);


%removes the beginning 'locus' part, so rpkm is just a number
for i=1:size(BLASTresults,1)
    if endsWith(BLASTresults{i,1},'_PRE')
        BLASTresults{i,'rpkm'} = extractBetween(BLASTresults{i,1}, "rpkm","_PRE");
    
    else
        BLASTresults{i,'rpkm'} = extractAfter(BLASTresults{i,1}, "rpkm");
    end

end

%find repeating ltr names, save the max rpkm
for i=1:size(BLASTresults,1)
    maxval = [];
    maxidx =[];
    idx = find(strcmp([BLASTresults{:,'LTRName'}], BLASTresults{i,'LTRName'}));
    [maxval, maxidx] = max(BLASTresults{idx,'rpkm'});
    BLASTresults{idx,'rpkm'} = maxval;
    BLASTresults{idx,'Transcriptome'} = BLASTresults{idx(maxidx),1};
    BLASTresults{idx,'AlignmentLength'} = BLASTresults{idx(maxidx),4};
    BLASTresults{idx,'AlignmentPercentIdentity'} = BLASTresults{idx(maxidx),3};
end

BLASTresults = unique(BLASTresults, 'stable');

%% add new columns to cdhit table
%add rpkm column
rpkm = zeros(height(cdhit),1);
cdhit = addvars(cdhit,rpkm);

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
    cdhit{idx,'rpkm'} = BLASTresults{i,'rpkm'};
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



%this while loop takes the mode of an rpkm value, if clusters don't have
%the same rpkm
data = [];
transcriptome = strings;
d=1; %cdhit increment
j=1; %data table increment
count = 0; 
while d < size(cdhit,1)
    
    if startsWith(cdhit{d,1},"Cluster")
        d=d+1; %increase increment at start of loop
        
        if d >= size(cdhit,1) %end statement to close loop, size is continuously decreasing so this is needed
            break;
        end

        while startsWith(cdhit{d,1},"scaffold") %once iteration hits an ltr
            count = count + 1; %increase cluster count
            data(j,1) = cdhit{d,'rpkm'}; %save rpkm values of the same cluster to a new data array
            transcriptome(j,1) = cdhit{d,'Transcript'}; %save the transcript names of one cluster to another array
            d=d+1;
            j=j+1;
        end
        
        
        modeval = mode(data); %find mode of cluster rpkm array
        modeidx = find(abs(data - modeval) < eps, 1, 'first'); %find index of the first time that mode appears
        
        for i = d - count:d-1 %the first ltr starts at d - cluster size. the cdhit increment (d) ends at one after the last ltr, so we need to do d-1 to get the last increment.
             cdhit{i,'rpkm'} = modeval; %update rpkm with mode
             if isinteger(transcriptome(modeidx))
                cdhit{i,'Transcript'} = transcriptome(modeidx); %update transcript
             end
             if ~isinteger(transcriptome(modeidx))
                break;
             end

        end
        
        count = 0;
        data =[];
        j=1;
    
    end
    

end

f = cdhit;

end