clc
clear all
close all

%list all genomes
genomes = {'neosp1', 'neolan1', 'caecom1','pirfi3','uh31','gfma','anasp1'};

%outer for loop is for iterating through the different databases
for i = 1:length(genomes)
    %define current database
    current_db = genomes{i};
    
    %define j (iteration for the queries)
    j=1;
    %while loop
    while j <= length(genomes)
        %define current query
        current_query = genomes{j};
        %skip when current db = curent query, because it will be 100%
        if strcmp(current_query, current_db)
            j=j+1;
            %will give an error if this isn't here, because when
            %anasp1=anasp1 it will increment to where j>length of genomes
            if j > length(genomes)
                break;
            end
        end
        %redefine current query
        current_query = genomes{j};
        %import blast results
        blast.(current_db).(current_query) = importblast([current_db,'_', current_query, '.txt']);
        %import ltr harvest sequence names
        ltrharvest.(current_db).(current_query) = importdata_ltrharvest(['ltrsequences_', current_query, '.txt']);
        %take unique LTRs in query column from blast file
        uniqueValues.(current_db).(current_query) = unique(blast.(current_db).(current_query)(:, 1));
        %calculate percent
        percent.(current_db).(current_query) = 100*size(uniqueValues.(current_db).(current_query),1)/size(ltrharvest.(current_db).(current_query),1);
        %increment j
        j=j+1;
    end
    %display results
    disp(current_db)
    disp(percent.(current_db))
end


%% simpler code, it does the same thing but easier to follow
% %import blast results
% blast = importblast('neolan1_neosp1.txt');
% 
% %import ltr harvest results
% ltrharvest = importdata_ltrharvest('ltrsequences_neosp1.txt');
% 
% %take unique values of from first column (i.e. the query)
% uniqueValues = unique(blast(:, 1));
% 
% 
% find percent
% percent = 100*size(uniqueValues,1)/size(ltrharvest,1)