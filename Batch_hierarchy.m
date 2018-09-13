%%% script run batch of all hierarchies of recordings starting from consensus
%%% clustering spike-train data
%%% Mark Humphries 27/11/2013

load(['DataSet_ConsensusClustering']); % original data-set: 


%% do hierarchy:

nRecordings = numel(Gcon_dataset);
Gh_All = cell(nRecordings,1);
Ghsizes = cell(nRecordings,1);

for iR = 1:nRecordings
    iR
    % get top-level of spike-train groups
    Gh_All{iR}{1} = Gcon_dataset{iR}.grps;   % stores groups as each node assignment...
    Ghsizes{iR}{1} = Gcon_dataset{iR}.grpsizes;
    ntrains = length(Gh_All{iR}{1});

    blnStop = 0; lvl = 0;
    while ~blnStop
        lvl = lvl + 1
        % build super-graph...
        nnodes = numel(Ghsizes{iR}{lvl}); % number of nodes is the prior number of groups
        A = zeros(nnodes); Cmap = cell(nnodes,1);
        for iN = 1:nnodes
            srcs = find(Gh_All{iR}{lvl}(:,2) == iN);
            Cmap{iN} = srcs;  % store indices for mapping
            for iC = iN+1:nnodes
                if iN ~= iC
                    tgts = find(Gh_All{iR}{lvl}(:,2) == iC);

                    tmp = Sxy_dataset{iR}(srcs,tgts);
                    % generate link: mean similarity between groups
                    A(iN,iC) =  sum(sum(tmp)) / (numel(srcs)*numel(tgts)); % mean similarity between the two groups       
                    A(iC,iN) = A(iN,iC);
                end
            end
        end
        
        % cluster super-graph....
        [C,Qmax,Ccon,Qc,N,Q] = allevsplitConTransitive(A,{'sqEuclidean'},200); 
        
        if isempty(C)
            % if no groups, then stop
            blnStop = 1;
        else
            % map indices back to original nodes...
            thisC = unique(Ccon); origG = zeros(ntrains,2); origG(:,1) = Gh_All{iR}{1}(:,1);
            for iC = 1:numel(thisC)
                merged = find(Ccon == thisC(iC));
                Ghsizes{iR}{lvl+1}(iC) = numel(merged); 
                for i = 1:numel(merged)
                    origG(Cmap{merged(i)},2) = thisC(iC);  % assign all nodes in old community to this community index
                end
                Gh_All{iR}{lvl+1} = origG;
            end
        end

        if numel(thisC) == 2
            % if reached two groups, then also stop
            blnStop = 1;
        end

    end
end

save Dataset_Consensus_Hierarchy Gh_All Ghsizes
