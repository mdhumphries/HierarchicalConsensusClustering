function [allC,allQ,allCn,allIters] = LouvainCommunityUDnondeterm(W,R,varargin)

% LOUVAINCOMMUNITYUDnondeterm Louvain community detection algorithm, undirected
% [C,Q,N,I] = LOUVAINCOMMUNITYUDnondeterm(W,R) for undirected, weighted network W; runs
% the Louvain community detection algorithm from Blondel et al (2008).
% Constructs a grouping at each hierarchical level.  Runs the algorithm R
% times, choosing a different node order each time. 
%
% Returns: 
% C, a cell array of length R: each cell contains a cell array of the community structure at each level; 
% Q, a cell array of length R: each cell contains an array of the modularity scores at each level;  
% N, a cell array of length R: each cell contains an cell array of group sizes at each level; 
% I, a cell array of length R: each cell contains an array of the number of iterations at each level
%
%... = LOUVAINCOMMUNITYUDnondeterm(...,1) finds the first level of the
%hierarchy only
%
% NOTES: this differs from the LouvainCommunityUD function in one small but
% crucial respect: whereas that function checks each node in the same order
% (from 1 to N) to do the iterative clustering, this version chooses a
% different node order each time. Thus, the algorithm is run many times,
% and the set of answers for each run is returned.
%
% TO DO: return struct...
%
% References:
% Blondel, V. D.; Guillaume, J.-L.; Lambiotte, R. & Lefebvre, E. Fast
% unfolding of communities in large networks J. Stat. Mech, 2008, P10008
% 
% Mark Humphries 27/09/2012

bln1 = 0; % do all levels of hierarachy

if nargin >= 3 & varargin{1} == 1;   % do first level only
    bln1 = 1;
end

threshold = 0;  % default threshold: any increase in Q

Worig = W;  % stupidly memory hungry; but keeps things simple
allC = cell(R,1);
allQ = cell(R,1);
allIters = cell(R,1);
allCn = cell(R,1);

%% main algorithm loop
for iR = 1:R
    [Norig c] = size(W); 
    blnAllGain = 1; % continuation flag
    nLevels = 0; % number of hierarchy levels
    W = Worig;  % on loop>=2, W now between-modules, so restore [stupidly memory hungry; but keeps things simple]
    allC{iR} = {}; allQ{iR} = []; allIters{iR} = []; allCn{iR} = {};
    while blnAllGain == 1
        nLevels = nLevels + 1;
        blnGain = 1;
        nIter = 0;
        % get properties of current W matrix
        [N,c] = size(W);    % number of nodes
        k = sum(W);         % degree: number of edges
        m = sum(k); 

        C = 1:N;    % initial community assignment

        % stage 1: cluster at this level of hierarchy
        while blnGain == 1
            blnGain = 0;    % set to no gain in modularity
            nIter = nIter + 1;
            list =  randperm(N); % current node order
            for iN = 1:N     
                i = list(iN);   % index of current node
                Ngbrs = find(W(:,i));   % idxes of neighbouring nodes [could find once, but needs lots of storage]
                deltaQ = zeros(numel(Ngbrs,1));

                % compute node i's modularity for current assignment
                Gcurr = find(C == C(i));    % all nodes in current group
                nullmodel = k(i) .* k(Gcurr) ./ m; % expected number of links in current assignment
                q_i_curr = sum(W(i,Gcurr) - nullmodel);      
                % compute gain for each new assignment
                for j = 1:numel(Ngbrs)
                    newC = C(Ngbrs(j)); % community of node j
                    if newC ~= C(i)  % only check if neighbour is in different group
                        Gnew = [i find(C == newC)];   % indices of all nodes in new community (including itself...)
                        % gain of moving i into C(j)
                        nullmodel = k(i) .* k(Gnew) ./ m; % expected number of links in new assignment
                        q_i_new = sum(W(i,Gnew) - nullmodel);   

                        deltaQ(j) = q_i_new - q_i_curr;     % gain in modularity due to move
                    end
                end

                % keep new assignment if greater than threshold
                ix = find(deltaQ > threshold & deltaQ == max(deltaQ));   
                if ~isempty(ix)
                    % if more than 1 with exact same gain then random assignment
                    slct = randsample(numel(ix),1);

                    % or just take first node! Get same results as other code
                    slct = 1;
                    newC = C(Ngbrs(ix(slct)));  % the community of max deltaQ
                    C(i) = newC;                % assign node i
                    blnGain = 1;    % gain in modularity
                end
            end
            % keyboard
        end

        allIters{iR} = [allIters{iR}; nIter];   

        % test for stopping
        if (nIter == 1 & blnGain == 0)    
            % only one iteration without gain: time to end!
            blnAllGain = 0; 
            break
        end

        % store results: unpack if above level 1 of hierarchy
        if nLevels == 1
            C = resortC(C);
            % keyboard
            allC{iR} = [allC{iR} C'];    % store in overall 
            allQ{iR} = [allQ{iR}; modularity(W,C)];    % compute overall modularity
            if bln1 == 1
                blnAllGain = 0; break
            end
        else
            % re-map from new community assignment to per-node assignment
            % before re-sorting!
            Corig = zeros(Norig,1);
            thisC = unique(C);
            for iC = 1:numel(thisC)
                merged = find(C == thisC(iC));  % communities with same new community assignment
                for i = 1:numel(merged)
                    Corig(Cmap{merged(i)}) = thisC(iC);  % assign all nodes in old community to this community index
                end
            end
            % keyboard
            C = resortC(Corig); % renumber 
            allC{iR} = [allC{iR} C];    % store in overall 
            allQ{iR} = [allQ{iR}; modularity(Worig,C)];    % compute overall modularity
            % keyboard
        end


        % stage 2: combine
        % make new W... and make corresponding smaller C...

        nC = numel(unique(C));
        Wnew = zeros(nC); Cmap = cell(nC,1);
        for iC = 1:nC
            srcs = find(C == iC); % nodes in this community
            Cmap{iC} = srcs;    % keep map of old node assignments
            for tC = iC:nC    % includes self-self connections when iC = tC
                tgts = find(C == tC); % nodes in target community
                Wnew(iC,tC) = sum(sum(Worig(srcs,tgts))); % assumes undirected 
                Wnew(tC,iC) = Wnew(iC,tC);
            end
        end
        W = Wnew; % remap and start again

    end % while still gaining in modularity...
    
    % count community sizes
    if nLevels > 1 
        nL = nLevels - 1; % last level did not get sorted 
    else
        nL = 1; 
    end
    for iL = 1:nL 
        nC = numel(unique(allC{iR}{iL})); Cs = zeros(nC,1);
        for iC = 1:nC
            Cs(iC) = sum(allC{iR}{iL} == iC); % number of nodes in this community   
        end
        allCn{iR} = [allCn{iR} Cs]; % store number of nodes
    end
    
end % end of loop over node-ordering

end  %%% end of main function

function newC = resortC(C)
    % at end of iterations through same hierarchy level, community number
    % assignment will be non-contiguous
    ixs = sort(unique(C));
    NG = numel(ixs);
    newC = C;
    for i = 1:NG
        newC(C == ixs(i)) = i;  
    end
end

function Q = modularity(W,C)
    nG = numel(unique(C));
    N = numel(C);
    % null model
    inks = sum(W);
    outks = sum(W');
    m = sum(inks);  % number of edges
    P = (outks' * inks) ./ m;   % null model
    
    % community matrix
    S = zeros(N,nG);

    for loop = 1:nG
        S(:,loop) = (C == loop);
    end
    
    Q = trace(S' * (W-P) * S) ./ m;   % scaled by 1/2m...
    
%     % for debugging...
%     rQ = 0;
%     k = sum(W);
%     for i = 1:N
%         % compute node i's modularity for current assignment
%         Gcurr = find(C == C(i));    % all nodes in current group
%         nullmodel = k(i) .* k(Gcurr) ./ m; % expected number of links in current assignment
%         q_i = sum(W(i,Gcurr) - nullmodel);     
%         rQ = rQ + q_i;
%     end
%     rQ = rQ ./ m;   % running Q
    % keyboard
end