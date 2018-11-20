function [grpscon,ctr,k,varargout] = ConsensusSweep(S,B,varargin)
 
% CONSENSUSSWEEP consensus partition using k-means sweep
%   [C,N,K] = CONSENSUSSWEEP(S,B) finds the consensus-clustering 
%   partition of the similarity matrix S, using k-means sweep between bounds in B = [l u]
%
%   Returns: 
%       C: column vector indicating group membership for the partition with maximum
%           modularity
%       N:  the number of iterations until consensus was reached. 
%       K: the number of clusters detected at each iteration (useful only
%       for 'each' option)
%
%   ...= CONSENSUSSWEEP(...,OPTS) 
%           OPTS.combine : 'each' (default), 'all'; 'each' produces a
%           consensus matrix per tested k, and chooses k as the least
%           ambiguous clustering; 'all' produces a consensus matrix over
%           all k 
%           OPTS.nreps : sets k-means to run N times for each specified group size (default is 50)
%           OPTS.dims : 'scaled' (default), 'all'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           OPTS.project : 'Laplacian' (default), 'Eigs'
%           OPTS.escape : number of consensus iterations before exiting
%           (default: 50)
%   
%   [...,bCON,CLU,Wk] = CONSENSUSSWEEP(...) where:
%       bCON : a Boolean flag indicating whether the consensus algorithm
%                converged or not
%       CLU:   a matrix of every single clustering of the similarity matrix S in the first
%              pass (i.e. before the consensus) - this is useful for further
%              post-processsing.
%
%   Notes: 
%   (0)For a similarity matrix, ensure: 
%       (i) no self-loops - diagonal of W is all zeros; 
%       (ii) it's a similiarity matrix, not a correlation matrix: no negative values
%   Warnings for both of these will be given
%
%   (1) Sequence of algorithm:
%       - Project data: Laplacian random walk of similarity matrix
%       (default), or eigenvectors of S
%       - use k-means on projected objects, per requested k; dimensions of
%           projection scale with k
%       - Make consensus: either per k; or over all k
%       - If made consensus per k, choose which to use going forward by 
%           - checking if any have converged
%           - if not, finding which is least ambiguous (greatest distance
%           between modes)
%       - Repeat from projection    
%
%   References: 
%   (1) Luxburg
%
%   20/11/2018 : initial version, derived from CONSENSUSSPECTRALCLUSTERING
%
%   Mark Humphries 

% defaults
options.nreps = 20;     % of each k
options.dims = 'scale';   % scale embedding dimensions to number of clusters for each k-means clustering
options.combine = 'each';   % consensus per k
options.project = 'Laplacian';
options.escape = 50;
%% check if the passed matrix is a graph: catch common errors when passing a similarity matrix

% (1) no self-loops allowed
if ~all(diag(S)==0) 
    warning('Results likely unreliable: similarity matrix has self-loops. Set diagonal to zero')
end

% (2) no negative links
x = sum(sum(S < 0));
if x > 0
    warning('Results likely unreliable: similarity matrix has negative values')
end

%% set up options
if nargin >= 3
    if isstruct(options) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            options = setfield(options,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end


% % set up saving of each iteration
% blnSave = 0; % internal flag for setting saving of data
% if blnSave  
%     fname = ['Consensus_Iterations_' datestr(now,30)];
%     save(fname,'dims','nreps','blnExplore');  % save initial data to allow -append to work below
% end

%% internal parameters

nIDs = size(S,1);     % number of nodes of the similarity matrix
blnConverged = 0;       % stopping flag
ctr = 0;                % iterations of consensus


%% loop: consensus matrix and its clustering until converged
ks = B(1):B(2);
k(1) = B(2);  % number of groups found
L = B(1); U = B(2);  % limits on k-means sweep
CCons = S;  % on first loop, use similarity matrix as the consensus matrix

while ~blnConverged
    ctr = ctr+1;  % increment consensus iteration count            
    if  ctr > options.escape
        % do escape if not converging
        warning(['Did not converge in ' num2str(options.escape) 'iterations - exiting without consensus answer'])
        grpscon = [];
        blnConverged = 0;
        return
    end

    % project the graph
    switch options.project
        case 'Laplacian'
            [D,~] = ProjectLaplacian(CCons,k(ctr));  % return all dimensions to upper bound
        case 'Eigs'
            [D,~] = ProjectEigs(CCons,k(ctr));  % return all dimensions to upper bound
        otherwise
            error('Unknown option for projecting data')
    end

    % get partitions using k-means
    C = kmeansSweep(D,L,U,options.nreps,options.dims);
    if ctr == 1  % store initial clusterings
        varargout{2} = C;
    end


    %% now make consensus
    switch options.combine
        case 'all'
            % make consensus over all sweep
            CCons = makeConsensusMatrix(C);

            % check convergence
            [blnConverged,grpscon] = CheckConvergenceConsensus(CCons);
            
            if ~blnConverged
                k(ctr+1) = B(2);  % max K for next time
                % to consider:  how to extend upper limit of k sweep if needed?
            else
                k(ctr+1) = max(grpscon);  % number of groups in converged answer
            end
            
        case 'each'
            % loop over tested k
            for iK = 1:numel(ks)
                % make consensus from corresponding part of C
                ix = 1+options.nreps*(iK-1):options.nreps*iK; % indices into C
                thisCCons = makeConsensusMatrix(C(:,ix));  % make consensus
                
                % test for convergence, returning theta
                [blnVec(iK),thisgrps{iK},theta(iK)] = CheckConvergenceConsensus(thisCCons);

                % find modes, and get distance between them
                try
                    allCs = thisCCons(triu(ones(size(thisCCons)),1));
                catch
                    keybaord
                end
                m1 = median(allCs <= theta);
                m2 = median(allCs > theta);
                d(iK) = m2 - m1;  % distance between modes of consensus entries
            end
            
            keyboard
            % if *any* converged, then stop
            ixConv = find(blnVec);
            if any(blnVec)
                blnConverged = 1;
                % if more than one converged
                if numel(ixConv) > 1
                   keyboard
                   % should be the same surely?
                end
                grpscon = thisgrps{ixConv};
                k(ctr+1) = max(grpscon);
            else
                % else choose consensus matrix to take forward...
                ixMax = find(d == max(d));
                
                % reconstruct one with best k
                ix = 1+options.nreps*(ixMax-1):options.nreps*ixMax;
                CCons = makeConsensusMatrix(C(:,ix));  % make consensus
            end
        otherwise 
            error('Unknown option for combining clusterings into consensus')
    end
    
    keyboard
    
%     if blnSave
%         eval(['CCons', num2str(ctr),' = CCons;']);  % store consensus matrix
%         eval(['Groups', num2str(ctr),' = grpscon;']); 
%         save(fname,['CCons', num2str(ctr)],['Groups', num2str(ctr)],'-append');
%     end
            
end

% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     
function D = ProjectEigs(S,M)

% find top M eigenvectors of similarity matrix S
    [V,egs] = eig(S,'vector');
    [~,ix] = sort(egs,'descend');       % make sure eigenvalues are sorted
    D = V(:,ix(1:M));  % return all eigenvectors to upper bound




  






