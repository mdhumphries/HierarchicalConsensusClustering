function [grpscon,ctr,k,varargout] = ConsensusSweepCommunityDetect(S,B,varargin)
 
% CONSENSUSSWEEPCOMMUNITYDETECT consensus community detection using k-means sweep, with multiple consensus matrices
%   [C,N,K] = CONSENSUSSWEEPCOMMUNITYDETECT (W,B) finds the consensus-clustering 
%   partition of the n*n weighted graph W, using k-means sweep between bounds in B = [l u]
%
%   Returns: 
%       C: column vector indicating community membership for the consensus
%       partition(s): (n rows, by c columns - one column per converged clustering)
%       N:  the number of iterations until consensus was reached. 
%       K: the number of clusters in each returned graph partition
%
%   ...= CONSENSUSSWEEP(...,gamma,OPTS) 
%           gamma: resolution parameter for the modularity matrix; set []
%           to omit (default is 1); gamma < 1 = larger communities; gamma >
%           1 = smaller communities
%           OPTS is a struct, with at least one of the fields:
%           OPTS.nreps : sets k-means to run N times for each specified group size (default is 50)
%           OPTS.dims : 'scaled' (default), 'all'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           OPTS.project : 'Laplacian' (default), 'Eigs'
%           OPTS.escape : number of consensus iterations before exiting
%           (default: 50)
%   
%   [...,bCON,CLU] = CONSENSUSSWEEP(...) where:
%       bCON : a Boolean flag indicating whether the consensus algorithm
%                converged or not
%       CLU:   a matrix of every single clustering of the similarity matrix S in the first
%              pass (i.e. before the consensus) - this is useful for further
%              post-processsing.
%
%   Notes: 
%   (0) Weighted graph should have:
%           - no self-loops (zeros on diagonal of W)
%           - only positive weights
%       Warnings for both of these are given, if detected
%
%   (1) Sequence of algorithm:
%       - Make modularity matrix Bmod from data
%       - Project Bmod
%       - use k-means on projected graph, per requested k; dimensions of
%           projection scale with k
%       - Make consensus per k [C(k)]
%       - Check if any C(k) have converged; if so, exit
%       - If not, for each k: 
%              - make consensus-specific modularity matrix Bcon(k)
%              - Project Bcon(k)
%              - Cluster Bcon(k) using k
%              - Make new consensus C(k)
%        Until converged
%
%
%   23/11/2018 : initial version, derived from CONSENSUSSWEEP
%
%   Mark Humphries 

% defaults
options.nreps = 50;     % of each k
options.dims = 'scale';   % scale embedding dimensions to number of clusters for each k-means clustering
options.project = 'Laplacian';
options.escape = 50;

gamma = 1;  % resolution parameter

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
if narargin >= 3  & ~isempty(varargin{1})
    gamma = varargin{1};
end    

if nargin >= 4
    if isstruct(options) 
        tempopts = varargin{2}; 
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

%% do initial projection and clustering
ks = B(1):B(2);
k(1) = B(2);  % number of groups found

% project the weighted graph using modularity
Bmod = W - gamma*expectedA(W);  % modularity matrix

D = ProjectEigs(Bmod,B(2));  % return all dimensions to upper bound

% get partitions using k-means
C = kmeansSweep(D,B(1),B(2),options.nreps,options.dims);
varargout{2} = C;  % store initial clusterings

%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    ctr = ctr+1  % increment consensus iteration count            
    if  ctr > options.escape
        % do escape if not converging
        warning(['Did not converge in ' num2str(options.escape) 'iterations - exiting without consensus answer'])
        grpscon = [];
        blnConverged = 0;
        return
    end

    %% now make consensus
    % loop over tested k
    allCCons = zeros(nIDs,nIDs,numel(ks));
    for iK = 1:numel(ks)
        % make consensus from corresponding part of C
        ix = 1+options.nreps*(iK-1):options.nreps*iK; % indices into C
        allCCons(:,:,iK) = makeConsensusMatrix(C(:,ix));  % make consensus

        % test for convergence, returning theta
        [blnVec(iK),thisgrps{iK},theta(iK)] = CheckConvergenceConsensus(squeeze(allCCons(:,:,iK)));   
    end

    % if *any* converged, then stop and store all converged answers
    ixConv = find(blnVec);
    if any(blnVec)
        blnConverged = 1;
        grpscon = zeros(nIDs,numel(ixConv));
        try
            for iC = 1:numel(ixConv)
                grpscon(:,iC) = thisgrps{ixConv(iC)};
            end
        catch
            keyboard
        end
        k = max(grpscon);  % return number of clusters in each partition
    end
    
    if ~blnConverged
        for iK = 1:numel(ks)
            % for each k
            % project the consensus matrix - *not* using gamma here (as
            % null model is *not* configuration model
            Bcon = allCCons(:,:,iK) - nullmodelConsensusSweep(ks(iK),options.nreps,nIDs);
            D = ProjectEigs(Bcon,ks(iK));  % return all dimensions for k

            % get partitions of that consensus matrix using k-means
            ix = 1+options.nreps*(iK-1):options.nreps*iK; % indices into C
            C(:,ix) = kmeansSweep(D,ks(iK),ks(iK),options.nreps,options.dims);  % repeat the same k
        end
    end
end

%     if blnSave
%         eval(['CCons', num2str(ctr),' = CCons;']);  % store consensus matrix
%         eval(['Groups', num2str(ctr),' = grpscon;']); 
%         save(fname,['CCons', num2str(ctr)],['Groups', num2str(ctr)],'-append');
%     end
            
% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     
function D = ProjectEigs(S,M)

% find top M eigenvectors of similarity matrix S
    [V,egs] = eig(S,'vector');
    [~,ix] = sort(egs,'descend');       % make sure eigenvalues are sorted
    D = V(:,ix(1:M));  % return all eigenvectors to upper bound




  







