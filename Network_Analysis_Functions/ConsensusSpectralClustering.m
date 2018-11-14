function [grpscon,ctr,k,varargout] = ConsensusSpectralClustering(W,varargin)
 
% CONSENSUSSPECTRALCLUSTERING consensus partition using spectral clustering
%   [C,N] = CONSENSUSSPECTRALCLUSTERING(S) finds the consensus-clustering 
%   partition of the similarity matrix S, using spectral clustering.
%
%   expected null model P, and minimum L and maximum M number of groups to detect.
%   
%   Returns: 
%       C: column vector indicating group membership for the partition with maximum
%           modularity
%       N:  the number of iterations until consensus was reached. 
%       K: the number of clusters detected at each iteration
%
%   ...= CONSENSUSSPECTRALCLUSTERING(...,L,U,N,DIMS,'explore') 
%           N : sets k-means to run  N times for each specified group size (default is 50)
%           DIMS : 'all' (default), 'scaled'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           'explore': if L==U, then set this option to let the consensus
%               algorithm explore greater numbers of groups
%   
%   [...,bCON,CLU,Wk] = CONSENSUSSPECTRALCLUSTERING(...) where:
%       bCON : a Boolean flag indicating whether the consensus algorithm
%                converged or not
%       CLU:   a matrix of every single clustering of the similarity matrix S in the first
%              pass (i.e. before the consensus) - this is useful for further
%              post-processsing.
%       Wk:     the similarity matrix resulting from the
%       k-nearest-neighbours step
%
%   Notes: 
%   (0)For a similarity matrix, ensure: 
%       (i) no self-loops - diagonal of W is all zeros; 
%       (ii) it's a similiarity matrix, not a correlation matrix: no negative values
%   Warnings for both of these will be given
%
%   (1) Spectral clustering version:
%       - k-nearest-neighbours of S
%       - Laplacian random walk of resulting matrix (Sk)
%       - pick number of clusters K using eigengap in eigenvalues of
%       Laplacian; project onto K top eigenvectors
%       - use k-means on projected objects
%
%   References: 
%   (1) Luxburg
%
%   14/11/2018 : initial version, derived from older code
%
%   Mark Humphries 

% defaults
nreps = 50;     % of each distance metric
dims = 'all';   % use all embedding dimensions for each k-means clustering
blnExplore = 1;

%% check if the passed matrix is a graph: catch common errors when passing a similarity matrix

% (1) no self-loops allowed
if ~all(diag(W)==0) 
    warning('Results likely unreliable: similarity matrix has self-loops. Set diagonal to zero')
end

% (2) no negative links
x = sum(sum(W < 0));
if x > 0
    warning('Results likely unreliable: similarity matrix has negative values')
end

%% set up options
if nargin >= 5 && ~isempty(varargin{1}) 
    nreps = varargin{1}; 
end

if nargin >= 6 && ~isempty(varargin{2}) 
    dims = varargin{2}; 
end    

if nargin >= 7 && strcmp(varargin{3},'explore') 
    blnExplore = 1;
end    

% % set up saving of each iteration
% blnSave = 0; % internal flag for setting saving of data
% if blnSave  
%     fname = ['Consensus_Iterations_' datestr(now,30)];
%     save(fname,'dims','nreps');  % save initial data to allow -append to work below
% end

%% internal parameters

nIDs = size(W,1);     % number of nodes of the similarity matrix
blnConverged = 0;       % stopping flag
ctr = 1;                % iterations of consensus

%% cluster similarity matrix

% if matrix is dense, make sparse using k-nearest-neighbours
Wk = KNearestNeighbours(Data.A,clusterpars.K);

% project that graph; finds number of groups using eigengap heuristic
[D,egs] = ProjectLaplacian(Wk);

% get partitions using k-means
k(ctr) = size(D,2);  % number of groups to find
C = kmeansSweep(D,k(ctr),k(ctr),clusterpars.nreps,'all');


%% check for exit if error on first step
varargout{1} = blnConverged;
varargout{2} = C;
varargout{3} = Wk;

if isempty(C) 
    % then no groups detected; return empty
    grpscon = zeros(nIDs,1);
    ctr = 0; 
    k = 0;
    blnConverged = 0;  % no answer
    return
end
%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    % make consensus
    CCons = makeConsensusMatrix(C);
    
    %% check convergence
    [blnConverged,grpscon] = CheckConvergenceConsensus(CCons);
    
    %% sort out what to do next
    if blnConverged
        Qcon = computeQ(grpscon,B,m);  % compute Q, and exit
    else
        ctr = ctr+1;  % increment consensus iteration counter
        if ctr > 50
            % do escape if not converging
            warning('Did not converge in 50 iterations - exiting without consensus answer')
            grpscon = [];
            blnConverged = 0;
            return
        else
            

            if blnExplore
                disp('I am exploring')
                [D,egs] = ProjectLaplacian(CCons);  % use fully-connected here?

                % get partitions using k-means
                k(ctr) = size(D,2);  % number of groups to find
                C = kmeansSweep(D,k(ctr),k(ctr),clusterpars.nreps,'all');


            else
                % just use the same number of groups as the original requested
                % set: find consensus in that space.
                disp('I am not exploring: to be finished')
                [D,egs] = ProjectLaplacian(CCons);  % use fully-connected here?

                % get partitions using k-means
                k(ctr) = k(1);  % number of groups to find
                C = kmeansSweep(D,k(ctr),k(ctr),clusterpars.nreps,'all');
            end
                                

            if isempty(C)
                warning('Consensus matrix projection is empty - exiting without consensus answer')
                grpscon = [];
                Qcon = 0;
                blnConverged = 0;               
                return
            end
        end
    end
end

% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     




  







