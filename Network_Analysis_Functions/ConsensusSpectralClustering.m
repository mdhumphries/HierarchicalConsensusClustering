function [grpscon,ctr,k,varargout] = ConsensusSpectralClustering(W,Knn,varargin)
 
% CONSENSUSSPECTRALCLUSTERING consensus partition using spectral clustering
%   [C,N,K] = CONSENSUSSPECTRALCLUSTERING(S,Knn) finds the consensus-clustering 
%   partition of the similarity matrix S, using spectral clustering. The
%   similarity matrix is turned into a K-nearest-neighbours graph using
%   Knn.
%
%   Returns: 
%       C: column vector indicating group membership for the partition with maximum
%           modularity
%       N:  the number of iterations until consensus was reached. 
%       K: the number of clusters detected at each iteration
%
%   ...= CONSENSUSSPECTRALCLUSTERING(...,N,DIMS,'fixed') 
%           N : sets k-means to run  N times for each specified group size (default is 50)
%           DIMS : 'all' (default), 'scaled'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           'fixed': set this option to force consensus algorithm to always
%           use initially identified number of clusters
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
blnExplore = 1;     % allow consensus algorithm to determine its own number of clusters

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

if nargin >= 7 && strcmp(varargin{3},'fixed') 
    blnExplore = 0;     % force consensus algorithm to 
end    

% % set up saving of each iteration
% blnSave = 0; % internal flag for setting saving of data
% if blnSave  
%     fname = ['Consensus_Iterations_' datestr(now,30)];
%     save(fname,'dims','nreps','blnExplore');  % save initial data to allow -append to work below
% end

%% internal parameters

nIDs = size(W,1);     % number of nodes of the similarity matrix
blnConverged = 0;       % stopping flag
ctr = 1;                % iterations of consensus

%% cluster similarity matrix

% if matrix is dense, make sparse using k-nearest-neighbours
Wk = KNearestNeighbours(W,Knn);

% project that graph; finds number of groups using eigengap heuristic
[D,egs] = ProjectLaplacian(Wk);

% get partitions using k-means
k(ctr) = size(D,2);  % number of groups to find
C = kmeansSweep(D,k(ctr),k(ctr),nreps,dims);


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
    
    keyboard
    
%     if blnSave
%         eval(['CCons', num2str(ctr),' = CCons;']);  % store consensus matrix
%         eval(['Groups', num2str(ctr),' = grpscon;']); 
%         save(fname,['CCons', num2str(ctr)],['Groups', num2str(ctr)],'-append');
%     end
            
    %% sort out what to do next
    if ~blnConverged 
        ctr = ctr+1;  % increment consensus iteration count            
        if  ctr > 50
            % do escape if not converging
            warning('Did not converge in 50 iterations - exiting without consensus answer')
            grpscon = [];
            blnConverged = 0;
            return
        end
        
        if blnExplore
            disp('I am exploring')
            [D,egs] = ProjectLaplacian(CCons);  % use fully-connected here?
            k(ctr) = size(D,2);  % number of groups to find defined by new projection
        else
            % just use the same number of groups as the original requested
            % set: find consensus in that space.
            % so return that number of eigenvectors
            disp('I am not exploring')
            k(ctr) = k(1);  % number of groups to find is fixed to the initial set
            [D,egs] = ProjectLaplacian(CCons,k(ctr));  % use fully-connected here?

        end
        % do k-means on projection of consensus matrix
        C = kmeansSweep(D,k(ctr),k(ctr),nreps,dims);


        if isempty(C)
            warning('Consensus matrix projection is empty - exiting without consensus answer')
            grpscon = [];
            blnConverged = 0;               
            return
        end
    end
end

% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     




  







