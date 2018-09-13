function [grps,Qmax,grpscon,Qcon,ctr,varargout] = ConsensusCommunityDetect(W,P,L,M,varargin)
 
% CONSENSUSCOMMUNITYDETECT partition signal network using eigenvectors of signal modularity matrix (with consensus)
%   [C,Qmax,Ccon,Qc,N,Q] = CONSENSUSCOMMUNITYDETECT(W,P,L,M) splits the
%   vertices of the nxn weighted, undirected signal network W into multiple groups, given
%   expected null model P, and minimum L and maximum M number of groups to detect.
%   
%   Returns: 
%       C: column vector indicating group membership for the partition with maximum
%           modularity
%       Qmax: the corresponding modularity score. 
%       Ccon: column vector indicating group membership for the consensus
%           clustering partition [see Notes]; 
%       Qcon: the corresponding is the modularity score for the consensus partition.
%       N:  the number of iterations until consensus was reached. 
%   
%   ...= CONSENSUSCOMMUNITYDETECT(...,N,DIMS,'explore') 
%           N : sets k-means to run  N times for each specified group size (default is 50)
%           DIMS : 'all' (default), 'scaled'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           'explore': if M==L, then set this option to let the consensus
%               algorithm explore greater numbers of groups
%   
%   [...,bCON,CLU] = CONSENSUSCOMMUNITYDETECT(...) where:
%   bCON : a Boolean flag indicating whether the consensus algorithm
%   converged or not
%   CLU is an optional output argument, 
%   returns every single clustering of the adjacency matrix A in the first
%   pass (i.e. before the consensus) - this is useful for further
%   post-processsing.
%
%   Notes: 
%   (0) When analysing time-series, W can be the similarity matrix. However, when
%   starting from a similarity matrix, ensure: 
%       (i) no self-loops - diagonal of W is all zeros; 
%       (ii) it's a similiarity matrix, not a correlation matrix: no negative values
%   Warnings for both of these will be given
%
%   (1) This is a one-step multiple partition method. The algorithm implemented takes the C such eigenvectors, and
%   uses k-means clustering on those eigenvectors to cluster the nodes into k = C+1 groups. 
%   A value for Q is computed for each k-means clustering (using the defined distance metrics).
%
%   Q = 1/2m * sum_ij(A_ij-P_ij) I(n_i,n_j) [where I(n_i,n_j) = 1 if nodes i,j
%   are in the same group, and 0 otherwise]
%   
%   (2) This is repeated for each C in L:M; set L=M for using just one set
%   of groups
%
%   (3) Consensus: this attempts to extract a stable set of groups that are robust to repeats
%   of the clustering process. All clusterings with Q>0 across all k-means variants and numbers of groups are
%   pooled. A consensus matrix is computed (Lancichinetti & Fortunato 2012): entry p_ij gives the
%   proportion of clusterings that placed nodes i and j in the same group.
%   A loop of cluster-then-consensus-matrix is repeated until the consensus matrix defines a unique grouping.
%
%   (5) For the community-detection algorithm, kmeans centres are initialised using the kmeans++ algorithm (Arthur & Vassilvitskii, 2007)
%
%   References: 
%   (1) Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
%   (2) Reichardt & Bornhaldt (2006) "Statistical mechanics of community detection".
%   Phys Rev E. 74, 016110
%   
%   (3) Lancichinetti, A. & Fortunato, S. (2012) Consensus clustering in complex networks.
%   Scientific Reports, 2, 336
%   
%   (4) Arthur, D. & Vassilvitskii, S. (2007) k-means++: the advantages of careful seeding. 
%   SODA '07: Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, Society for Industrial and Applied Mathematics, 1027-1035
%
%   2/3/2017 : initial version
%   16/07/2018: made consensus exploration optional
%   17/07/2018: fixed bugs in forced consensus
%   Mark Humphries 

% defaults
nreps = 50;     % of each distance metric
dims = 'all';   % use all embedding dimensions for each k-means clustering
blnExplore = 0;

%% check if the passed matrix is a graph: catch common errors when passing a similarity matrix

% (1) no self-loops allowed
if ~all(diag(W)==0) 
    warning('Results likely unreliable: adjacency matrix has self-loops. Set diagonal to zero if no self-loops are needed.')
end

% (2) no negative links
x = sum(sum(W < 0));
if x > 0
    warning('Results likely unreliable: adjacency matrix has negative values')
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

nIDs = size(W,1);     % number of nodes of the weight matrix
m = sum(sum(W))/2;    % number of unique links (or total unique weights)

blnConverged = 0;       % stopping flag
ctr = 1;                % iterations of consensus

%% cluster signal network
B = W - P;          % initial modularity matrix, given data matrix W and specified null model P
[V,egs] = eig(B,'vector');
[~,ix] = sort(egs,'descend');    % sort into descending order
V = V(:,ix);                       % ditto the eigenvectors 

C = kmeansSweep(V(:,1:M-1),L,M,nreps,dims);  % find groups in embedding dimensions: sweep from L to M
% C = kmeansSweep(V(:,1:M),L,M,nreps,dims);  % find groups in embedding dimensions: sweep from L to M

for iQ = 1:size(C,2)
    Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering
end

%% check for exit or store Qmax results
varargout{1} = blnConverged;
if isempty(C) || all(Q(:) <= 0)
    % then no groups detected; return empty
    grps = zeros(nIDs,1);
    Qmax = 0; 
    grpscon = zeros(nIDs,1);
    Qcon = 0;
    ctr = 0; 
    varargout{2} = C;
    blnConverged = 0;  % no answer
    return
else
    Qmax = max(Q);
    ix = find(Qmax == Q);
    grps = C(:,ix(1));  
    varargout{2} = C;
end
%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    % make consensus
    Allowed = (Q > 0);       % only take consensus using Q-positive clusterings...
    CCons = makeConsensusMatrix(C(:,Allowed));
    
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
            Qcon = 0;
            blnConverged = 0;
            return
        else
            
            % if a single size of groups was requested, and not exploring
            % just use the same number of groups as the original requested
            % set: find consensus in that space.
            if L == M && ~blnExplore
                disp('I am not exploring')
                T = sum(reshape(Allowed,nreps,1+M-L));  % count how many at each K were retained
                [D,~,Mcons] = EmbedConsensusNull(CCons,'sweep',L:M,T);  % option 'expected' available as well as 'sweep'
                % do k-means sweep using D, restricted to original M groups
                % (thus M-1 dimensions)
                if Mcons >= M
                    C = kmeansSweep(D(:,1:M-1),L,M,nreps,dims);  % find groups in embedding dimensions
                elseif isempty(D)
                    C = [];
                else
                    % use all of D if less than M returned
                    C = kmeansSweep(D,L,M,nreps,dims);  % find groups in embedding dimensions 
                end
            end
            
%             % find upper limit of groups - replicate this code when using
%             null model for consensus 
%             [D,~,Mcons] = EmbedConsensusWishart(CCons);
%              % do k-means sweep using found M
%             C = kmeansSweep(D,L,Mcons,nreps,dims);  % find groups in embedding dimensions
           
            % keyboard
            if L~=M || blnExplore
                disp('I am exploring')
                T = sum(reshape(Allowed,nreps,1+M-L));  % count how many at each K were retained
                [D,~,Mcons] = EmbedConsensusNull(CCons,'sweep',L:M,T);  % option 'expected' available as well as 'sweep'
                M = Mcons;
                 % do k-means sweep using found M
                C = kmeansSweep(D,L,M,nreps,dims);  % find groups in embedding dimensions
            end
            
%             % do Laplacian on consensus matrix, using original M
%             D = ProjectLaplacian(CCons,M);
%             keyboard
%             C = kmeansSweep(D,M,M,nreps,dims);  % find groups in embedding dimensions

            if isempty(C)
                warning('Consensus matrix projection is empty - exiting without consensus answer')
                grpscon = [];
                Qcon = 0;
                blnConverged = 0;               
                return
            else
             
                % compute Q
                Q = zeros(size(C,2),1);
                for iQ = 1:size(C,2)
                    Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering using original modularity matrix
                end
            end
        end
    end
end

% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     




  







