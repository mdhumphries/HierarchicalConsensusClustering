function [grpscon,ctr,k,varargout] = ConsensusLouvain(W,varargin)
 
% CONSENSUSLOUVAIN consensus partition using Louvain algorithm
%   [C,N] = CONSENSUSLOUVAIN(S) finds the consensus-clustering 
%   partition of the similarity matrix S, using the Louvain algorithm for
%   community detection
%
%   Returns: 
%       C: column vector indicating group membership for the partition with maximum
%           modularity
%       N:  the number of iterations until consensus was reached. 
%
%   [...] = CONSENSUSLOUVAIN(...,R)
%        R: set the number of repeats of the Louvain algorithm (default: 20)      
%   
%   [...,bCON,CLU] = CONSENSUSLOUVAIN(...) where:
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
%
%   References: 
%   (1) Blondel et al 2008
%
%   19/11/2018 : initial version, derived from older code
%
%   Mark Humphries 

% defaults
nreps = 20;     % of Louvain

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

%% options
if nargin > 1
    nreps = varargin{1};
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

[allC,allQ,allCn,allIters] = LouvainCommunityUDnondeterm(W,nreps,1); % first level of hierarchy


keyboard

% create C


%% check for exit if error on first step
varargout{1} = blnConverged;
varargout{2} = C;
varargout{3} = Wk;

if isempty(C) 
    % then no groups detected; return empty
    grpscon = zeros(nIDs,1);
    ctr = 0; 
    blnConverged = 0;  % no answer
    return
end
%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    % make consensus
    CCons = makeConsensusMatrix(C);
    
    %% check convergence
    [blnConverged,grpscon] = CheckConvergenceConsensus(CCons);
    
    
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
        
        % do k-means on projection of consensus matrix
        [allC,allQ,allCn,allIters] = LouvainCommunityUDnondeterm(W,nreps,1);
        

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
    
     




  







