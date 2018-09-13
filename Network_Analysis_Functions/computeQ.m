function Q = computeQ(C,B,varargin)

% COMPUTEQ modularity score of a multiple-partition 
% Q = COMPUTEQ(C,B) computes the modularity score Q, given:
%   C: a column vector of group membership IDs (integer values, one per
%   node)
%   B: the modularity matrix
% Q=COMPUTEQ(...,M) normalises Q by 2*M  
%   M: the number of unique edges or sum of unique weights (i.e. that undirected networks have had each link counted once only)
%
%   Q = sum_ij(A_ij-P_ij) I(n_i,n_j) [where I(n_i,n_j) = 1 if nodes i,j
%   are in the same group, and 0 otherwise]
%
%   Reference:  Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
% Mark Humphries 2/3/2017

n = size(B,1);
grpIDs = unique(C);
ngrps = numel(grpIDs);

if nargin > 2
    m = varargin{1};
else
    m = 1;
end

% construct S matrix of group membership: each column is a group
% See: Newman (2006) Eq 31            
S = zeros(n,ngrps);

for loop = 1:ngrps
    S(:,loop) = (C == grpIDs(loop));
end

% compute modularity
Q = trace(S' * B * S) / (2*m);  % note: assumes here that m is the sum of unique weights (i.e. that undirected networks have had each link counted once only)
