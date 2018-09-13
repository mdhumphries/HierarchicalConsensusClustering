function [D,B,M,varargout] = EmbedConsensus(C)

% EMBEDCONSENSUS low-dimensional projection and N groups in consensus matrix
% [D,B,M] = EMBEDCONSENSUS(C) for the n*n consensus matrix C, with values
% in [0,1], produces a low-dimensional projection and estimates the maximum
% number of groups within it.
%   
% Returns:
%   D: the n*d array of n-length vectors corresponding to the d embedding
%   dimensions (i.e. the top d eigenvectors)
%   B: the n*n modularity matrix
%   M: the estimated number of groups in C
%
% Notes:
%   (1) Assumes the null model that the consensus matrix is a random matrix
%   if the clustered objects are unrelated. Thus we use random matrix
%   theory to find embedding dimensions that exceed this null model
%
% References:  
%   MacMahon, M. & Garlaschelli, D. (2013) Community detection for correlation
%   matrices. arxiv 1311.1924
%   Edelman, A. & Wang, Y. (2013) Random Matrix Theory and Its Innovative Applications Advances in 
%   Applied Mathematics, Modeling, and Computational Science, Springer Nature, 91-116
%
% Mark Humphries 2/3/2017

N = size(C,1);                     % number of objects 
C = (C - mean(C(:))) ./ var(C(:)); % normalise to mean 0 and var=1
H = (C+C')/2;                    % make Hermitian

[V,egs] = eig(H,'vector');         % get eigenvalues and vectors of C
[egs,ix] = sort(egs,'descend');    % sort into descending order
V = V(:,ix);                       % ditto the eigenvectors 

egsNorm = egs ./ sqrt(N/2);        % normalise eigenvalues to upper limit of Wishaw

% NB: later version will use Tracy-Widom distribution here for max eigenvalue
ixE = find(egsNorm > 2);             % eigenvalues that exceed upper limit of Wishaw

D = V(:,ixE);                      % embedding dimensions
M = numel(ixE)+1;                  % number of groups

% construct modularity matrix: this is "signal" reconstructed from retained
% eigenvectors
B = zeros(N);
for i=1:numel(ixE)
    B = B + egs(ixE(i)) * V(:,ixE(i))*V(:,ixE(i))';   
end

varargout{1} = egsNorm;

                
