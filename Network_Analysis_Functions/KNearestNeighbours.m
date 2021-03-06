function Wk = KNearestNeighbours(W,K,varargin)

% KNEARESTNEIGHBOURS create a K-nearest-neighbours graph
% WK = KNEARESTNEIGHBOURS(W,K) turns the NxN similarity matrix W into a
% K-nearest neighbour graph WK. 
%   K can be:
%       integer > 1: take exactly that many neighbours 
%       scalar in the range (0,1] : take round(KN) neighbours
%
% In WK, entry (i,j) and (j,i) = W(i,j) if connected, 0 otherwise
%
% 14/11/2018: initial version
% Mark Humphries

N = size(W,1);

if K <= 1
    K = round(N*K);
end

% find neighbours
[~,I] = sort(W,'descend');

Wk = zeros(N);
for iN = 1:N
    Wk(I(1:K,iN),iN) = W(I(1:K,iN),iN);
    Wk(iN,I(1:K,iN)) = W(I(1:K,iN),iN); % make undirected
end





