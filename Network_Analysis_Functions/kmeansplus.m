function m = kmeansplus(X,k)
% KMEANPLUS centroid initialisation for k-means
% C = KMEANSPLUS(X,K) given mxd data matrix X (m objects in d dimensions to be clustered), 
% and K clusters, find the most-dispersed set of centroids C (a K*d matrix
% of centroid positions).
%
% From: Arthur, D. & Vassilvitskii, S. (2007) k-means++: the advantages of careful seeding. 
%   SODA '07: Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, Society for Industrial and Applied Mathematics, 1027-1035
%
% Mark Humphries 2/3/2017

[n,d] = size(X); 
m = zeros(k,d);
v = inf(n,1); % current minimum distance of data-point to closest centre
m(1,:) = X(ceil(n*rand),:); % put 1st co-ordinates into 

for i = 2:k
    % find distance from all data-points to closest center
    D = sqrt(sum((X - repmat(m(i-1,:),n,1)).^2,2));  % distance of all data-points from last picked center
    v = min(v,D);  % update minimum distance to closest centre
    p = cumsum(v/sum(v));  % ECDF
    m(i,:) = X(find(rand < p,1),:);  % pick next point proportional to p
end