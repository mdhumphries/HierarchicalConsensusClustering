function [E,V] = nullmodelConsensusSweep(K,T,N)

% NULLMODELCONSENSUSSWEEP k-means sweep null model for consensus matrices
% [E,V] = NULLMODELCONSENSUSSWEEP(K,T,N) constructs the expected null model for a
% consensus matrix of random clusterings from a k-means sweep, given:
%       K: a M-length vector giving a list of the number of clusters tested 
%       T: a M-length vector, giving the number of clusterings for each entry of K.
%       N: the number of objects clustered
%
% Returns:
%   E: the N*N matrix of expected proportions of clusterings
%   V: the NxN matrix of variance in the proportion of clusterings
%
% Mark Humphries 7/3/2017

if numel(T) ~= numel(K)
    error('nullmodelConsensusSweep:parameters','T and K need to be the same length');
end

% sum over binomial distribitions B(T_k,p_k)
for iK = 1:numel(K)
    P = 1/K(iK);                 % for each draw, the uniformly random probability of being in the same cluster
    E(iK) = T(iK)*P; % sum(T);      % how many draws?
    V(iK) = P*T(iK)*(1-P);          % variance of draws
    % keyboard
end
E = zeros(N) + sum(E)/sum(T);   % expectation
V = zeros(N) + sum(V)/sum(T);

E(eye(N)==1) = 0;           % no self-connections
V(eye(N)==1) = 0;           % no self-connections