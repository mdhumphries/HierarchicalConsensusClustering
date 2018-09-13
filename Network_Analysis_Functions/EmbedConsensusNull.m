function [D,B,M,varargout] = EmbedConsensusNull(C,Model,varargin)

% EMBEDCONSENSUSNULL low-dimensional projection and N groups using null model for consensus matrix
% [D,B,M] = EMBEDCONSENSUSNULL(C,MODEL) for the n*n consensus matrix C, with values
% in [0,1], produces a low-dimensional projection and estimates the maximum
% number of groups within it.
%           MODEL is a string that sets the null model:
%           'sweep': if the consensus matrix was obtained from a k-means
%           sweep. Specify additional arguments (...,K,T): K is a list of
%           cluster sizes tried in the sweep; T is the number of
%           clusterings tried at each size [both M-length vectors].
%
%           'expect': to derive the null model as the expectiation of the
%           consensus matrix
%   
% Returns:
%   D: the n*d array of n-length vectors corresponding to the d embedding
%   dimensions (i.e. the top d eigenvectors)
%   B: the n*n modularity matrix
%   M: the estimated number of groups in C (= d+1)
%
% Notes:
%   (1) 'sweep' uses NULLMODELCONSENSUSSWEEP
%   (2) 'expect' uses NULLMODELCONSENSUSEXPECTATION
%
% Mark Humphries 7/3/2017

N = size(C,1);                     % number of objects 
egmin = 1e-5;

if any(strfind(Model,'sweep')) && nargin < 4
    error('EmbedConsensusNull:Parameter','Missing arguments for sweep model')
end

if nargin > 2
    K = varargin{1};
    T = varargin{2};
end

switch Model
    case 'sweep'
        [E,~] = nullmodelConsensusSweep(K,T,N);
%        [E,Var] = nullmodelConsensusSweep(K,T,N);
%         CI = CIfromSEM(sqrt(Var),zeros(N)+T*numel(K),0.95);
%         E95 = E + CI;
    case 'expect'
        E = nullmodelConsensusExpectation(C);
    otherwise
        error('EmbedConsensusNull:Parameter','Unknown null model option')
end

% make modularity matrix using null model
B = C - E;

% get eigenvectors
[V,egs] = eig(B,'vector');         % get eigenvalues and vectors of C
[egs,ix] = sort(egs,'descend');    % sort into descending order
V = V(:,ix);                       % ditto the eigenvectors 

ixE = find(egs > egmin);             % eigenvalues that exceed model
% ixE = knee_pt(egs,1:N);            % eigenvalues beyond knee point 
D = V(:,ixE);                      % embedding dimensions
M = numel(ixE)+1;                  % number of groups


varargout{1} = egs;

                
