function [D,egs,varargout] = ProjectLaplacian(W,varargin)

% PROJECTLAPLACIAN low-dimensional projection using Laplacian
% D = PROJECTLAPLACIAN(W) for the n*n symmetric similarity matrix W, 
% creates a low-dimensional projection of the data using 
% the random-walk Laplacian matrix.
% 
% Finds the number of eigenvectors to retain using the eigengap heuristic
%
% [D,E] = PROJECTLAPLACIAN(W) will also return the n-length vector of
% eigenvalues
%
% ... = PROJECTLAPLACIAN(...,M) sets limits on the number of dimensions, if M is:
%           a vector [L,U], will compute the eigengap between eigenvalues L and U
%           a scalar, will return exactly M leading eigenvectors, ignoring the eigengap
%
% Returns:
%   D: the n*M array of n-length vectors corresponding to the M embedding
%   dimensions (i.e. the M eigenvectors corresponding to the M *smallest* eigenvalues)
%   
% Notes:
%   (1) Computes the random-walk Laplacian by default
%
% References:  
% von Luxburg, U. (2007) A tutorial on spectral clustering. Statistics and Computing, 17, 395-416
%
% 7/3/2017: initial version
% 13/11/2018: added eigengap heuristic 
%
% Mark Humphries
N = size(W,1);                     % number of objects 

Lbnd = 1; Ubnd = N;   % default: check all eigenvalues

if nargin > 1
    if numel(varargin{1}) == 1 % passing a scalar
        Lbnd = varargin{1}; Ubnd = Lbnd;
    else
        Lbnd = varargin{1}(1); Ubnd = varargin{1}(2);
    end
end

D = diag(sum(W));             % degree matrix: count of unique edges on the diagonal
L = D-W;
Lrw = D' * L;                   % Laplacian random walk
% keyboard


%% FIX: double eigenvector computation
[V,egs] = eig(Lrw,'vector');         % get eigenvalues of Laplacian (RW)
[egs,ix] = sort(egs,'ascend');    % sort into *ascending* order: we want the *smallest* eigenvalues here [Luxburg pg 3, last line]
V = V(:,ix);

if Ubnd==Lbnd
    k = Lbnd;
else
    % eigengap heuristic, between specified bounds
    k = (Lbnd-1) + find(max(diff(egs(Lbnd:Ubnd))));    % generally fails miserably
end

keyboard

[~,D] = eigs(Lrw,k);                      % embedding dimensions

varargout{1} = egs;     % for debugging

                
