function [D,varargout] = ProjectLaplacian(W,M)

% PROJECTLAPLACIAN low-dimensional projection using Laplacian
% D = PROJECTLAPLACIAN(W,M) for the n*n symmetric matrix W, creates a low-dimensional
% projection using the M top eigenvectors of the random-walk Laplacian
% matrix.
%   
% Returns:
%   D: the n*d array of n-length vectors corresponding to the M embedding
%   dimensions (i.e. the M eigenvectors corresponding to the M *smallest* eigenvalues)
%   
% Notes:
%   (1) Computes the random-walk Laplacian by default
%
% References:  
% von Luxburg, U. (2007) A tutorial on spectral clustering Statistics and Computing, 17, 395-416
%
% Mark Humphries 7/3/2017

N = size(W,1);                     % number of objects 

D = diag(sum(W));             % degree matrix: count of unique edges on the diagonal
Lrw = eye(N) - D' * W;             % Laplacian random walk 

% keyboard

[V,egs] = eig(Lrw,'vector');         % get eigenvalues and vectors of C
[egs,ix] = sort(egs,'ascend');    % sort into *ascending* order: we want the *smallest* eigenvalues here [Luxburg pg 3, last line]
V = V(:,ix);                       % ditto the eigenvectors 

D = V(:,1:M);                      % embedding dimensions

varargout{1} = egs;     % for debugging

                
