
function [W,ixRetain,Comps,Comp_sizes] = prep_A(W)
% PREP_A get final weighted adjacency matrix
% [W*,I,C,S] = PREP_A(W) given the undiredted weighted adjacency matrix W,
% finds the largest component, and returns that as W*, the matrix to be
% used in all further analyses. Also returns:
%   I: the index of retained nodes in W*
%   C: the complete index list of components (1,...,C components)
%   S: the sizes of each component
%
% NOTE:
% Uses GET_COMPONENTS from the Brain Connectivity Toolbox
%
% Mark Humphries 9/3/2017

A = double(W>0);
[Comps,Comp_sizes] = get_components(A);

% Keep only the largest component
ixC = find(Comp_sizes == max(Comp_sizes));
ixRetain = find(Comps == ixC);
W = W(ixRetain,ixRetain);

% Remove diagonal elements
W((eye(size(W,1)))==1) = 0;



