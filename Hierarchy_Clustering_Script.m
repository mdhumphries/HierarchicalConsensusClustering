%% template script for running hierarchical consensus clustering
% Mark Humphries 13/9/2018

clearvars; close all

addpath('Network_Analysis_Functions/')
addpath('Helper_Functions/')

% network to analyse
fname = 'Allen_Gene_Leaf'; 


% clustering parameters
clusterpars.nreps = 100;
clusterpars.explore = 'explore';  % allow consensus to use more groups than specified by spectral rejection


%% load data-file
load(['Networks/' fname]); 

% deal with individual differences in data formats:
% A : full weight matrix
% nodelabels : a string array, one string per row
if strfind(fname,'cosyne')
    A = adjMatrix;
    m = cellfun('length',cosyneData.authorHash);
    nodelabels = [];
    for i = 1:numel(cosyneData.authorHash)
        l = numel(cosyneData.authorHash{i});
        nodelabels = [nodelabels; cosyneData.authorHash{i} blanks(max(m) - l)];
    end
    nodelabels = nodelabels;
% elseif strfind(fname,'CosyneYear')
%     A = adjMatrix;
%     nodelabels = nodelabel;
elseif exist('Problem')
    A = full(Problem.A);
    % Generate node labels for later visualisation to work
    if isfield(Problem,'aux')
        nodelabels = Problem.aux.nodename;
    else
        nodelabels = string(1:size(A,1))';
    end
end

% make undirected if necessary
A = (A + A') / 2; % make undirected

% clean-up A, get largest component, and store as basis for all further analysis
% all indices are with reference to Data.A
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);
Data.nodelabels = nodelabels(Data.ixRetain,:);   % update the node labels

%% make first level of hierarchy

% TO SOLVE HERE:
% (1) how to get limits on cluster numbers? 
%     (a) Could use spectral rejection of course; but then this would mean using all that code and explaining it all in the report... 
%    (b) or just a fixed number...
% (2) use a different consensus 

% simple eigen-gap solution

P = expectedA(Data.A);
B = A - P;
Edata = sort(eig(B),'descend');


% (1) spectral clustering

% (2) make consensus


% [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
%           ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn,clusterpars.nreps,[],clusterpars.explore);
       
       

%% create hierarchy

%% view some results