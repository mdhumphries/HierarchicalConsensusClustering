%% template script for running hierarchical consensus clustering
% Mark Humphries 13/9/2018

clearvars; close all

addpath('Network_Analysis_Functions/')
addpath('Helper_Functions/')

% network to analyse
fname = 'Allen_Gene_Leaf'; 


% clustering parameters
clusterpars.nreps = 100;        % of k-means

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
%     (b) original approach of just +ve eigs (but too many)
%     (c) bisection method
%     (d) or just a fixed number...
% (2) use a different clustering methods to get initial set for
% constructing consensus matrix

% [grpscon,ctr,k] = ConsensusSpectralClustering(Data.A,clusterpars.K);

% Or do something else entirely? (e.g. Louvain)       
%[grpscon,ctr] = ConsensusLouvain(Data.A);
clusterpars.project = 'Laplacian';
[grpscon,ctr,k] = ConsensusSweep(Data.A,[2,20],clusterpars);

clusterpars.project = 'Eigs';
[grpscon,ctr,k] = ConsensusSweep(Data.A,[2,20],clusterpars);



%% create hierarchy

%% view some results