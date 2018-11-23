%% template script for running hierarchical consensus clustering on data
% Mark Humphries 13/9/2018

clearvars; close all

addpath('Network_Analysis_Functions/')
addpath('Helper_Functions/')

% network to analyse
fname = 'Allen_Gene_Leaf'; 

% clustering parameters
clusterpars.nreps = 100;        % of k-means

%% load data-file
load(['Datasets/' fname]); 

%% cluster
clusterpars.project = 'Laplacian';
[grpscon,ctr,k] = ConsensusSweep(A,[2,20],clusterpars);

clusterpars.project = 'Eigs';
[grpscon,ctr,k] = ConsensusSweep(A,[2,20],clusterpars);


%% view some results

