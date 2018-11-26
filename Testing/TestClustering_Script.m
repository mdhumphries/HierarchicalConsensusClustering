%% template script for running hierarchical consensus clustering on data
% Mark Humphries 13/9/2018

clearvars; close all

addpath('../Network_Analysis_Functions/')
addpath('../Helper_Functions/')

% get test data
W.N = 4;
W.size = 100;
W.within1.m = 0.5;
W.within1.s = 0.1;
W.within2.m = 0.25;
W.within2.s = 0.1;
W.between.m = 0.1;
W.between.s = 0.1;

S = MakeTestData(W);
figure
imagesc(S); colormap(gray)

% clustering parameters
clusterpars.nreps = 100;        % of k-means

%% cluster
clusterpars.project = 'Laplacian';
[grpscon,ctr,k,bln,Clu] = ConsensusSweep(S,[2,10],clusterpars);



%% test parameter options
clusterpars.project = 'Eigs';
[grpscon,ctr,k] = ConsensusSweep(S,[2,10],clusterpars);

options.escape = 3;
[grpscon,ctr,k] = ConsensusSweep(S,[2,20],clusterpars);


%% view some results

% clustermap