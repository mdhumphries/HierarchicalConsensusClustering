%% template script for running hierarchical consensus clustering on data
% Mark Humphries 13/9/2018

clearvars; close all

addpath('Network_Analysis_Functions/')
addpath('Helper_Functions/')

% clustering parameters
clusterpars.nreps = 100;        % of k-means

clster = 'RGC'; % {'Allen','RGC','white wine')

% load data-file
switch clster
    case 'Allen'
        load('Datasets/Allen_Gene_Leaf'); 
         
    case 'RGC'
        load('C:\Users\lpzmdh\Dropbox\Analyses\Networks\datasets\Retinal Ganglion Cell gene expression\RGC_R_Matrix.mat');
        A = RGC_R_Matrix;
    case 'red wine'
        load('C:\Users\lpzmdh\Dropbox\Analyses\Networks\datasets\Wine quality\Wine_red');
        A = Data.Pearson;
        
    case 'white wine'
        load('C:\Users\lpzmdh\Dropbox\Analyses\Networks\datasets\Wine quality\Wine_white');
        A = Data.Pearson;
        
end

%% cluster
% clusterpars.project = 'Laplacian';
% [grpscon,ctr,k] = ConsensusSweep(A,[2,20],clusterpars);

clusterpars.project = 'Eigs';
[grpscon,ctr,k] = ConsensusSweep(A,[2,50],clusterpars);


%% view some results

switch clster
    case {'red wine','white wine'}
        Tgts = Quality - min(Quality)+1; % numbering 1...n (assumes 
        if diff(unique(Tgts)) > 1 keyboard; end
        VI = zeros(numel(k),1); VIn = VI;
        for iC = 1:numel(k)
            [VI(iC),VIn(iC)] = VIpartitions(grpscon(:,iC),Tgts);
        end
end