function Results = clusterLowDNetwork(W,P,L,M,nreps,group_membership)

% CLUSTERLOWDNETWORK clusters low-dimensional projection of network
% R = CLUSTERLOWDNETWORK(W,P,N,M,nreps,G) clusters the (n*n) network W using
% its projection into a low-dimensonal space, to detect 
% between L and M groups (often L==M). Also needs:
%       P : expected null model used to define the low-D space
%       nreps : number of repeats of k-means at each number of dimensions
%       G : ground-truth clustering as an n-length array
% 
%  Note: will return groups = 1 and VI = 0 if either:
%       L = 1    % there is nothing to cluster: 1 group requested
%       n < M+1  % network is too small to cluster
%
% 9/7/2018 : initial version
% 19/07/2018: min network size
%
% Mark Humphries 
 
n = size(W,1);  % size of network

 if L > 1 && n > M+1      % if dimensions exist to cluster in, and enough nodes exist to make clusters...
    [Results.QmaxCluster,Qmax,Results.ConsCluster,ConsQ,~] = ...
              ConsensusCommunityDetect(W,P,L,M,nreps);
    % quality of estimation of retained communities
    if ~isempty(Results.QmaxCluster)
        [~,Results.normVIQmaxSpectra]=VIpartitions(Results.QmaxCluster,group_membership);
        Results.nGrpsQmaxSpectra = max(Results.QmaxCluster);
    else
        Results.normVIQmaxSpectra= 0;
        Results.nGrpsQmaxSpectra = nan;                
    end

    if ~isempty(Results.ConsCluster)
        [~,Results.normVIConsensusSpectra]=VIpartitions(Results.ConsCluster,group_membership);
        Results.nGrpsConsensusSpectra = max(Results.ConsCluster);
    else
        Results.normVIConsensusSpectra= 0;
        Results.nGrpsConsensusSpectra = nan;                
    end
 else
    % no modules detected: just one group!
    Results.QmaxCluster = []; Results.ConsCluster = [];
    Results.normVIQmaxSpectra=0;
    Results.normVIConsensusSpectra=0;
    Results.nGrpsQmaxSpectra = 1;
    Results.nGrpsConsensusSpectra = 1;
end
