function Sc = makeConsensusMatrix(Grps)

% MAKECONSENSUSMATRIX make a consensus matrix from a set of clusterings
% C = MAKECONSENSUSMATRIX(G) given n*c matrix of c clusterings of n objects
% computes the consensus matrix C. Entry C_ij is the proportion of times
% object (i,j) were in the same cluster
%
% Mark Humphries 2/3/2017

[nIDs,nreps] = size(Grps);
Sc = zeros(nIDs);

% pair-wise similarity matrix
for nr = 1:nIDs
    for nC = nr:nIDs
        if nr ~= nC
            Gi = Grps(nr,:); Gj = Grps(nC,:);
            Sc(nr,nC) = sum(Gi == Gj) / nreps;
            Sc(nC,nr) = Sc(nr,nC);
        end
    end
end
     