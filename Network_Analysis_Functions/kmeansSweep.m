function allgrps = kmeansSweep(D,LowerK,UpperK,Treps,dims)

% KMEANSSWEEP sweep of kmeans clustering across multiple K
% C = KMEANSSWEEP(D,L,U,R,DIMS) performs a k-means sweep of a data matrix, given:  
% D: set of data points (nxd): n objects in d dimensions  (e.g. low-D
% projection)
% L,U: (L)ower and (U)pper number of groups to test. Tested range K = U-L (set U=L to do one value of K)
% R: repeats of each clustering at each number of groups K
% DIMS = 'all' or 'scale': 'all' uses all of D, irrespective of K; 'scale'
% uses only the i-1 set of dimensions
%
% Returns C, a m x [(K-1)xR] matrix of every k-means clustering, in order of low to high K 
%   (columns 1:K-1 are for K=L; columns K to 2*K-1 are for K=L+1 etc) 
%
% Mark Humphries 6/3/2017

if LowerK < 2 
    error('kmeansSweep:parameter','Specify at least L=2 as minimum number of groups')
end

if LowerK > UpperK
    error('kmeansSweep:parameter','Lower-bound greater than upper-bound')
end

[n,d] = size(D);
if any(strfind(dims,'scale')) && d < UpperK-1
    error('kmeansSweep:parameter','Not enough embedding dimensions to scale to upper bound')
end

K = 1+ UpperK - LowerK;
ixNow = 1;
allgrps = zeros(n,(K-1)*Treps);
for ngrps = LowerK:UpperK
    switch dims
        case 'all'
            thisVector = D;
        case 'scale'
            thisVector = D(:,ngrps-1); % subset of dimensions
        otherwise
            error('Unknown options for DIMS')
    end
    
    for rep = 1:Treps

        cpos = kmeansplus(thisVector, ngrps); % initialise centers (chooses starting position at random)

        try            
            allgrps(:,ixNow) = kmeans(thisVector,ngrps,'Distance','sqeuclidean','Start',cpos);

        catch
            % if kmeans throws a wobbly, set to "no groups"...
            warning('kmeans wobbly')
            keyboard
        end
        ixNow = ixNow + 1;
    end % end k-means repeats loop
end % end groups loop