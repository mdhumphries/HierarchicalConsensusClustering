function [blnTrans,grpscon,varargout] = CheckConvergenceConsensus(C)

% CHECKCONSENSUSCONVERGENCE determine if consensus clustering has converged
% [B,G] =  CHECKCONSENSUSCONVERGENCE(C) checks if the n*n consensus matrix C
% contains a single clustering of the n objects.
%
% Returns: B = {0,1} a Boolean flag for whether the matrix converged (1) or
% not (0); G, the consensus clustering that resulted: a n-element column vector 
%
% [...,T] =  CHECKCONSENSUSCONVERGENCE also returns the detected threshold
% T between the two modes
%
% Change log:
% 02/03/2017: initial version
% 21/05/2018: Otsu test, and option to return threshold
%
% Mark Humphries

nIDs = size(C,1);  % number of objects
% get upper-triangular entries of consensus matrix
idx = find(triu(ones(size(C)),1)); % upper triangular above diagonal;
allCs = C(idx);

binwidth = 0.01;  % histogram for Otsu

% convergence of weights?
% if all(allCs == 1)
%     theta = -inf;   % converged to all 1s, so no need to check for bimodality
% else
%     try
%         % initialised centroids dynamically according to spread of data
%         s =  prctile(allCs,[5,95]);
%         idx = kmeans(allCs,2,'Start',s'); % use k-means to separate into mode groups
%     catch
%         keyboard
%     end
%     
%     % work out which is lower and which is higher distribution, and set
%     % dividing line of modes (theta)
%     m1 = mean(allCs(idx==1));  mn1 =  min(allCs(idx==1));
%     m2 = mean(allCs(idx==2));  mn2 =  min(allCs(idx==2));        
%     if m1 > m2  
%         theta = mn1;
%     elseif m2 > m1 
%         theta = mn2;
%     end
% end

% use Otsu's method to divide bimodal histogram in two
% key choice: histogram bin widths
if all(allCs == 1)
    theta = -inf;   % converged to all 1s, so no need to check for bimodality
else
    [hst,edges] = histcounts(allCs,'BinWidth',binwidth); % ,'BinLimits',[0 1]);
    bintheta = otsu1D(hst);
    theta = edges(bintheta);
end

A_upper = C;
A_upper(C < theta) = 0;  % remove all of lower mode links from matrix

% find groups deterministically...
grpscon = [0 0]; grpctr = 1; blnTrans = 1;
for iN = 1:nIDs
    if ~any(iN == grpscon(:,1))  % if not already sorted
        thisgrp = [iN; find(A_upper(iN,:) > 0)']; % all nodes connected to this one are in the same group
        % if any of this group are already in the storage matrix,
        % then this is not transitive...
        for iT = 1:numel(thisgrp)
            if any(thisgrp(iT) == grpscon(:,1))
                blnTrans = 0; break
            else
                blnTrans = 1;
            end
        end
        if blnTrans == 0
            break
        else
            grpscon = [grpscon; thisgrp ones(numel(thisgrp),1)*grpctr]; % save this...
            grpctr = grpctr + 1;
        end
    end
end

grpscon(1,:) = []; % remove padding zeros

% sort into ID order, and return just grouping
[~,ixsrt] = sort(grpscon(:,1),'ascend');
grpscon = grpscon(ixsrt,2); % return just the grouping IDs

varargout{1} = theta;