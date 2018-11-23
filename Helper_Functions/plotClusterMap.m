function [h,hc,I] = plotClusterMap(W,G,varargin)

% PLOTCLUSTERMAP sorted heatmap of detected clusters
% [H,C,I] = PLOTCLUSTERMAP(W,G) for a set of N data objects, plots the NxN comparison matrix W
% (e.g. correlation or similarity matrix) as a heatmap grouped into
% clusters G, an N-element array of group membership indexes.
%
% ... = PLOTCLUSTERMAP(...,C,L,FLAG) are a set of optional arguments:
%       C = nx3 array of colormap values for the cluster map; set to [] to omit. 
%       L = 1x3 array of colormap values for the lines boxing each cluster. 
%       FLAG : a string of options:
%           'S' : sorts the groups by intrinsic similarity before plotting
%           'N' : turn off the colorbar 
%
% Returns:
%    H: the axis handle
%    C: the colorbar handle
%    I: the sorted indices of the nodes in left-to-right (& top-to-bottom)
%
% 28/7/2016: initial version
% 30/8/2016: added colormaps for different situations; plot sorted by
%            similarity
% 23/9/2016: plots within groups in order of similarity; increased resolution of heatmap   
% 15/8/2017: returned handles, add line colour option, updated Help
%
% Mark Humphries 28/7/2016

lineColor = [0 0 0];
blnCBar = 1;

if nargin > 2 && ~isempty(varargin{1})
    % specify colormap
    cmap = varargin{1}; 
else
    % black and white colormap for binary matrices
    if numel(unique(W)) == 2
        cmap = [1 1 1; 0 0 0];
    else
        % make a colormap from C Brewer schemes
        cmap = brewermap(15,'*Blues');  % reversed schemes: light = high
    end
end

if nargin > 3 && ~isempty(varargin{2})
    % specified linecolor
    lineColor = varargin{2}; 
else
    if numel(unique(W)) == 2
        % use red to outline binary matrix clusters
        lineColor = [0.8 0.3 0.4];
    else
        % use white to complement default blue colormap
        lineColor = [1 1 1];
    end
end


% sorting of nodes, and other options
[srt,I] = sort(G,'descend'); % sort into group index order... (group all 1s, 2s etc)
if nargin > 4
     if strfind(varargin{3},'S')
        G = [[1:numel(G)]' G];
        [newG,Sgrp,Sin,~] = sortbysimilarity(G,W); % sort clusters by similarity (new G)
     %  keyboard  
        % sort by group similarity, and within that by node similarity
        maxG = max(newG(:,2));
        I = []; srt = [];
        for iG = 1:numel(Sgrp)
            % get each group, starting with most similar
            ixNodes = find(newG(:,2) == maxG - (iG-1)); % nodes in this group (original indexes)
            % sort in similarity order
            theseSin = Sin(ixNodes);
            [~,ixSin] = sort(theseSin,'descend');
            % add to group array
            srt = [srt; maxG-(iG-1)*ones(numel(ixNodes),1)];  % sorted descending index of groups
            % store remapped indices
            I = [I; ixNodes(ixSin)]; % original node indices in order
            % keyboard
        end
        
%         % then sort into group order
%         [srt,I] = sort(newG(:,2),'descend'); % so that the most similar group is first

        % keyboard
     end
     if strfind(varargin{3},'N')
            blnCBar = 0;
     end
end

lines = [0; find(abs(diff(srt))==1); numel(srt)]+0.5;  % absolute difference, so it doesn't account

figure
h = imagesc(W(I,I));
colormap(cmap);
axis square
if blnCBar 
    hc = colorbar; 
else
    hc = [];
end

% keyboard
% draw outline box around each cluster
for i=2:numel(lines)
    line([lines(i-1) lines(i) lines(i) lines(i-1); lines(i) lines(i) lines(i-1) lines(i-1)],...
         [lines(i-1) lines(i-1) lines(i) lines(i); lines(i-1) lines(i) lines(i) lines(i-1)],...,
         'Color',lineColor,'LineWidth',2)
end