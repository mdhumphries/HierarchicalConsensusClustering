function [grpscon,ctr,k,varargout] = ConsensusSweep(S,B,varargin)
 
% CONSENSUSSWEEP consensus partition using k-means sweep, with multiple consensus matrices
%   [C,N,K] = CONSENSUSSWEEP(S,B) finds the consensus-clustering 
%   partition of the n*n affinity matrix S, using k-means sweep between bounds in B = [l u]
%
%   Returns: 
%       C: column vector indicating group membership for the consensus
%       partition(s): (n rows, by c columns - one column per tested K)
%       N:  the number of iterations until consensus was reached. 
%       K: the number of clusters in each returned partition, per K tested
%
%   ...= CONSENSUSSWEEP(...,OPTS) 
%           OPTS.nreps : sets k-means to run N times for each specified group size (default is 50)
%           OPTS.dims : 'scaled' (default), 'all'; the embedding dimensions for each tested K. See KMEANSSWEEP
%           OPTS.project : 'Laplacian' (default), 'Eigs'
%           OPTS.stop : 'first' (default),'all' ; stop when first consensus
%           matrix converges; or stop when all have converged (or escape
%           reached)
%           OPTS.escape : number of consensus iterations before exiting
%           (default: 10)
%   
%   [...,bCON,CLU] = CONSENSUSSWEEP(...) where:
%       bCON : a Boolean flag indicating whether the consensus algorithm
%                converged or not
%       CLU:   a matrix of every single clustering of the similarity matrix S in the first
%              pass (i.e. before the consensus) - this is useful for further
%              post-processsing.
%
%   Notes: 
%
%   (1) Sequence of algorithm:
%       - Project data: Laplacian random walk of similarity matrix
%       (default), or eigenvectors of S
%       - use k-means on projected objects, per requested k; dimensions of
%           projection scale with k
%       - Make consensus per k [C(k)]
%       - Check if any/all C(k) have converged (depending on option flags); if so, exit
%       - If not: 
%              - Project each C(k)
%              - Cluster each C(k) using k
%              - Make new consensus C(k)
%        Until converged
%
%   References: 
%   (1) Luxburg
%
%   20/11/2018 : initial version, derived from CONSENSUSSWEEPSINGLE
%
%   Mark Humphries 

% defaults
options.nreps = 50;     % of each k
options.dims = 'scale';   % scale embedding dimensions to number of clusters for each k-means clustering
options.project = 'Laplacian';
options.stop = 'first';
options.escape = 10;


% check for errors
if sum(S) ~= sum(S')
    error('Similarity matrix is not symmetric');
end

%% set up options
if nargin >= 3
    if isstruct(options) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            options = setfield(options,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end


% % set up saving of each iteration
% blnSave = 0; % internal flag for setting saving of data
% if blnSave  
%     fname = ['Consensus_Iterations_' datestr(now,30)];
%     save(fname,'dims','nreps','blnExplore');  % save initial data to allow -append to work below
% end

%% internal parameters
ks = B(1):B(2);        % range of cluster sizes to test
kept_ks = ks;          % set of kept cluster sizes (dropping out as each converges) 

nIDs = size(S,1);     % number of nodes of the similarity matrix
blnConverged = 0;       % stopping flag
grpscon = zeros(nIDs,numel(ks));        % keep converged grps per K
ctr = 0;                % iterations of consensus

%% do initial projection and clustering

% project the affinity matrix
switch options.project
    case 'Laplacian'
        [D,~] = ProjectLaplacian(S,B(2)+1);  % return all dimensions to upper bound (add one more dimension for first unit vector)

    case 'Eigs'
        D = ProjectEigs(S,B(2));  % return all dimensions to upper bound
    otherwise
        error('Unknown option for projecting data')
end

% get partitions using k-means
C = kmeansSweep(D,B(1),B(2),options.nreps,options.dims);
varargout{2} = C;  % store initial clusterings

%% loop: consensus matrix and its clustering until converged

while ~blnConverged
    ctr = ctr+1  % increment consensus iteration count            
    if  ctr > options.escape
        % do escape if not converging
        warning(['Did not converge in ' num2str(options.escape) 'iterations - exiting without consensus answer'])
        blnConverged = 0;
        return
    end

    %% now make consensus
    % loop over tested k
    allCCons = zeros(nIDs,nIDs,numel(kept_ks)); % empty consensus matrix
    blnVec = zeros(numel(kept_ks),1); theta = blnVec; thisgrps = cell(numel(kept_ks),1);
    for iK = 1:numel(kept_ks)
        % make consensus from corresponding part of C
        ix = 1+options.nreps*(iK-1):options.nreps*iK; % indices into C
        % thisS = allSum(ix);
        allCCons(:,:,iK) = makeConsensusMatrix(C(:,ix));  % make consensus

        % test for convergence, returning theta
        [blnVec(iK),thisgrps{iK},theta(iK)] = CheckConvergenceConsensus(squeeze(allCCons(:,:,iK)));   
    end

    % check what we should do if anything converged
    ixConv = find(blnVec);  % index of converged K (in set of kept K)
    % keyboard
    if any(blnVec)
        % get converged answers
        for iC = 1:numel(ixConv)
            try
            if max(thisgrps{ixConv(iC)}) ~= kept_ks(ixConv(iC))
                warning(['Converged groups do not match number of tested groups at k=' num2str(kept_ks(ixConv(iC)))]);
                % keyboard
            end
            catch
                keyboard
            end
            
            ixK = kept_ks(ixConv(iC)) - B(1) + 1;   % map this K value onto the index for storing groups
            grpscon(:,ixK) = thisgrps{ixConv(iC)};  % store converged group

        end
        k = max(grpscon);  % return number of clusters in each partition
        
        % check if we're stopping yet
        kept_ks(ixConv) = [];  % update set of retained Ks, remove all that converged
        switch options.stop
            case 'first'
                blnConverged = 1;       % if any converged, then stop
            case 'all'                  % wanting all to converge
                if isempty(kept_ks)
                    blnConverged = 1;
                else
                    blnConverged = 0;
                end
        end

    end
    
    if ~blnConverged
        % keyboard
        C = zeros(nIDs,numel(kept_ks) * options.nreps); % clear C!
        for iK = 1:numel(kept_ks) % checking all kept K, shifting starting index to 1, so that iK counts within the complete range of K
            % for each value of K we're still checking...
            % project the consensus matrix
            switch options.project
                case 'Laplacian'
                    [D,~] = ProjectLaplacian(squeeze(allCCons(:,:,iK)),kept_ks(iK)+1);  % return all dimensions for k

                case 'Eigs'
                    try
                    D = ProjectEigs(squeeze(allCCons(:,:,iK)),kept_ks(iK));  % return all dimensions for k
                    catch
                        keyboard
                    end
               
                otherwise
                    error('Unknown option for projecting data')
            end

            if any(imag(D(:)))
                keyboard
            end

            % get partitions of that consensus matrix using k-means
            ix = 1+options.nreps*(iK-1):options.nreps*iK; % indices into C
            C(:,ix) = kmeansSweep(D,kept_ks(iK),kept_ks(iK),options.nreps,options.dims);  % repeat the same k
        end
    end
end

%     if blnSave
%         eval(['CCons', num2str(ctr),' = CCons;']);  % store consensus matrix
%         eval(['Groups', num2str(ctr),' = grpscon;']); 
%         save(fname,['CCons', num2str(ctr)],['Groups', num2str(ctr)],'-append');
%     end
            
% overwrite convergence flag if successful
varargout{1} = blnConverged;
    
     
function D = ProjectEigs(S,M)

% find top M eigenvectors of similarity matrix S
    [V,egs] = eig(S,'vector');
    [~,ix] = sort(egs,'descend');       % make sure eigenvalues are sorted
    D = V(:,ix(1:M));  % return all eigenvectors to upper bound




  







