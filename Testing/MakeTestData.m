function S = MakeTestData(W)

% MAKETESTDATA make test similarity matrix
% S = MAKETESTDATA(W) creates a synthetic similarity matrix S 
% which has a two-layer hierarchy, given the parameters in struct W:
%       W.N   : number of groups at the bottom of the hierarchy (even number)
%       W.size : size of each group at bottom of hierarchy (S is (size  x N))
%       W.within1.m : mean weight within the 1st layer groups 
%       W.within1.s : SD weight within the 1st layer groups
%       W.within2.m : mean weight within the 2nd layer groups (and between
%       the 1st layer groups)
%       W.within2.s : SD weight within the 2nd layer groups
%       W.between.m : mean weight between the 2nd layer groups
%       W.between.s : sd weight between the 2nd layer groups
%        
% Constructs a two-layer hierarchy:
%       - W.N groups at the bottom, connected within by .within1.
%       - 2 groups above (each W.N/2 in size), connected within by
%       .within2.
%       
%
% 26/11/2018: initial version
% Mark Humphries

% limits on entries - here assuming similarity in [0,1]
l = 0;
u = 1;

% catch errors
if mod(W.N,2) ~=0 
       error('Even number of groups needed')
end

% warn of unusual choices
if W.within1.m < W.within2.m
    warning('Internal weight of layer 1 smaller the between weight of layer 1')
end

if W.within1.m < W.between.m
    warning('Internal weight of layer 1 smaller the between weight layers')
end

%% initialise
n = W.size * W.N;
S = zeros(n);

%% build layer 1
maskWithin1 = zeros(n); 
% indices of within layer 1
for i = 1:W.N
    ixs = (W.size *(i-1)+1):i*W.size;
    maskWithin1(ixs,ixs) = 1; 
end
maskWithin1up = triu(maskWithin1,1);

% sample truncated Gaussian weights
% X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
L = (l - W.within1.m) / W.within1.s; U = (u - W.within1.m) / W.within1.s;
Nnz = sum(maskWithin1up(:));
z = trandn(ones(Nnz,1)*L,ones(Nnz,1)*U);
Z = W.within1.m +  W.within1.s * z;

S(maskWithin1up == 1) = Z;

%% build layer 2
maskWithin2 = zeros(n); 
for i = 1:2
    ixs = n/2 * (i-1)+1:i*n/2;
    maskWithin2(ixs,ixs) = 1; 
end
maskBetween = ones(n) - maskWithin2; % between is everything not within
maskWithin2 = maskWithin2 - maskWithin1; % now substract off within

% upper triangular only
maskWithin2up = triu(maskWithin2,1);

L = (l - W.within2.m) / W.within2.s; U = (u - W.within2.m) / W.within2.s;
Nnz = sum(maskWithin2up(:));
z = trandn(ones(Nnz,1)*L,ones(Nnz,1)*U);
Z = W.within2.m +  W.within2.s * z;

S(maskWithin2up == 1) = Z;

%% between layer 2
maskBetweenup = triu(maskBetween,1);

L = (l - W.between.m) / W.between.s; U = (u - W.between.m) / W.between.s;
Nnz = sum(maskBetweenup(:));
z = trandn(ones(Nnz,1)*L,ones(Nnz,1)*U);
Z = W.between.m +  W.between.s * z;

S(maskBetweenup == 1) = Z;

%% make symmetric, and put 1s on diagonal!
S = S + S';  % make symmetric

S(eye(n)==1) = 1;