function [thresh,eta_max] = otsu1D(P)

% OTSU1D Otsu's method for finding a threshold on bimodality in histograms
% T = OTSU1D(P) given the L-length histogram of discrete probabilities in P,
% iteratively determines the bin number T at which the between-class variance
% is maximised. T is then the most likely threshold between the two modes
% (if the modes are clearly separated): bins(1:T) are mode 1; bins(T+1:L)
% are mode 2.
% 
% Reference:
% Otsu, N (1979) A Threshold Selection Method from Gray-Level Histograms, 
%   IEEE Trans Sys Man Mach, 9, 62-66 
%
% Change log:
% 18/05/2018 Initial version
% 21/05/2018 Directly implemented from Otsu's 1979 paper
%
% Mark Humphries 

T = sum(P);
L = numel(P);
if T ~= 1 P = P ./ sum(P); end
T = sum(P);

% global mean  - Otsu Eq 8
muT = sum((1:L).*P);

% global variance - Otsu Eq 15
var_T = sum(((1:L)-muT).^2 .* P);

maximum = 0;
thresh = 1;
w_k = 0; 
mu_k = 0;

for t = 1:L-1    % advance threshold
    w_k = w_k + P(t);   % cumulative probability of class 1 (Eq 6)
    mu_k = mu_k + t*P(t);   % cumulative mean of class 1 (Eq 7)
    
    var_b = (muT * w_k - mu_k).^2 / (w_k*(1-w_k));      % between class variance (Eq 18)
     
    if var_b > maximum
        thresh = t;
        maximum = var_b;
        eta_max = var_b / var_T;
    end
end
