function CI = CIfromSEM(SD,N,I)

% CIfromSEM computes confidence intervals from SEM estimates
% CI = CIfromSEM(S,N,I) computes confidence intervals for the mean, given:
%       S: an n-length array of standard deviations
%       N: an n-length array of the number of samples
%       I: the set of percentiles for the confidence intervals to calculate (length m)
%               (e.g. I = [0.95,0.99]);
% Returns:
%       CI: an n x m array, giving the confidence interval for each 
%           supplied standard deviation, at each requested interval percentile
%
% Notes:
%       Computes interval using t-distribution, assuming SD is a sample 
%
% Mark Humphries 1/8/2016

if length(SD) ~= length(N)
    error('Standard deviation and sample number arrays do not match')
end

Vs = 1-(1-I)/2;  % confidence interval is symmetric around the mean

tvalue = tinv(Vs,N);

CI = tvalue .* SD ./ sqrt(N);
