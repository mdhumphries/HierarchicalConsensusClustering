function X = discreteinvrnd(p,m,n)

% p = vector of discrete probabilities; m,n = size of returned matrix

X = zeros(m,n); % Preallocate memory
cum_p = cumsum(p);
for i = 1:m*n
    u = rand;
    X(i) = find(u < cum_p,1);
    % X(i) = min(I);
end
