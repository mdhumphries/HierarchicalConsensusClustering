function E = nullmodelConsensusExpectation(C)

% NULLMODELCONSENSUSEXPECTATION expectation null model for consensus matrices
% E = NULLMODELCONSENSUSEXPECTATION(C) constructs the expected null model for a
% consensus matrix of random clusterings from its expectation, given:
%       C: the consensus matrix (each entry a proportion)
%
% Returns:
%   E: the N*N matrix of expected proportions of clusterings
%
% Reference:
%   Betzel, R. F.; Medaglia, J. D.; Papadopoulos, L.; Baum, G. L.; Gur, R.; Gur, R.; Roalf, D.; 
%   Satterthwaite, T. D. & Bassett, D. S. (2017) 
%   The modular organization of human anatomical brain networks: Accounting for the cost of wiring 
%   Network Neuroscience, in press.
%
% Mark Humphries 7/3/2017

N = size(C,1);

T = 2*sum(sum(triu(C,1))); % total proportion of times each pair was clustered together

E = zeros(N) + T/(N*(N-1));  % expectation of the proportion
E(eye(N)==1) = 0;           % no self-connections
