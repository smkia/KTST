function [s] = statfun_ktst(cfg, dat, design)

% STATFUN_ktst is a function for computing a statistic for the relation
% between biological data and a design vector containing trial
% classifications or another independent variable
% It performs kernel two-sample test for the given two conditions.
%
% This function is called by STATISTICS_RANDOM, where you can specify
% cfg.statistic = 'ktst' which will be evaluated as statfun_ktst.
%
% The external interface of this function has to be
%   [s] = statfun_ktst(cfg, dat, design);
% where
%   dat    contains the biological data, Nvoxels x Nreplications
%   design contains the independent variable,  1 x Nreplications
%
% Additional settings can be passed through to this function using
% the cfg structure.
%     cfg: config structure that could have the following fields:
%           - iterations: Specifies number of iterations for permutation test. 
%             the default value is 10000.
%           - kernelType: string that specifies type of kernel, could be either
%             'linear' for linear kernel or 'gaussian' for gaussian kernel.
%             the default is 'gaussian'. Using linear kernel is not recommended. 
%           - kernelParam: Parameter needed for the specified kernel. In the case of
%             gaussian kernel the user can specify the sigma^2. For
%             example cfg.kernelPraram = 1. The default value is
%             equal to median of squared oaired distances.
%                          
% The output:
%               s.stat: Monte carlo p_value resulted from permutation
%               test.
% 
% Developed by Emanuele Olivetti (olivetti@fbk.eu) and Seyed Mostafa Kia
% (m.kia83@gmail.com), Thomas Hartmann (thomas.hartmann@th-ht.de), November, 2013.
% References: "The Kernel Two-Sample Test vs. Brain Decoding", 2013, 
% IEEE Proceedings of 3rd International Workshop on Pattern Recognition in Neuroimaging
% DOI: 10.1109/PRNI.2013.41

% Initialization
ft_defaults;
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig

% Setting the default configs
cfg.iterations = ft_getopt(cfg,'iterations',10000);
cfg.kernelType = ft_getopt(cfg,'kernelType','gaussian');

% Employing information in the design matrix
selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  error('inappropriate design, should contain 1''s and 2''s');
end
X = dat(:,selA);
Y = dat(:,selB);
% Here on, X stands for trials related to first condition and Y stands
% for trials of second condition.
m = size(X,1);
n = size(Y,1);
if size(X,2) ~= size(Y,2)
    error('two classes smust have same number of features.')
else
    d = size(X,2);
end

XY = [X;Y];
switch cfg.kernelType
    case 'linear'
        K = XY*XY';
    case 'gaussian'
        squared_distance_matrix = pdist2(XY,XY,'euclidean');
        if ~isfield(cfg,'kernelParam')
            cfg.kernelParam = median(squared_distance_matrix(:));
        end;
        sigma2 = cfg.kernelParam;
        K = exp(-squared_distance_matrix / sigma2);
    otherwise
        error('Wrong kernel type.');
end
% Computing MMD statistics
mmd2u = MMD2u(K, m, n);
% Computing null distribution
mmd2u_null = compute_null_distribution(K, m, n, cfg.iterations);

p_value = max(1/cfg.iterations, sum(mmd2u_null > mmd2u)/ cfg.iterations);
s.stat = mmd2u;
s.prob = p_value;
end

function [mmd] = MMD2u(K,m,n)
% This subfunction computes the MMD. See the reference for more
% information.
%
% Inputs:
%   - K: the kernel
%   - m: Number of trials in first condition
%   - n: Number of trials in second condition
% Output: 
%   - mmd: the computed MMD

Kx = K(1:m, 1:m);
Ky = K(m+1:end, m+1:end);
Kxy = K(1:m, m+1:end);
mmd =  1 / (m * (m - 1)) * (sum(sum(Kx)) - sum(diag(Kx))) + ...
    1 / (n * (n - 1)) * (sum(sum(Ky)) - sum(diag(Ky))) - ...
    2 / (m * n) * (sum(sum(Kxy)) - sum(diag(Kxy)));
end

function [mmd2u_null] = compute_null_distribution (K, m, n, iterations)
% This subfunction computes the null distribution using permutation.
%
% Inputs:
%   - K: the kernel
%   - m: Number of trials in first condition
%   - n: Number of trials in second condition
%   - iterations: number of iterations for computing null distribution.
% Output: 
%   - mmd2u_null: the computed null-distribution

mmd2u_null = zeros(1,iterations);
for i=1 : iterations
    idx = randperm(m+n);
    K_i = K(idx, idx);
    mmd2u_null(i) = MMD2u(K_i,m, n);
end
end