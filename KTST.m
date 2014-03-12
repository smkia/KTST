function [ p_value ] = KTST(X,Y,cfg)
% KTST performs kernel two-sample test for the giving two datasets X and Y.
% The inputs:
%               X: the first dataset with observation in rows and variables in
%               columns.
%               Y: the second dataset with observation in rows and variables in
%               columns. Number of variables should be equal to the number
%               of variables in X.
%               cfg: config structure that could have the following fields:
%                   - iterations: Specifies number of iterations for permutation test.
%                   the default value is 10000.
%                   - kernelType: string that specifies type of kernel, could be either
%                   'linear' for linear kernel or 'gaussian' for gaussian kernel.
%                   the default is 'gaussian'. Using linear kernel is not recommended.
%                   - kernelParam: Parameter needed for the specified kernel. In the case of
%                   gaussian kernel the user can specify the sigma^2. For
%                   example cfg.kernelPraram = 1. The default value is
%                   equal to median of squared oaired distances.
%                   - verbose: could be 0 or 1.
%
% The output:
%               p_value: Monte carlo p_value resulted from permutation
%               test.
% example: [ p_value ] = KTST(X,Y,10000,'gaussian',0)
% Developed by Emanuele Olivetti (olivetti@fbk.eu) and Seyed Mostafa Kia
% (m.kia83@gmail.com), July, 2013.
% References: "The Kernel Two-Sample Test vs. Brain Decoding", to appear,
% IEEE Proceedings of 3rd International Workshop on Pattern Recognition in Neuroimaging
% DOI: 10.1109/PRNI.2013.41

if nargin < 3
    error('Not enough inputs.');
end
if ~isfield(cfg,'iterations')
    cfg.iterations = 10000;
end;
if ~isfield(cfg,'kernelType')
    cfg.kernelType = 'gaussian';
end;
if ~isfield(cfg,'verbose')
    cfg.verbose = 0;
end;

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
            dists = squared_distance_matrix.^2;
<<<<<<< HEAD
            dists = tril(dists,-1);
            dists = reshape(dists,size(XY,1)^2,1);
            cfg.kernelParam = sqrt(0.5*median(dists(dists>0)));
        end;
        sigma = cfg.kernelParam;
        K = exp(-(squared_distance_matrix.^2)/2 / sigma^2);
=======
            dists = dists-tril(dists);
            dists=reshape(dists,size(XY,1)^2,1);
            cfg.kernelParam = sqrt(0.5*median(dists(dists>0)));
            %cfg.kernelParam = median(squared_distance_matrix(:));
        end;
        sigma = cfg.kernelParam;
        K = exp(-(squared_distance_matrix.^2)/2 / (sigma^2));
>>>>>>> 9f11cd5ad3947ccfaf7f8ac55d2073575345ca24
    otherwise
        error('Wrong kernel type.');
end

mmd2u = MMD2u(K, m, n);
if strcmp(cfg.nullDist,'permutation')
    mmd2u_null = compute_null_distribution(K, m, n, cfg.iterations);
elseif strcmp(cfg.nullDist,'estimation')
    [mmd2u_null] = nullDistEst(K,m,n,cfg.iterations);
end
if cfg.verbose
    figure;
    hist(mmd2u_null, 50);
end
p_value = max(1/cfg.iterations, sum(mmd2u_null > mmd2u)/ cfg.iterations);
end

function [mmd] = MMD2u(K,m,n)
Kx = K(1:m, 1:m);
Ky = K(m+1:end, m+1:end);
Kxy = K(1:m, m+1:end);
mmd =  (1 / (m * (m - 1)) * (sum(sum(Kx)) - sum(diag(Kx)))) + ...
    (1 / (n * (n - 1)) * (sum(sum(Ky)) - sum(diag(Ky)))) - ...
<<<<<<< HEAD
    (2 / (m * n) * sum(sum(Kxy)));
=======
    (2 / (m * n) * (sum(sum(Kxy))));
>>>>>>> 9f11cd5ad3947ccfaf7f8ac55d2073575345ca24
end

function [mmd2u_null] = compute_null_distribution (K, m, n, iterations)
mmd2u_null = zeros(1,iterations);
for i=1 : iterations
    idx = randperm(m+n);
    K_i = K(idx, idx);
    mmd2u_null(i) = MMD2u(K_i,m, n);
end
end

function [mmd2u_null] = nullDistEst(K,m,n,numNullSamp)
numEigs = m+n-2;
H = eye(m+n) - 1/(m+n)*ones(m+n,m+n);
Kz = H*K*H;
kEigs = eigs(Kz,numEigs);
kEigs = abs(kEigs)/(m+n); 
mmd2u_null = zeros(1,numNullSamp);
for i = 1 : numNullSamp
  mmd2u_null(i) = 2*sum(kEigs.*(randn(length(kEigs),1).^2));
end
mmd2u_null = mmd2u_null./(m+n); 
end