function [ p_value, h ] = smk_freq_KTST(cfg,TF1,TF2)
% smk_freq_KTST performs kernel two-sample test for the giving two time-frequency 
% representation of MEG data and it plots the topographic map of significant 
% channels if needed.
%
%
% TH: As pointed out in nathan's email, convert the function to fieldtrip's
% statistics framework. This would, for instance, include, that the data is
% provided in a different way. In fact, you could have the data either in
% only one structure (if you want to compare trials within a subject) our
% have (e.g. an averaged) freqanalysis structure for each subject. In the
% later case, the function would accept a variable number of input
% arguments. you can refer to the function "ft_statfun_depsamplesT" for a
% simple example, how to do it. This would allow us to use your function
% within the fieldtrip statistic framework, making it infinitly more
% usefull to us.
%
% TH: please add a definition of the function's signature also in the
% documentation. this would be the perfect place to do it. the problem is
% that matlab only shows what is in the comment although it would be very
% logical to include the function's signature as well.
%
% The inputs:
%   TF1: the first condition data. Actually this input is the output of 
%        ft_freqanalysis for the first condition (eg. faces). 
%   TF2: the second condition data. Actually this input is the output of 
%        ft_freqanalysis for the second condition (eg. houses).
%
%   cfg: config structure that could have the following fields:
%       - iterations: Specifies number of iterations for permutation test. 
%         the default value is 10000.
%       - kernelType: string that specifies type of kernel, could be either
%         'linear' for linear kernel or 'gaussian' for gaussian kernel.
%         the default is 'gaussian'. Using linear kernel is not recommended. 
% TH: what do you mean by "not recommended"? if people should not use it
% at all, i do not see any reason to include it here. i think, it would
% help people greatly evaluating what kernel to use if you provided some
% use cases for choosing, e.g. the linear kernel.
%       - kernelParam: Parameter needed for the specified kernel. In the case of
%         gaussian kernel the user can specify the sigma^2. For
%         example cfg.kernelPraram = 1. The default value is
%         equal to median of squared paired distances (the default value is recommended).
%       - alpha: critical alpha value. The default is 0.05.
%       - verbose: could be 0 or 1.
%       - freqLim: 1x2 vector that specifies frequencies of
%         interest. eg. [5,45].Default is all frequency bins ('all'). 
%       - timeLim: 1x2 vector that specifies timebins (in second) of
%         interest. eg. [-0.1,1.5]; Default is all timebins ('all').
% TH: do these two parameters above really accept 'all' as an input. this
% would be a great feature so i would love to see it added if it does not
% exist by now ;-) but it is not necessary.
%       - plot: could be 0 or 1. If 1 it plots topographic map. (Default is
%       0);
%       - layout: topographic map layout (it is required if plot is 1).
%       - parallel: could be 0 or 1. if 1 the function is executing in 
%         parallel manner on several cores. and then it is faster.   
%         The default is 0;
% TH: please add one sentence, how this is implemented. we use quite a
% variety of parallel environments here. so if you require a specific one
% it should be mentioned here.
%
% The output:
%   p_value: Monte carlo estimated p_value resulted from permutation
%            test.
%   h: FDR corrected significant channels.
%
% example: [ p_value, h ] = smk_freq_KTST(cfg,TF1,TF2) where cfg contains:
%           cfg.iteragtions = 10000;
%           cfg.kernelType = 'gaussian';
%           cfg.alpha = 0.05;
%           cfg.layout = 'CTF275.lay';
% TH: an example is always great. however, in this particular one i do not understand why you included the layout option
% as it you write earlier that it is only needed when plotting is turned
% on.
%
% Developed by Emanuele Olivetti (olivetti@fbk.eu) and Seyed Mostafa Kia
% (m.kia83@gmail.com), Thomas Hartmann (thomas.hartmann@th-ht.de), September, 2013.
% References: "The Kernel Two-Sample Test vs. Brain Decoding", to appear,
% IEEE Proceedings of 3rd International Workshop on Pattern Recognition in Neuroimaging
% DOI: 10.1109/PRNI.2013.41
%
% TH: please add a revision history here including the initials and date of
% changes.

% Initialization
ft_defaults;
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig

% Checking the inputs
if nargin < 3
    error('Not enough inputs.');
end
[TF1] = ft_checkdata(TF1,'feedback','yes','datatype',{'freq'},'dimord','rpt_chan_freq_time');
[TF2] = ft_checkdata(TF2,'feedback','yes','datatype',{'freq'},'dimord','rpt_chan_freq_time');

% Setting the default configs
cfg.iterations = ft_getopt(cfg,'iterations',10000);
cfg.kernelType = ft_getopt(cfg,'kernelType','gaussian');
cfg.verbose = ft_getopt(cfg,'verbose',0);
cfg.alpha = ft_getopt(cfg,'alpha',0.05);
cfg.parallel = ft_getopt(cfg,'parallel',0);
cfg.freqLim = ft_getopt(cfg,'freqLim','all');
cfg.timeLim = ft_getopt(cfg,'timeLim','all');
cfg.plot = ft_getopt(cfg, 'timeLim', 0);
if cfg.plot == 1 & ~isfield(cfg,'layout')
    error('Please specify the layout for topographic map.');
end

TF1.powspctrm(isnan(TF1.powspctrm)) = 0;
TF2.powspctrm(isnan(TF2.powspctrm)) = 0;
% TH: are you sure that you want to treat nans as 0? i think, you assume
% here that the nans in the two TF structures match. what, if they do not?

[trialNum1, channelNum, freqNum, timeNum] = size(TF1.powspctrm);
[trialNum2] = size(TF2.powspctrm,1);

% TH: please add a check and a gracefull error in case the dimensions of
% TF1.powspctrm and TF2.powspctrm do not match.

if ischar(cfg.freqLim) % TH: please also test for the right string and give an error if it is not 'all'.
     cfg.minFreqIdx = 1;
     cfg.maxFreqIdx = freqNum;
else
    cfg.minFreqIdx = find(round(TF1.freq) == cfg.freqLim(1));
    cfg.maxFreqIdx = find(round(TF1.freq) == cfg.freqLim(2));
end
if ischar(cfg.timeLim) % TH: see above
     cfg.minTimeIdx = 1;
     cfg.maxTimeIdx = timeNum;
else
    cfg.minTimeIdx = find(TF1.time == cfg.timeLim(1));
    cfg.maxTimeIdx = find(TF1.time == cfg.timeLim(2));
end

% TH: please do not add anything to the cfg structure. it is ok to
% substitute missing fields by defaults but adding does not help. fieldtrip
% "tracks" the cfg structure, i.e. it adds it to the processed data. so, it
% is important to keep it clean.

p_value = zeros(1,channelNum);

if cfg.parallel == 1
    matlabpool;
end

for i = 1 : channelNum
    if isfield(cfg,'kernelParam')
        sigma2 = cfg.kernelParam;
    end
    
    % Here on, X stands for trials related to first condition and Y stands
    % for trials of second condition.
    X = zeros(trialNum1,(cfg.maxFreqIdx-cfg.minFreqIdx+1)*(cfg.maxTimeIdx-cfg.minTimeIdx+1));
    Y = zeros(trialNum2,(cfg.maxFreqIdx-cfg.minFreqIdx+1)*(cfg.maxTimeIdx-cfg.minTimeIdx+1));
    for j = 1 : trialNum1
        temp = squeeze(TF1.powspctrm(j,i,cfg.minFreqIdx:cfg.maxFreqIdx,cfg.minTimeIdx:cfg.maxTimeIdx));
        X(j,:) = temp(:);
    end
    % TH: i think, these two loops can be vectorized. but that is not high
    % priority...
    for j = 1 : trialNum2
        temp = squeeze(TF2.powspctrm(j,i,cfg.minFreqIdx:cfg.maxFreqIdx,cfg.minTimeIdx:cfg.maxTimeIdx));
        Y(j,:) = temp(:);
    end
    XY = [X;Y];
    switch cfg.kernelType
        case 'linear'
            kernel = XY*XY';
        case 'gaussian'
            squared_distance_matrix = pdist2(XY,XY,'euclidean');
            if ~isfield(cfg,'kernelParam')
                sigma2 = median(squared_distance_matrix(:));
            end % TH: what happens if there is a cfg.kernelParam field? i think you
                % handle that case above. i would suggest doing it here for
                % readability....
            kernel = exp(-squared_distance_matrix / sigma2);
        otherwise
            error('Wrong kernel type.');
    end
    % Computing MMD
    mmd2u = MMD2u(kernel, trialNum1, trialNum2);
    % Computing null distribution
    mmd2u_null = compute_null_distribution(kernel, trialNum1, trialNum2, cfg.iterations);
%     if cfg.verbose
%         figure;
%         hist(mmd2u_null, 50);
%     end
    p_value(i) = max(1/cfg.iterations, sum(mmd2u_null > mmd2u)/ cfg.iterations);
    if cfg.verbose
        disp(strcat('Channel number:',num2str(i),'/',num2str(channelNum)));
    end
end
% Correcting for MCP
[h] = FDR_BH(p_value,cfg.alpha);
% Plotting topographic map
if cfg.plot == 1
    try
        topo_plot(TF1,TF2,h,cfg);
    catch
        disp('Plot error: Please add fieldtrip to your matlab paths.');
    end
end
if cfg.parallel == 1
    matlabpool close;
end

ft_postamble provenance
ft_postamble trackconfig
end

function topo_plot(TF1,TF2,h,c)
% Plots the topographic map. The values of map are based on the raweffect
% computed by subtracting result of ft_freqdeccriptives function for both
% conditions. The significant channels are highlighted by * signs.
cfg = [];
cfg.foilim = c.freqLim;
cfg.toilim = c.timeLim;
TF1_desc = ft_freqdescriptives(cfg, TF1);
TF2_desc  = ft_freqdescriptives(cfg, TF2);
raweffect = TF1_desc.powspctrm - TF2_desc.powspctrm;
plotFormat.label = TF1.label;
plotFormat.freq = TF1.freq(c.minFreqIdx:c.maxFreqIdx);
plotFormat.dimord = 'chan_freq_time';
plotFormat.time = TF1.time(c.minTimeIdx:c.maxTimeIdx);
plotFormat.cfg = TF1.cfg;
plotFormat.grad = TF1.grad;
plotFormat.powspctrm = raweffect;
plotFormat.raweffect = raweffect;
plotcfg = [];
plotcfg.interpolation = 'v4';
plotcfg.layout = c.layout;
plotcfg.comment = 'no';
plotcfg.parameter = 'raweffect';
plotcfg.highlight = 'on';
plotcfg.highlightchannel =  plotFormat.label(logical(h),1);
plotcfg.highlightsymbol = '*';
plotcfg.highlightcolor =  [0 0 0]; %(black))
plotcfg.highlightsize = 6;
plotcfg.highlightfontsize = 8;
plotcfg.marker = 'off';
ft_topoplotTFR(plotcfg,plotFormat);
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
parfor i=1 : iterations
    idx = randperm(m+n);
    K_i = K(idx, idx);
    mmd2u_null(i) = MMD2u(K_i,m, n);
end
end

function [h] = FDR_BH(pValues,alpha)
% This function receives the p_values and critival alpha and returns MCP corrected results h.
% Here the BH (Benjamini and Hochberg) method is used for controlling FDR. For reference see:
% Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
% rate:Apractical and powerful approach tomultiple testing. Journal of
% the Royal Statistical Society. Series B (Methodological), 57, 289–300.

dim = size(pValues);
pValues = pValues(:);
[pValuesSorted, sortedIndices] = sort(pValues);
testNum = length(pValues);
thresh = ((1:testNum)/testNum)  * alpha;
h = (pValuesSorted<=thresh');
[~, unsort] = sort(sortedIndices);
h = h(unsort);
h = reshape(h, dim);
h = double(h);
end