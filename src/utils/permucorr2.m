function [r,p,ci,stats,dist] = permucorr2(x,varargin)
%PERMUCORR  Linear or rank permutation-based correlation (Pearson, Spearman, Rankit, Kendall τ-b).
%   R = PERMUCORR(X) returns a matrix containing the pairwise linear
%   correlation coefficients between each pair of columns in X based on
%   Pearson's r. For nonlinear correlations, set 'type' to 'spearman' or
%   'rankit' to transform raw data into ranks or rankits before applying
%   the Pearson formula. Set 'type' to 'kendall' to compute Kendall's tau-b.
%
%   R = PERMUCORR(X,Y) returns the pairwise correlation coefficient between
%   vectors X and Y. X and Y must have the same length. If X and Y are
%   matrices, the correlation coefficients between each corresponding pair
%   of columns in X and Y are calculated, and a vector of results is
%   returned.
%
%   [R,P] = PERMUCORR(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distribution that
%   the data come from. When only one sample is entered in X, two-tailed
%   permutation tests are automatically used. P-values are automatically
%   adjusted for multiple comparisons using the max correction method.
%
%   [R,P,CI] = PERMUCORR(...) returns a 100*(1-ALPHA)% confidence interval
%   (CI) for each coefficient. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [R,P,CI,STATS] = PERMUCORR(...) returns a structure with:
%       'df'        -- degrees of freedom (Pearson/Spearman/Rankit: n-2;
%                      Kendall: n*(n-1)/2 pairs)
%
%   [R,P,CI,STATS,DIST] = PERMUCORR(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUCORR(...,'PARAM1',VAL1,...) parameters:
%       'alpha'     Significance level (default=0.05).
%       'dim'       Work along columns (1, default) or rows (2).
%       'tail'      'both'(default) | 'right' | 'left'
%       'type'      'pearson'(default) | 'spearman' | 'rankit' | 'kendall'
%       'nperm'     Number of permutations (default=10000 or all perms if n<8).
%       'correct'   FWER control via max correction, true/false (default=true).
%       'rows'      'all' (default) | 'complete'
%       'seed'      RNG seed (default: time-based).
%       'verbose'   1 to print progress, 0 otherwise (default=1).
%
%   See also CORR CORRCOEF PARTIALCORR TIEDRANK.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS
%
%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%       [2] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [4] Bishara AJ, Hittner JB, (2012) Testing the Significance of a
%           Correlation With Nonnormal Data: Comparison of Pearson,
%           Spearman, Transformation, and Resampling Approaches. Psychol
%           Methods, 17(3):399-417.
%       [5] Bishara AJ, Hittner JB, (2017) Confidence intervals for
%           correlations when data are not normal. Behav Res, 49:294-309.
%
%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin<2 || ischar(varargin{1})
    y = [];
else
    y = varargin{1};
    varargin = varargin(2:end);
end

% Parse input arguments
arg = ptparsevarargin(varargin);

% Validate input parameters (ensure your ptvalidateparamin accepts 'kendall')
ptvalidateparamin(x,y,arg)

% Orient data column-wise
if arg.dim==2 || isrow(x)
    x = x';
end
if ~isempty(y) && (arg.dim==2 || isrow(y))
    y = y';
end

% Set up comparison
if isempty(y)
    warning('Comparing all columns of X in a correlation matrix...')
    [x,y] = ptpaircols(x);
    arg.tail = 'both';
    arg.mat = true;
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(y),2),:);
        y = y(~any(isnan(y),2),:);
        y = y(~any(isnan(x),2),:);
        x = x(~any(isnan(x),2),:);
    case 'all'
        if any(isnan(x(:))) || any(isnan(y(:)))
            error('X or Y is missing values. Set ROWS to ''complete''.')
        end
end

% Get data dimensions
[nobs,nvar] = size(x);

% Degrees of freedom placeholder
if nargout > 3
    df = nobs-2; % Pearson/Spearman/Rankit default
end

% Data transform (Kendall: no transform; ties handled internally)
switch lower(arg.type)
    case 'spearman'
        x = tiedrank(x);
        y = tiedrank(y);
    case 'rankit'
        x = norminv((tiedrank(x)-0.5)/nobs);
        y = norminv((tiedrank(y)-0.5)/nobs);
    case 'kendall'
        % no transform
    case 'pearson'
        % no transform
    otherwise
        error('Unknown correlation type: %s', arg.type)
end

% Flag: Kendall τ-b?
isKendall = strcmpi(arg.type,'kendall');

% Compute observed statistic r
if ~isKendall
    % Pearson/Spearman/Rankit -> Pearson formula on transformed data if any
    sdxy = sqrt((sum(x.^2)-(sum(x).^2)/nobs).*(sum(y.^2)-(sum(y).^2)/nobs));
    mu   = sum(x).*sum(y)/nobs;
    r    = (sum(x.*y)-mu)./sdxy;
else
    % Kendall τ-b (ties corrected), by column
    r = zeros(1,nvar);
    for k = 1:nvar
        r(k) = corr(x(:,k), y(:,k), 'type','Kendall', 'rows','complete');
    end
    if nargout > 3
        df = nobs*(nobs-1)/2; % pair count
    end
end

if nargout > 1
    % Permutations
    rng(arg.seed);
    if nobs < 8
        arg.nperm = factorial(nobs);
        idx = perms(1:nobs)';   % each column is a permutation
        if arg.verbose
            warning('Computing all possible permutations due to small N.')
            fprintf('Number of permutations used: %d\n', arg.nperm);
        end
    else
        [~,idx] = sort(rand(nobs,arg.nperm)); % columns are permutations
    end

    % Progress display controls
    showProgress = (isfield(arg,'verbose') && ~isequal(arg.verbose,0));
    reportStep   = max(1, floor(arg.nperm/20));  % print ~every 5%

    % Estimate sampling distribution dist (rows: perms; cols: variables)
    dist = zeros(arg.nperm,nvar);

    if ~isKendall
        % Recompute for safety (x,y unchanged across perms)
        sdxy = sqrt((sum(x.^2)-(sum(x).^2)/nobs).*(sum(y.^2)-(sum(y).^2)/nobs));
        mu   = sum(x).*sum(y)/nobs;
        for i = 1:arg.nperm
            dist(i,:) = (sum(x(idx(:,i),:).*y)-mu)./sdxy;
            if showProgress && mod(i,reportStep)==0
                fprintf('Permutation %d / %d (%.1f%%)\n', i, arg.nperm, i/arg.nperm*100);
            end
        end
    else
        for i = 1:arg.nperm
            xi = x(idx(:,i),:);
            for k = 1:nvar
                dist(i,k) = corr(xi(:,k), y(:,k), 'type','Kendall', 'rows','complete');
            end
            if showProgress && mod(i,reportStep)==0
                fprintf('Permutation %d / %d (%.1f%%)\n', i, arg.nperm, i/arg.nperm*100);
            end
        end
    end

    if showProgress
        fprintf('Permutation complete. Computing p-values...\n');
    end

    % Max correction if specified
    if arg.correct
        switch arg.tail
            case 'both'
                [~,midx] = max(abs(dist),[],2);
                csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
                dist = dist';
                dist = dist(midx+csvar); % vector of max |stat| per perm
            case 'right'
                dist = max(dist,[],2);
            case 'left'
                dist = min(dist,[],2);
        end
    end

    % p-values & CI from empirical distribution
    switch arg.tail
        case 'both'
            p = 2*(min(sum(r<=dist),sum(r>=dist))+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha/2));
                ci = [max(-1,r-crit); min(1,r+crit)];
            end
        case 'right'
            p = (sum(r<=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha));
                ci = [max(-1,r-crit); Inf(1,nvar)];
            end
        case 'left'
            p = (sum(r>=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(-dist,100*(1-arg.alpha));
                ci = [-Inf(1,nvar); min(1,r+crit)];
            end
        otherwise
            error('Unknown tail option: %s', arg.tail)
    end
end

% Arrange results in a matrix if requested (pairwise corr matrix mode)
if isfield(arg,'mat') && arg.mat
    r = ptvec2mat(r);
    if nargout > 1
        p = ptvec2mat(p);
    end
    if nargout > 2
        ciLwr = ptvec2mat(ci(1,:));
        ciUpr = ptvec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
    end
    if nargout > 3
        df = ptvec2mat(df);
    end
end

% Stats struct
if nargout > 3
    stats.df = df;
end
end
