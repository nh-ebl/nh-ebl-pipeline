% Function takes in a 1-D 'x' and 1-D 'y' with an optional 'w' weights vector of the same size as 'x' and 'y'.
% Returns 'fitter' struct that has the fit coefficients as well as quality-of-fit parameters.

function fitter = linear_fit(x,y,w)

n = length(x); % Get the length of the data
if( ~exist('w','var') )
    w = ones(n,1); % Set weights to 1
end

% [fitobj,gof,output] = fit(x,y,'poly1','Weights',w);

% Weighted least squares linear regression follows eqs from https://ms.mcmaster.ca/canty/teaching/stat3a03/Lectures7.pdf
% xWM = sum(w.*x)/sum(w); % Weighted mean of x
% yWM = sum(w.*y)/sum(w); % Weighted mean of y
% m = sum(w.*(x-xWM).*(y-yWM))/sum(w.*(x-xWM).^2); % Calc slope
% b = yWM-m*xWM; % Calc intercept

% Matrix version of weighted least squares linear regression
% From this guy who explained it all https://www.youtube.com/watch?v=aMJOgopkLeY
% Confirmed identical to above scalar method
X = [ones(n,1),x]; % X matrix, 1's column represents the constant b, other column represents m (b/c m*x I guess)
W = eye(n).*w; % Weight matrix
beta = (X'*W*X)^-1*X'*W*y; % Calculate beta estimators (holds m and b) ( from https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation )
m = beta(2); % Get out slope which is beta(2)
b = beta(1); % Get out intercept which is beta(1)
% Full derivation here https://stats.stackexchange.com/a/52712

% General fit quality parameters that apply to both weighted and non-weighted fits
%--- residuals ---
residuals = y - (m*x + b); % Deviation from true y, with the order being Real - Guessed

%--- Degrees of Freedom ---
DoF = n-2; % Degrees of freedom, which is just the length of the data minus the number of variables used in the estimator (2)

%--- Variance of residuals ---
% Un-weighted variance calcs here for reference, checked and if w is all 1's they match the weighted versions
% variB_r = sum((residuals - mean(residuals)).^2)/n;% Biased sample variance ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance )
% vari_r = sum((residuals - mean(residuals)).^2)/(n-1);% Unbiased sample variance [equiv to var(residuals)] ( https://en.wikipedia.org/wiki/Variance#Sample_variance )
% Used reliability weights because that's what the weighting is I think?
residualsWM = sum(w.*residuals)/sum(w); % Weighted mean of residuals
variBw_r = sum(w.*(residuals - residualsWM).^2)/sum(w); % Biased weighted sample variance [equiv to var(residuals,w)] ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance )
variw_r = variBw_r/(1 - (sum(w.^2)/sum(w).^2));% Unbiased weighted sample variance ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights ) & backed up here, confirmed identical ( https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic#Geochronology )
% Note frequency weights exist, it would be weighted by the number of counts ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Frequency_weights )

%--- Variance of y ---
% Un-weighted variance calcs here for reference, checked and if w is all 1's they match the weighted versions
% variB_y = sum((y - mean(y)).^2)/n;% Biased sample variance ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance )
% vari_y = sum((y - mean(y)).^2)/(n-1);% Unbiased sample variance [equiv to var(residuals)] ( https://en.wikipedia.org/wiki/Variance#Sample_variance )
% Used reliability weights because that's what the weighting is I think?
yWM = sum(w.*y)/sum(w); % Weighted mean of residuals
variBw_y = sum(w.*(y - yWM).^2)/sum(w); % Biased weighted sample variance [equiv to var(residuals,w)] ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance )
variw_y = variBw_y/(1 - (sum(w.^2)/sum(w).^2));% Unbiased weighted sample variance ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights ) & backed up here, confirmed identical ( https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic#Geochronology )
% Note frequency weights exist, it would be weighted by the number of counts ( https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Frequency_weights )


%--- Standard Deviation of residuals ---
% stdevBw = sqrt(variBw_r); % Weighted standard deviation [equiv to std(residuals,w)]
stdevw_r = sqrt(variw_r); % Weighted unbiased standard deviation
% Theoretically this should be equal to rmse_r == se_r, but it isn't while rmse_r == se_r is true
% For _m and _b, stdevw_m == se_m is true, rmse_m uncalc'd same for b

%--- Chi-Squared of residuals ---
% !! this chiSq seems to be wrongly using the variw instead of some other variance that it should be using (since Sw below does not involve a variance)
% chiSq_r = sum(residuals.^2/variw_r); % Chi-squared value ( https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic#Definition )
% chiSqReduc_r = chiSq_r/DoF; %Reduced chi-squared statistic
Sw_r = residuals'*W*residuals; % Minimum value of the weighted residuals ( https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation )
chiSqReducw_r = Sw_r/DoF; % Reduced chi squared ( https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic / https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation ) 

%--- Mean Squared Error of residuals ---
% Un-weighted MSE/MSE used in regression here for reference
% RSS_r = sum(residuals.^2); % residuals sum of squares (https://en.wikipedia.org/wiki/residuals_sum_of_squares)
% MSE_r = RSS_r/n; % Mean Squared Error ( https://en.wikipedia.org/wiki/Mean_squared_error )
% MSE_lsq_r = RSS_r/DoF; % Mean Squared Error referenced in regression - uses DoF instead of n ( https://en.wikipedia.org/wiki/Mean_squared_error#In_regression )
RSSw_r = sum(w.*residuals.^2); % Weighted residuals sum of squares ( https://ms.mcmaster.ca/canty/teaching/stat3a03/Lectures7.pdf )
MSEw_r = RSSw_r/DoF; % Weighted mean square error for regression using DoF ( https://ms.mcmaster.ca/canty/teaching/stat3a03/Lectures7.pdf )
% per https://en.wikipedia.org/wiki/Mean_squared_error
% for unbiased estimator, MSE is the variance

%--- Root Mean Squared Error of residuals ---
% per https://en.wikipedia.org/wiki/Mean_squared_error
% for unbiased RMSE == standard error ==
% standard deviation == square root of the unbiased variance
RMSEw_r = sqrt(MSEw_r); % Just the root of the MSEw, this is very close to the SE as well as the standard deviation of the residuals

%--- Jacobian ---
J = [x,ones(n,1)]; % Jacobian is based on the partials for y = m*x+b, with each column being a partial and each row being a data point
% dy/dm = x so column 1 is just x | dy/db = 1 so column 2 is just 1's

%--- Weighted Variance-Covariance Matrix of Estimated Parameters ---
% def wrong Mbeta = (X'*W*X)^-1*X'*W*M*W'*X*(X'*W'*X)^-1; % Estimated variance-covariance matrix of the parameter estimates ( https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation )
Mbeta = chiSqReducw_r*(X'*W*X)^-1; % Variance-Covariance matrix for the estimated parameters ( https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation )
% NOTE that chiSqReducw is the effective variance of y if Var(y) = variance*W^-1 per https://stats.stackexchange.com/a/52712
% all from https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation
variw_m = Mbeta(2,2); % Variance of slope, 2,2 since beta(2) = m
stdevw_m = sqrt(variw_m); % Standard deviation of slope, 2,2 since beta(2) = m
variw_b = Mbeta(1,1); % Variance of intercept, 1,1 since beta(1) = b
stdevw_b = sqrt(variw_b); % Standard deviation of intercept, 1,1 since beta(1) = b
corrcoef = Mbeta(1,2)/(stdevw_m*stdevw_b); % Correlation coefficient

%--- Weighted Variance-Covariance Matrix of residuals ---
% all from https://en.wikipedia.org/wiki/Weighted_least_squares#residuals_values_and_correlation
% not sure how to make this useful
% r = y - H*y; % residuals vector, same as var residuals

%!! This should be right, but it is commented out b/c it is very slow at the moment
% Hat = X*(X'*W*X)^-1*X'*W; % "hat matrix"
% M = W^-1; % Variance-covariance matrix for the observations, assuming W = M^-1 right now ( https://en.wikipedia.org/wiki/Weighted_least_squares#Parameter_errors_and_correlation )
% Mr = (eye(n)-Hat)*M*(eye(n)-Hat)';
%!! This should be right but it's unused and slow so it's commented out
% Note that X'*W*residuals ~= 0 which means that W = M^-1 is true according to link above

%--- Standard Error --- !! May not be weighted right
% SE_r = stdevw_r/sqrt(n); % Estimate of the standard error ( https://en.wikipedia.org/wiki/Standard_error#Estimate ) <- w/ DoF instead of n very close to below
SE_r = sqrt(Sw_r/DoF); % SE of residuals from https://learnche.org/pid/least-squares-modelling/least-squares-model-analysis
% def wrong SE_r = sqrt(chiSqReduc); % Estimate of the standard error for regression ( https://en.wikipedia.org/wiki/Standard_error ) 
% SE_m = sqrt(RSSw_r/DoF*Mbeta(2,2)); % Estimate of standard error of slope ( https://stats.stackexchange.com/a/72062 )
% SE_b = sqrt(RSSw_r/DoF*Mbeta(1,1)); % Estimate of standard error of intercept ( https://stats.stackexchange.com/a/72062 )
% These seem to be wrong, according to matlab confint.m, SE_m should be stdevw_m and SE_b should be sstdevw_b
SE_m = stdevw_m; % Estimate of the standard error of slope
SE_b = stdevw_b; % Estimate of the standard error of intercept

%--- Sum of Squared [Stuff] ----
% Also known as residuals sum of squares (RSS), sum of squared residuals (SSR), or sum of squared estimate of errors (SSE)
SSE = RSSw_r; % Had it as RSSw_r AND Sw_r above, weighted of course!

%--- 95% Confidence interval ---
conf = .95; % 95% confidence
alpha = 1 - conf; % Based on https://www.mathworks.com/help/stats/tinv.html
crit_lo = tinv(alpha/2,DoF); % Lower bound is negative of upper
crit_up = tinv(1-alpha/2,DoF); % Upper one is the one to use for general b/c positive

% Based on confint.m (built into Matlab) it says:
% and the confidence intervals are
%    b +/- sqrt(diag(V)) * tinv(1-(1-level)/2, dfe)
% where V is the coefficient covariance matrix 'Mbeta' and dfe is 'DoF'
conf_m = [SE_m*crit_lo+m,SE_m*crit_up+m]; % Confidence interval of slope ( https://stats.stackexchange.com/a/72062 & confint.m )
conf_b = [SE_b*crit_lo+b,SE_b*crit_up+b]; % Confidence interval of intercept ( https://stats.stackexchange.com/a/72062 & confint.m )

% Load everything into the structure for returning
fitter = struct; % Prep the structure
fitter.m = m; % Record the slope
fitter.b = b; % Record the intercept
fitter.residuals = residuals; % Record the residuals
fitter.Jacobian = J; % Record the Jacobian
fitter.beta = beta; % Record the b/m vector for its order
fitter.covMat = Mbeta; % Record the Variance-Covariance matrix for the estimated parameters
fitter.variancem = variw_m; % Record the variance of the slope [part of Mbeta]
fitter.varianceb = variw_b; % Record the variance of the intercept [part of Mbeta]
fitter.variancer = chiSqReducw_r; % Record the estimated variance of the residuals [possibly chiSqReducw]
fitter.corrcoef = corrcoef; % Record the correlation coefficient [derived from Mbeta]
fitter.confm = conf_m; % Record 95% confidence interval for slope
fitter.confb = conf_b; % Record 95% confidence interval for intercept
fitter.variance = variw_r; % Record the unbiased weighted variance of the residuals
fitter.se = SE_r; % Record standard error of the residuals
fitter.rmse = RMSEw_r; % Record Root Mean Squared Error of residuals
fitter.sse = SSE; % Record Sum of Squared Estimate of Errors
fitter.dof = DoF; % Record the degrees of freedom

end