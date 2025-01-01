clc;
clear;

%% Import Data

f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 2);
datee= data.Date;                        % Date
m_stream = data.Streamflow_m_3_s_ ;           % Monthly Streamflow
m_rainfall = data.Precipitation_mm_ ;         % Monthly Precipitation


%% Dataset of Different Temporal Scale (e.g., Jan, Feb,... Winter, Spring,..., Annual)

Months = month(datee); un_Months = unique(month(datee));
Years = year(datee); un_Years = unique(year(datee));

% each month
for i = 1:length(un_Months)
    id = find(Months==un_Months(i));
    months_rainfall{i} = m_rainfall(id);
    months_stream{i} = m_stream(id);
end

% Seasonal Scale

winter_id = [1 2 3];
spring_id = [4 5 6];
summer_id = [7 8 9];
autumn_id = [10 11 12];

for i = 1:length(un_Years)
    id = find(Years==un_Years(i));
    x=m_rainfall(id);
    y=m_stream(id);
    seasons_rainfall(i,1) = sum(x(1:3));     % Winter
    seasons_rainfall(i,2) =sum(x(4:6));      % Spring
    seasons_rainfall(i,3) = sum(x(7:9));     % Summer
    seasons_stream(i,1) = sum(y(1:3));       % Winter
    seasons_stream(i,2) =sum (y(4:6));       % Spring
    seasons_stream(i,3) = sum(y(7:9));       % Summer
    if i~=length(un_Years)
    seasons_rainfall(i,4) = sum(x(10:12));   % Autumn
    seasons_stream(i,4) = sum(y(10:12));     % Autumn
    else
    seasons_rainfall(i,4) = sum(x(10:11));   % Autumn
    seasons_stream(i,4) = sum(y(10:11));     % Autumn
    end
end

% Annual scale

for i = 1:length(un_Years)
    id = find(Years == un_Years(i));
    annual_rainfall(i) = sum(m_rainfall(id));
    annual_stream(i) = sum(m_stream(id));
end

for i=1:12
    sum_rainfall(i) = sum(months_rainfall{i});
    mean_stream(i) = mean(months_stream{i});
end
%------------- Modified Mannlendall_test (Trend Analysis)------------------

% Trends of months

for i = 1:size(months_stream,2)
    in = months_stream{i};
    in2 = months_rainfall{i};
    [tau, z_score, p_value, H] = Modified_MannKendall_test(1:length(in), in, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in); in']', 0.05);
    Trends_streamflow_months (i,:) = [tau, z_score, p_value, H, beta, HH]; 
    [tau2, z_score2, p_value2, H2] = Modified_MannKendall_test(1:length(in2), in2, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in2); in2']', 0.05);
     Trends_rainfall_months (i,:) = [tau2, z_score2, p_value2, H2, beta, HH];   
end

% Seosonal Trend

for i=1:size(seasons_stream,2)
    in = seasons_stream(:,i);
    in2 = seasons_rainfall(:,i);
    [tau, z_score, p_value, H] = Modified_MannKendall_test(1:length(in), in, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in); in']', 0.05);
    Trends_streamflow_seasons (i,:) = [tau, z_score, p_value, H, beta, HH];
    [tau2, z_score2, p_value2, H2] = Modified_MannKendall_test(1:length(in2), in2, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in2); in2']', 0.05);
    Trends_rainfall_seasons  (i,:) = [tau2, z_score2, p_value2, H2, beta, HH];      
end

% Annual Trend

    in= annual_stream;
    in2= annual_rainfall;
    [tau, z_score, p_value, H] = Modified_MannKendall_test(1:length(in), in, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in); in]', 0.05);
    Trends_streamflow_annual = [tau, z_score, p_value, H, beta, HH];
    [tau2, z_score2, p_value2, H2] = Modified_MannKendall_test(1:length(in2), in2, 0.05, 0.05);
    [HH, beta] = Sen_Slope([1:length(in2); in2]', 0.05);
    Trends_rainfall_annual = [tau2, z_score2, p_value2, H2, beta, HH]; 

%------------------------ Statistical Summary--------------------------

for i = 1:17
    if i<= 12
[Max Min Mean Std Kurt Skew CV] = statistical_summary(months_stream{i});
        summary_discription_streamflow(i,:)=[Max Min Mean Std Kurt Skew CV];
[Max Min Mean Std Kurt Skew CV] = statistical_summary(months_rainfall{i});
        summary_discription_rainfall(i,:) = [Max Min Mean Std Kurt Skew CV];
    elseif i>12 && i<=16
        [Max Min Mean Std Kurt Skew CV] = statistical_summary(seasons_stream(:,i-12));
        summary_discription_streamflow(i,:)=[Max Min Mean Std Kurt Skew CV];
        [Max Min Mean Std Kurt Skew CV] = statistical_summary(seasons_rainfall(:,i-12));
        summary_discription_rainfall(i,:)=[Max Min Mean Std Kurt Skew CV];
    else
        [Max Min Mean Std Kurt Skew CV]= statistical_summary(annual_stream);
        summary_discription_streamflow(i,:)=[Max Min Mean Std Kurt Skew CV];
        [Max Min Mean Std Kurt Skew CV]= statistical_summary(annual_rainfall);
        summary_discription_rainfall(i,:) = [Max Min Mean Std Kurt Skew CV];  
    end
end

%--------------------Correlation------------------------------------

for i=1:17
    if i<= 12
        [rho,pval]= corr(months_stream{i},months_rainfall{i},'Type','Pearson');
        Corrs_pearson(i,:) = [rho,pval];
        [rho,pval] = corr(months_stream{i},months_rainfall{i},'Type','Kendall');
        Corrs_kendal(i,:) = [rho,pval];
        [rho,pval] = corr(months_stream{i},months_rainfall{i},'Type','Spearman');
        Corrs_spearman(i,:) = [rho,pval];

    elseif i>12 && i<=16
        [rho,pval] = corr(seasons_stream(:,i-12),seasons_rainfall(:,i-12),'Type','Pearson');
        Corrs_pearson(i,:) = [rho,pval];
        [rho,pval] = corr(seasons_stream(:,i-12),seasons_rainfall(:,i-12),'Type','Kendall');
        Corrs_kendal(i,:) = [rho,pval];
        [rho,pval] = corr(seasons_stream(:,i-12),seasons_rainfall(:,i-12),'Type','Spearman');
        Corrs_spearman(i,:) = [rho,pval];
    else
        [rho,pval]= corr(annual_stream' , annual_rainfall','Type','Pearson');
        Corrs_pearson(i,:) = [rho,pval];
        [rho,pval] = corr(annual_stream' , annual_rainfall','Type','Kendall');
        Corrs_kendal(i,:) = [rho,pval];
        [rho,pval] = corr(annual_stream' , annual_rainfall','Type','Spearman');
        Corrs_spearman(i,:) = [rho,pval];
    end
end

clearvars -except Corrs_pearson  Corrs_kendal Corrs_spearman summary_discription_streamflow ...
    summary_discription_rainfall Trends_streamflow_annual Trends_rainfall_annual ... 
    Trends_streamflow_seasons Trends_rainfall_seasons Trends_streamflow_months Trends_rainfall_months

%% Return Level Quantification of extreme events

f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 1);
data.Date = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
streamflow = data.Streamflow_m_3_s_ ;
precip=data.Precipitation_mm_;

% Adjust to Water Year (October 1 - September 30)
data.WaterYear = year(data.Date); % Initialize WaterYear
isOctToDec = month(data.Date) >= 10; % Identify dates from October to December
data.WaterYear(isOctToDec) = data.WaterYear(isOctToDec) + 1; % Assign the next year for water year
years=unique((data.WaterYear));

% Get daily streamflow for the current year
for i = 1:length(years)
    s = streamflow(data.WaterYear == years(i));
    r = precip(data.WaterYear == years(i));
    % Find the maximum flow for the year
    annualMaxFlow(i) = max(s);
    annualMaxRainfall(i) = max(r);
end

%-------------------Maximum Annual Precipitation and streamflow-----------

%% maximum correlation analysis
figure;
plot(years,annualMaxFlow);
hold on
plot(years,annualMaxRainfall);
xlabel('Annual Maximum Rainfall (mm)');
ylabel('Annual Maximum Flow Discharge (m^3/s)');
title('Correlation Between Extreme Rainfall and Flow Discharge');
grid on;
axis square
% Plot Relationship Between Extreme Rainfall and Flow
figure;
scatter(annualMaxRainfall, annualMaxFlow, 'filled');
xlabel('Annual Maximum Rainfall (mm)');
ylabel('Annual Maximum Flow Discharge (m^3/s)');
title('Correlation Between Extreme Rainfall and Flow Discharge');
grid on;
axis square
[rho,pval] = corr(annualMaxFlow',annualMaxRainfall','Type','Pearson');
Corrs_extreme(1,:) = [rho,pval];
[rho,pval] = corr(annualMaxFlow',annualMaxRainfall','Type','Kendall');
Corrs_extreme(2,:) = [rho,pval];
[rho,pval] = corr(annualMaxFlow',annualMaxRainfall','Type','Spearman');
Corrs_extreme(3,:) = [rho,pval];

% ------------------lag time calculation------------------------


precip_norm = (precip - mean(precip)) / std(precip);
streamflow_norm = (streamflow - mean(streamflow)) / std(streamflow);
% Compute cross-correlation
[maxLag] = 100; 
[ccf, lags] = xcorr(streamflow_norm, precip_norm, maxLag, 'coeff');
% Find the lag with the maximum cross-correlation
[~, maxIdx] = max(ccf);
optimalLag = lags(maxIdx); 


figure;
stem(lags, ccf, 'filled');
xlabel('Lag (timesteps)');
ylabel('Cross-correlation coefficient');
title('Cross-correlation between Precipitation and Streamflow');
grid on;

% Display result
fprintf('The general lag time between extreme precipitation and streamflow is %d timesteps.\n', optimalLag);

%------------------GEV Distribution Fitting --------------
gevParams_flow = fitdist(annualMaxFlow', 'gev');
gevParams_rainfall = fitdist(annualMaxRainfall', 'gev');

figure;

subplot(1,2,1)
histogram(annualMaxFlow, 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on
axis square
% Generate x values for GEV PDF
x_values = linspace(min(annualMaxFlow) - 10, max(annualMaxFlow) + 10, 1000);
gevPDF_flow = pdf(gevParams_flow, x_values);
plot(x_values, gevPDF_flow, 'r-', 'LineWidth', 2);
xlabel('Annual Maximum Streamflow');
ylabel('Probability Density');
title('Fitted GEV Distribution with Histogram');
legend('Histogram', 'Fitted GEV PDF');
grid on;
disp('Fitted GEV Parameters:');
disp(gevParams_flow);

subplot(1,2,2)
histogram(annualMaxRainfall, 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on
axis square
% Generate x values for GEV PDF
x_values = linspace(min(annualMaxRainfall) - 10, max(annualMaxRainfall) + 10, 1000);
gevPDF_rainfall = pdf(gevParams_rainfall, x_values);
plot(x_values, gevPDF_rainfall, 'r-', 'LineWidth', 2);
xlabel('Annual Maximum Rainfall');
ylabel('Probability Density');
title('Fitted GEV Distribution with Histogram');
legend('Histogram', 'Fitted GEV PDF');
grid on;
disp('Fitted GEV Parameters:');
disp(gevParams_rainfall);

%%-----------------Fitted GEV Validation--------------------------------

figure;
subplot(1,2,1)
qqplot(annualMaxFlow, gevParams_flow);
title('Q-Q Plot for Streamflow GEV Fit');

% KS Test for GEV-Streamflow
[h, p_value_gev_flow] = kstest(annualMaxFlow, 'CDF', gevParams_flow);
disp(['KS Test for GEV-Streamflow: p-value = ', num2str(p_value_gev_flow)]);
axis square

subplot(1,2,2)
qqplot(annualMaxRainfall, gevParams_rainfall);
title('Q-Q Plot for Rainfall GEV Fit');
axis square
% KS Test for GEV-Streamflow
[h, p_value_gev_rainfall] = kstest(annualMaxRainfall, 'CDF', gevParams_rainfall);
disp(['KS Test for GEV-Rainfall: p-value = ', num2str(p_value_gev_rainfall)]);
%------------------Return level calculation using GEV method--------------

T = [2, 5, 10, 50, 100]; % Return periods in years

% Return Levels using GEV
Return_GEV_flow = gevinv(1 - 1 ./ T, gevParams_flow.k, gevParams_flow.sigma, gevParams_flow.mu);
Return_GEV_rainfall = gevinv(1 - 1 ./ T, gevParams_rainfall.k, gevParams_rainfall.sigma, gevParams_rainfall.mu);


figure;
subplot(1,2,1)
plot(T, Return_GEV_flow, 'o-', 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.8]);
xlabel('Return Period (years)');
ylabel('Return Level (m^3/s)');
title('Return Levels - Streamflow');
set(gca, 'XScale', 'log'); % Set X-axis to logarithmic scale
set(gca, 'XTick', T, 'XTickLabel', string(T)); % Set custom X-ticks
set(gca, 'YTick', Return_GEV_flow, 'YTickLabel', string(Return_GEV_flow)); % Set custom X-ticks
grid on;
axis square
subplot(1,2,2)
plot(T, Return_GEV_rainfall, 'o-', 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.8]);
xlabel('Return Period (years)');
ylabel('Return Level (mm)');
title('Return Levels - Rainfall');
set(gca, 'XScale', 'log'); % Set X-axis to logarithmic scale
set(gca, 'XTick', T, 'XTickLabel', string(T)); % Set custom X-ticks
set(gca, 'YTick', Return_GEV_rainfall, 'YTickLabel', string(Return_GEV_rainfall)); % Set custom X-ticks
grid on;
axis square




%% -----------------Functions--------------------

% Modified Mannlendall_test
function [tau, z, p, H] = Modified_MannKendall_test(t, X, alpha, alpha_ac)
    
    %% FUNCTION INPUTS AND OUTPUTS
    
    % INPUTS
    % t         - time vector corresponding to the timeseries.
    % X         - The timeseries for which kendall's tau value is to be evaluated.
    % alpha     - significance level above which kendall tau values are statistically significant. One-tailed test.
    % alpha_ac  - predetermined level for selecting only those autocorrelation lag values that are statistically significant. Two-tailed test.

    % OUTPUTS
    % tau   - The kendall rank correlation coefficient for the timeseries.
    % z     - z-score of the obtained tau value.
    % p     - p-value of the obtained tau value.
    % H     - Denotes whether to reject the null hypothesis or not. Null Hypothesis: There is no trend in the data. 0 -> retain the null hypothesis, 1 -> reject and increasing trend, -1 -> reject and decreasing trend, 2 -> if variance turns out to be negative.


    %% REFERENCES

    % 1) Kendall, M. G. (1938). A new measure of rank correlation. Biometrika, 30(1/2), 81-93.
    % 2) Kendall, M. G. (1948). Rank correlation methods.
    % 3) Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of hydrology, 204(1-4), 182-196.
    % 4) Hamed, K. H. (2008). Trend detection in hydrologic data: the Mann–Kendall trend test under the scaling hypothesis. Journal of hydrology, 349(3-4), 350-363.
    % 5) Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water resources management, 18(3), 201-218.


    %% CONVERT DATA INTO ROW VECTOR
    tic
    % Convert input time and timeseries vectors to row vectors
    X = reshape(X, 1, length(X));
    t = reshape(t, 1, length(t));
    n = length(X);

    a = toc;


    %% CALCULATE THE KENDALL TAU VALUE
    tic
    % Notation taken from: 3) Hamed, K. H., & Rao, A. R. (1998)
    S = 0;
    
    for i = 1: n-1
        for j = i+1: n
            S = S + sign(X(j) - X(i));
        end
    end
    
    % Calculate kendall's rank correlation coefficient - tau
    tau = S / (n * (n-1) / 2);

    b = toc;


    %% CALCULATE VARIANCE OF KENDALL TAU UNCORRECTED FOR AUTOCORRELATION
    tic
    % Correction for tied ranks taken from:
    % 2) Kendall, M. G. (1948) edition 5, chapter 4, page 66
    % 4) Hamed, K. H. (2008)

    % Variance of kendall tau when the following conditions hold: 1) No autocorrelation, 2) No ties
    var_S_notie = n * (n - 1) * (2*n + 5) / 18;
    
    % Correcting variance in the case of tied ranks    
    tied_ranks = [];
    Y = sort(X);
    % Add dummy value to Y. Otherwise loop doesn't count tied ranks that exist at the end of Y.
    Y = [Y, Y(end) + 1];
    prev_counter = 1;
    current_counter = 1;
    for i_Y = 2: n+1
        if Y(i_Y) == Y(i_Y - 1)
            current_counter = current_counter + 1;
        else
            prev_counter = current_counter;
            current_counter = 1;
        end

        if current_counter == 1 && prev_counter ~= 1
            tied_ranks = [tied_ranks, prev_counter];
        end

    end

    var_S_tie_correction = sum(tied_ranks .* (tied_ranks - 1) .* (2*tied_ranks + 5) / 18);

    % Calculate variance corrected for ties but without considering autocorrelation
    var_S_noAC = var_S_notie - var_S_tie_correction;

    c = toc;

    
    %% REMOVE SEN TREND ESTIMATE FROM THE DATA
    tic
    % Procedure and rationale for trend removal taken from:
    % 5) Yue, S., & Wang, C. (2004)
    % Wikipedia: Theil–Sen estimator, https://en.wikipedia.org/wiki/Theil-Sen_estimator
    % Basically, the presence of a trend leads to wrong estimation of the actual autocorrelation present. Therefore, the trend must first be removed before estimating the autocorrelation.
    m_list = zeros(1, (n * (n-1) / 2));
    b_list = [];
    
    element_counter = 1;
    for i = 1: n-1
        for j = i+1: n
            m_list(element_counter) = (X(j) - X(i)) / (t(j) - t(i));
            element_counter = element_counter + 1;
        end
    end
    
    m_sen = median(m_list);

    b_list = X - m_sen*t;
    b_sen = median(b_list);
    
    % Remove sen trend estimate from the data
    X = X - m_sen*t - b_sen;

    d = toc;


    %% CALCULATE AUTOCORRELATION VALUES FOR STATISTICALLY SIGNIFICANT LAGS
    tic
    X_rank_order = tiedrank(X);
    z_ac = abs(norminv(alpha_ac / 2));  % norminv() is the inverse of the normcdf() function
    [acf, ~, acf_bounds] = autocorr(X_rank_order, NumLags = n - 1, NumSTD = z_ac);
    
    % Retain only those lags for which the autocorrelation value is statistically significant
    rho_lags = [];
    rho = [];

    for i = 2: n
        if acf(i) > acf_bounds(1) || acf(i) < acf_bounds(2)
            rho = [rho, acf(i)];
            rho_lags = [rho_lags, i-1];
        end
    end

    e = toc;
    

    %% CALCULATE AUTOCORRELATION CORRECTED VARIANCE OF KENDALL TAU
    tic
    % Calculate variance correction factor
    const_factor = 2 / (n * (n-1) * (n-2));
    rho_factor_sum = 0;
    for i = 1: length(rho)
        rho_factor_sum = rho_factor_sum + (n - rho_lags(i)) * (n - rho_lags(i) - 1) * (n - rho_lags(i) - 2) * rho(i);
    end

    var_AC_correction_factor = 1 + const_factor * rho_factor_sum;

    % Calculate variance corrected for autocorrelation
    var_S = var_S_noAC * (var_AC_correction_factor);
    
    f = toc;


    %% CHECK FOR STATISTICAL SIGNIFICANCE OF THE KENDALL TAU VALUE
    tic
    % Since the correction factor for the true variance is an approximation, in rare cases it may turn out to be negative. In that scenario abort the function and return H = 2 as an exception value.
    if var_S < 0
        z = 0;
        p = 0.5;
        H = 2;
        return
    end

    % Calculate z-score.
    % z-score = (value - mean) / std_dev. That is, how far is the value from mean as a multiple of the standard deviation.
    % Below the +1 and -1 are for "continuity correction".
    % Continuity correction taken from: 2) Kendall, M. G. (1948), edition 5, chapter 4, page 65
    if S > 0
        z = (S - 1) / sqrt(var_S);
    elseif S == 0
        z = 0;
    elseif S < 0
        z = (S + 1) / sqrt(var_S);
    end

    % Calculate p-value.
    % p-value is the probability of obtaining values further from the obtained value, on either the negative or positive tail.
    % p-value found here is for one-tailed test. Right tail values denote positive trend probabilities and left tail values represent negative trend probabilities
    if z >= 0
        p = 1 - normcdf(z);   % normcdf() returns the cdf value from the standard normal.
    elseif z < 0
        p = normcdf(z);
    end
    
    % Set H0 (Null Hypothesis) rejection value.
    % H = 0 -> retained. H = 1 -> rejected and increasing. H = -1 -> rejected and decreasing.
    if p <= alpha
        % That is there is a trend
        if z > 0
            % Positive trend
            H = 1;
        elseif z < 0
            % Negative trend
            H = -1;
        end
    else
        % That is there is no trend
        H = 0;
    end

    g = toc;

    MMK_total_time = a + b + c + d + e + f;
    
%     fprintf('Unoptimized Mann-Kendall Function\n');
%     fprintf("a = %f\n", a);
%     fprintf("b = %f\n", b);
%     fprintf("c = %f\n", c);
%     fprintf("d = %f\n", d);
%     fprintf("e = %f\n", e);
%     fprintf("f = %f\n", f);
%     fprintf("g = %f\n", g);
%     fprintf("MMK_total_time = %f\n", MMK_total_time);

    % time_taken = [time_taken; a, b, c, d, e, f];


end
% Sen Slope estimator
function [H, beta1, LL, UL, beta0, H05, H01] = Sen_Slope(TSD, alpha)
%% Calculate slopes and beta1
[~, m]=size(TSD);
if m>2; error('The input time serise is defined for one variable; TSD = (n x 2) double'); end
NaNData=sum(sum(isnan(TSD)));if NaNData>0; display(NaNData); disp('NaN(s) data is found'); error('d'); end
%%
y = TSD(:,2);
n = length(y);
N = n*(n-1)/2; % Number of slope estimates
Q = zeros(1, N);
c=0;
for j=2:n
    for i=1:j-1
        c=c+1;
        Q(c)=(y(j)-y(i))/(j-i);
    end
end
Q = sort(Q);
beta1 = median(Q);
%% Nonparametric Intercept, beta0, of the Linear Trend
y_MED = median(y);
t_MED = median(TSD(:,1));
beta0 = y_MED-beta1*t_MED;
%% Approximate Confidence Limits for the Nonparametric Slope Estimate beta1
% This method assumes that n is large, say n>10.
if n<=10; warning('The upper and lower limits on the estimated slope are approximate when n<10'); end
Z_alpha_p_2 = norminv(1-alpha/2, 0, 1);
Var_S = n*(n-1)*(2*n+5)/18;
C_a = Z_alpha_p_2*sqrt(Var_S);
L = round((N-C_a)/2);
LL = Q(L); % Lower limits for the slope
U = round((N+C_a)/2+1);
UL = Q(U); % Upper limits for the slope
%% Test hypothesis H0 for the beta1/slope at the alpha significance level
if LL<=0 && UL>=0; H=0; end
if LL<0 && UL<0;   H=1; end
if LL>0 && UL>0;   H=1; end
%% Test hypothesis H0 for the beta1/slope at the 0.05 and 0.01 (alpha = 0.05 & 0.01) significance level
Z_05_p_2 = norminv(1-.05/2, 0, 1);
Z_01_p_2 = norminv(1-.01/2, 0, 1);
C_a_05 = Z_05_p_2*sqrt(Var_S);
C_a_01 = Z_01_p_2*sqrt(Var_S);
L_05 = round((N-C_a_05)/2);
LL_05 = Q(L_05); % Lower limits for the slope
U_05 = round((N+C_a_05)/2+1);
UL_05 = Q(U_05); % Upper limits for the slope
L_01 = round((N-C_a_01)/2);
LL_01 = Q(L_01); % Lower limits for the slope
U_01 = round((N+C_a_01)/2+1);
UL_01 = Q(U_01); % Upper limits for the slope
if LL_05<=0 && UL_05>=0; H05=0; StarH=0;   %display('Trend is NOT significant  at the 95% confidence level');
end
if LL_05<0 && UL_05<0;   H05=1; StarH='$^{(*)}$'; %display('Trend is significant  at the 5% significance level; Values tend to decrease with time');
end
if LL_05>0 && UL_05>0;   H05=1; StarH='$^{(*)}$'; %display('Trend is significant  at the 5% significance level; Values tend to increase with time');
end
if LL_01<=0 && UL_01>=0; H01=0;   %display('Trend is NOT significant  at the 99% confidence level');
end
if LL_01<0 && UL_01<0; H01=1; StarH='$^{(**)}$'; %display('Trend is significant  at the 1% significance level; Values tend to decrease with time');
end
if LL_01>0 && UL_01>0; H01=1; StarH='$^{(**)}$'; %display('Trend is significant  at the 1% significance level; Values tend to increase with time');
end

end
% Statistical Descriprion
function [Max, Min, Mean, Std, Kurt, Skew, CV]= statistical_summary(input)
Max=max(input);
Min=min(input);
Mean=mean(input);
Std=std(input);
Kurt=kurtosis(input);
Skew=skewness(input);
CV=Std/Mean*100;
end