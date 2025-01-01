clc;
clear;

%% ----------import Data (Prec, Temp, Initial cond, lower and upper bound of parameter, and...)

global m_temp m_prec L_b U_b initial_c Temp Prec months stream_obs stream_obs_valid stream_obs_calib Prec_valid Prec_calib Temp_valid Temp_calib validation_id calibration_id yc Q_obs_ac

f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 1);
datee = data.Date;

stream_obs = data.Streamflow_m_3_s_;  % daily streamflow from 1/1/1981 to 12/03/2024 (m^3/s)
calibration_id = 1827:11322;  % From 1/1/1986 to 12/31/2011
validation_id = 11323:16043;  % From 1/1/2012 to 12/03/2024
stream_obs_valid = stream_obs(validation_id);
stream_obs_calib = stream_obs(calibration_id);
y = year(data.Date);  
yc = y(y > 1985 & y < 2012);
yv = y(y > 2011);
Q_obs_ac = groupsummary(stream_obs_calib, yc, 'mean');
Q_obs_av = groupsummary(stream_obs_valid, yv, 'mean');
months = month(datee);
Prec = data.Precipitation_mm_;  % daily precipitation from 1/1/1981 to 12/03/2024 (mm)
Prec_valid = Prec(validation_id);
Prec_calib = Prec(calibration_id);
Temp = data.Temperature_mean_; % daily temperature from 1/1/1981 to 12/03/2024   (C)
Temp_valid = Temp(validation_id);
Temp_calib = Temp(calibration_id);
%---------
f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 6);
m_temp = data.Temp(1:12);  % monthly mean of temperature
m_prec = data.Prec(1:12);  % monthly mean of precipitation
%---------
initial_c = data.InitialConditions(1:5);
%---------
L_b = data.LB(1:11); %lower bound
U_b = data.UB(1:11); %upper bound

%% User choice
choice = input('Enter 1 to see calibrated results or 2 to start calibration process: ');
switch choice
    case 1
        % Load and display calibrated results
        fvals = load('glue_fval.mat'); fvals = -fvals.glue_fval;                    % NSE values of acceptable combinations
        histories = load('glue_history.mat'); histories = histories.glue_history;   % History of each combination; from first sample to optimized values
        xsols = load('glue_xsol.mat'); xsols = xsols.glue_xsol;                     % Parameters combinations which simulate runoff with acceptable accuracy (NSE>=0.6)

    case 2
        % Start new calibration process
        x = input('Enter the number of iterations for the calibration process: ');
        [xsols, fvals, histories] = calibration_GLUE_fmincom(x);

        % Save the results
        % save('glue_xsol.mat', 'glue_xsol');
        % save('glue_fval.mat', 'glue_fval');
        % save('glue_history.mat', 'glue_history');

        fprintf('Calibration process completed and results saved.\n');

    otherwise
        fprintf('Invalid choice. Please enter 1 or 2.\n');
end

clearvars -except xsols fvals histories

%% --------------------------Start to Run the model with GLUE Uncertainty Analysis--------------------


function [glue_xsol,  glue_fval, glue_history] = calibration_GLUE_fmincom (N)

global L_b U_b
ub = U_b'; lb = L_b'; numDim = length(ub);


x0 = lhsdesign(N, numDim,'criterion','maximin', 'iterations',100 );
id=0;
for j = 1 : N
    x0(j, :) = x0(j, :) .* ( ub - lb) + lb;
    [xsol,fval,history] = calibrate_HBV(x0(j, :));

if -fval >= 0.6
    id=id+1;
    glue_fval(id) = fval;
    glue_history{id} = history;
    glue_xsol(id,:) = xsol;

end
end
end

%% ------------------------------------------Calibration Function--------------------------------------
function [xsol,fval,history,searchdir] = calibrate_HBV(x0)
global L_b U_b
ub = U_b'; lb = L_b'; 
% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];
 
% Call optimization

options = optimoptions('fmincon','OutputFcn',@outfun,'Display','iter', 'Algorithm','sqp');
[xsol,fval] = fmincon(@eval_HBV,x0,[],[],[],[],lb,ub,[],options);
 
 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         case 'done'
             hold off
         otherwise
     end
 end
end


%% -----------Calibration Error Calculation Function-----------

function [errfun]  = eval_HBV ( par_values )

global stream_obs_calib calibration_id yc Q_obs_ac

[ Q_sim ] = run_HBV( par_values);


metric = 'NSE'; % options: NSE, NSE_log, BIAS, ME, MAE, MSE, RMSE, ...

Q_obs = stream_obs_calib ;
Q_sim = Q_sim ( calibration_id );

switch metric
    case 'NSE'
        mean_obs = mean ( Q_obs );
        denominator = mean ( ( Q_obs - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( Q_obs - Q_sim ) .^ 2 );
        errfun = 1 - numerator / denominator;
       errfun = -1 * errfun; % used for optimization: to change problem to minimization
    case 'NSE_log'
        mean_obs = mean ( log(Q_obs) );
        denominator = mean ( ( log(Q_obs) - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( log(Q_obs) - log(Q_sim) ) .^ 2 );
        errfun = 1 - numerator / denominator;
        errfun = -1 * errfun; % used for optimization: to change problem to minimization
    case 'BIAS'
        errfun = sum ( (Q_sim - Q_obs) ) / sum ( Q_obs );
    case 'ME'
        errfun = mean ( Q_obs - Q_sim );
    case 'MAE'
        errfun = mean ( abs ( Q_obs - Q_sim ) );
    case 'MSE'
        errfun = mean ( ( Q_obs - Q_sim ) .^ 2 );
    case 'RMSE'
        errfun = sqrt ( mean ( ( Q_obs - Q_sim ) .^ 2 ) );
    case 'RSR'
        RMSE = sqrt ( mean ( ( Q_obs - Q_sim ) .^ 2 ) );
        SD_obs = std(Q_obs);
        errfun=  RMSE / SD_obs;

end
end
%% -----------Validation Error Calculation Function-----------

function errfun  = eval_HBV_valid ( Q_obser,par_values, ini_values_valid )

global stream_obs_valid validation_id yc Q_obs_ac

[ Q_sim , ET_sim , SMS_sim ] = run_HBV_valid( par_values,ini_values_valid);


metric = 'RMSE'; % options: NSE, NSE_log, BIAS, ME, MAE, MSE, RMSE, ...

Q_obs = Q_obser ;

switch metric
    case 'NSE'
        mean_obs = mean ( Q_obs );
        denominator = mean ( ( Q_obs - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( Q_obs - Q_sim ) .^ 2 );
        errfun = 1 - numerator / denominator;
       errfun = -1 * errfun; % used for optimization: to change problem to minimization
    case 'NSE_log'
        mean_obs = mean ( log(Q_obs) );
        denominator = mean ( ( log(Q_obs) - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( log(Q_obs) - log(Q_sim) ) .^ 2 );
        errfun = 1 - numerator / denominator;
%         errfun = -1 * errfun; % used for optimization: to change problem to minimization
    case 'BIAS'
        errfun = sum ( Q_sim - Q_obs ) / sum ( Q_obs );
    case 'ME'
        errfun = mean ( Q_obs - Q_sim );
    case 'MAE'
        errfun = mean ( abs ( Q_obs - Q_sim ) );
    case 'MSE'
        errfun = mean ( ( Q_obs - Q_sim ) .^ 2 );
    case 'RMSE'
        errfun = sqrt ( mean ( ( Q_obs - Q_sim ) .^ 2 ) );
end
end
%% Run_calibration

function [ Q_sim , ET_sim , SMS_sim , ini_values_valid, output_state] = run_HBV ( par_values )
% This function receives the parameter values of the HBV Hydrologic Model as input

global initial_c months Prec Temp m_prec m_temp calibration_id

watershed_area = initial_c(1); 
ini_values = [initial_c(1), initial_c(2), initial_c(3), initial_c(4)];
forcing = [Prec(1:calibration_id(end)), Temp(1:calibration_id(end)), months(1:calibration_id(end))];
long_term=[m_temp, m_prec];

[ output_flux , output_state ] = HBV_SASK ( forcing , long_term , par_values , ini_values , watershed_area );
Q_sim = output_flux ( 1 : end, 1 );
ET_sim = output_flux ( 1 : end, 7 );
SMS_sim = output_state( 1 : end, 2 );
ini_values_valid= output_state(end,:) ;
end
%% Run_validation

function [ Q_sim , ET_sim ] = run_HBV_valid ( par_values, ini_values_valid )
% This function receives the parameter values of the HBV Hydrologic Model as input

global initial_c months Prec Temp m_prec m_temp validation_id

watershed_area = initial_c(1); 
ini_values = ini_values_valid;
forcing = [Prec(validation_id), Temp(validation_id), months(validation_id)];
long_term=[m_temp, m_prec];

[ output_flux] = HBV_SASK ( forcing , long_term , par_values , ini_values , watershed_area );
Q_sim = output_flux ( 1 : end, 1 );
ET_sim = output_flux ( 1 : end, 7 );
end
%% --------------------------------HBV Model--------------------------

function [ output_flux , output_state ] = HBV_SASK ( forcing , long_term , par_values , ini_values , watershed_area )

TT = par_values(1); C0 = par_values(2); ETF = par_values(3);
LP = par_values(4); FC = par_values(5); beta = par_values(6);
FRAC = par_values(7); K1 = par_values(8); alpha = par_values(9); 
K2 = par_values(10);

if length(par_values) >= 11; UBAS = par_values(11);
else UBAS = 1;
end
if length(par_values) == 12; PM = par_values(12); 
else PM = 1;
end

LP = LP * FC;

initial_SWE = ini_values(1); initial_SMS = ini_values(2); 
initial_S1= ini_values(3);   initial_S2 = ini_values(4);

P = PM * forcing(:,1); T = forcing(:,2);
month_time_series = forcing(:,3);
monthly_average_T = long_term(:,1);
monthly_average_PE = long_term(:,2);

SWE(1) = initial_SWE; SMS(1) = initial_SMS;
S1(1) = initial_S1; S2(1) = initial_S2;
period_length = size(P , 1);

for t = 1 : period_length

    [ SWE(t+1), ponding(t) ] = precipitation_module ...
                                            ( SWE(t), P(t), T(t), TT, C0);
                                        
    [ AET(t) , PET(t) ] = evapotranspiration_module ...
                                            ( SMS(t), T(t), ...
                                            month_time_series(t) , ...
                                            monthly_average_T, ...
                                            monthly_average_PE, ETF, LP);
                                        
    [ SMS(t+1), S1(t+1), S2(t+1), Q1(t), Q2(t) ] = ...
                                            soil_storage_routing_module ...
                                            ( ponding(t), SMS(t), ...
                                            S1(t), S2(t), AET(t), ...
                                            FC, beta, FRAC, K1, alpha, K2);
end

[Q1_routed] = triangle_routing(Q1, UBAS);
Q = Q1_routed + Q2;
Q_cms = ( Q * watershed_area * 1000 ) /( 24 * 3600 );


output_flux = [ Q_cms', Q', Q1', Q2', ponding', PET', AET' ];
output_state = [ SWE', SMS', S1', S2'  ];
end

%% -------------------------Module Functions (precipitation, evapotranspiration and ....)-----

function [Q_routed] = triangle_routing(Q, UBAS)
UBAS = max(UBAS, 0.1);
length_triangle_base = ceil(UBAS);
if UBAS == length_triangle_base
    x = [0, 0.5 * UBAS, length_triangle_base];
    v = [0,      1      ,            0        ];
else
    x = [0, 0.5 * UBAS, UBAS, length_triangle_base];
    v = [0,      1      ,   0   ,            0        ];
end
for i = 1 : length_triangle_base
    if (i-1) < (0.5 * UBAS) && i > (0.5 * UBAS)
        weight(i) = 0.5 * ( interp1(x,v, i - 1) + interp1(x,v, 0.5 * UBAS) ) * ( 0.5 * UBAS - i + 1) +  0.5 * ( interp1(x,v, 0.5 * UBAS) + interp1(x,v, i) ) * ( i - 0.5 * UBAS ) ;
    elseif i > UBAS
        weight(i) = 0.5 * ( interp1(x,v, i - 1) ) * ( UBAS - i + 1);
    else
        weight(i) = interp1(x,v, i - 0.5);
    end
end
weight = weight / sum(weight);

for i = 1 : length(Q)
    temp = 0;
    for j = 1 : min( i, length_triangle_base )
        temp = temp + weight(j) * Q(i - j + 1);
    end
    Q_routed(i) = temp;
end
end

function [ SWE_new, ponding ] = precipitation_module ( SWE, P, T, TT, C0)
% *****  TT : Temperature Threshold or melting/freezing point - model parameter *****
% *****  C0: base melt factor - model parameter ***** 
% *****  P: Precipitation - model forcing ***** 
% *****  T: Temperature - model forcing ***** 
% *****  SWE: Snow Water Equivalent - model state variable ***** 

if T >= TT
    rainfall = P;
    potential_snow_melt  = C0 * (T - TT);
    snow_melt = min ( potential_snow_melt  , SWE );
    ponding = rainfall + snow_melt; % Liquid Water on Surface
    SWE_new = SWE - snow_melt; % Soil Water Equivalent - Solid Water on Surface
else
    snowfall = P;
    snow_melt = 0;
    ponding = 0; % Liquid Water on Surface
    SWE_new = SWE + snowfall;  % Soil Water Equivalent - Solid Water on Surface
end
end

function [ AET , PET ] = evapotranspiration_module ( SMS, T, ...
                            month_number, monthly_average_T, ...
                            monthly_average_PE, ETF, LP )
% *****  T: Temperature - model forcing *****
% *****  month_number: the current month number - for Jan=1, ..., Dec=12 *****
% *****  SMS: Soil Moisture Storage - model state variable *****
% *****  ETF - This is the temperature anomaly correction of potential evapotranspiration - model parameters
% *****  LP: This is the soil moisture content below which evaporation becomes supply-limited - model parameter
% *****  PET: Potential EvapoTranspiration - model parameter
% *****  AET: Actual EvapoTranspiration - model 
PET = ( 1 + ETF * ( T - monthly_average_T(month_number) ) ) * monthly_average_PE(month_number);  % Potential Evapotranspiration
PET = max(PET, 0);
PET = min(PET, 2*monthly_average_PE(month_number));

if SMS > LP
    AET = PET;
else
    AET = PET * ( SMS / LP );
end
AET = min (AET, SMS); % to avoid evaporating more than water available
end

function [ SMS_new, S1_new, S2_new, Q1, Q2 ] = soil_storage_routing_module ...
                ( ponding, SMS, S1, S2, AET, FC, beta, FRAC, K1, alpha, K2)


% *****  T: Temperature - model forcing *****
% *****  month_number: the current month number - for Jan=1, ..., Dec=12 *****
% *****  SMS: Soil Moisture Storage - model state variable *****
% *****  ETF - This is the temperature anomaly correction of potential evapotranspiration - model parameters
% *****  LP: This is the soil moisture content below which evaporation becomes supply-limited - model parameter
% *****  PET: Potential EvapoTranspiration - model parameter


% *****  FC: Field Capacity - model parameter ---------
% *****  beta: ????? - model parameter ---------
% This controls the relationship between soil infiltration and soil water release.
% The default value is 1. Values less than this indicate a delayed response, while higher
% values indicate that runoff will exceed infiltration.

if SMS < FC
    soil_release = ponding * (( SMS / FC )^beta); % release of water from soil
else
    soil_release = ponding; % release of water from soil
end

SMS_new = SMS - AET + ponding - soil_release;
soil_release_to_fast_reservoir = FRAC * soil_release;
soil_release_to_slow_reservoir = ( 1 - FRAC ) * soil_release;

Q1 = K1 * S1 ^ alpha;
if Q1 > S1
    Q1 = S1;
end

S1_new = S1 + soil_release_to_fast_reservoir - Q1;

Q2 = K2 * S2;

S2_new = S2 + soil_release_to_slow_reservoir - Q2;
end

function [NSE, NSE_log, PBIAS, RMSE, CC, KGE_d, RSR]  = error_model ( Q_obs, Q_sim )


% 
% Q_sim = Q_sim ( calibration_id );
% 
% Q_sim_ac = groupsummary(Q_sim, yc, 'mean');
mu_obs = mean(Q_obs);
mu_sim = mean(Q_sim);
sigma_obs = std(Q_obs);
sigma_sim = std(Q_sim);
r_d = corr(Q_obs, Q_sim);
beta_d = mu_sim / mu_obs; % Bias ratio
gamma_d = sigma_sim / sigma_obs; % Variability ratio

KGE_d = 1 - sqrt((r_d - 1)^2 + (beta_d - 1)^2 + (gamma_d - 1)^2);
% 
% mu_obs = mean(Q_obs_ac);
% mu_sim = mean(Q_sim_ac);
% sigma_obs = std(Q_obs_ac);
% sigma_sim = std(Q_sim_ac);
% r_a = corr(Q_obs_ac, Q_sim_ac);
% beta_a = mu_sim / mu_obs; % Bias ratio
% gamma_a = sigma_sim / sigma_obs; % Variability ratio
% KGE_a = 1 - sqrt((r_a - 1)^2 + (beta_a - 1)^2 + (gamma_a - 1)^2);
% 
% errfun = sqrt((1 - KGE_d)^2 + (1 - KGE_a)^2);

        mean_obs = mean ( Q_obs );
        denominator = mean ( ( Q_obs - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( Q_obs - Q_sim ) .^ 2 );
        NSE = 1 - numerator / denominator;
        % NSE = -1 * errfun; % used for optimization: to change problem to minimization

        mean_obs = mean ( log(Q_obs) );
        denominator = mean ( ( log(Q_obs) - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
        numerator = mean ( ( log(Q_obs) - log(Q_sim) ) .^ 2 );
        NSE_log = 1 - numerator / denominator;
        % NSE_log = -1 * errfun; % used for optimization: to change problem to minimization

        PBIAS = sum ( abs(Q_sim - Q_obs) ) / sum ( Q_obs )*100;

        RMSE = sqrt ( mean ( ( Q_obs - Q_sim ) .^ 2 ) );
        SD_obs = std(Q_obs);
        RSR=  RMSE / SD_obs;
        CC=(corr(Q_obs,Q_sim)).^2;
end

