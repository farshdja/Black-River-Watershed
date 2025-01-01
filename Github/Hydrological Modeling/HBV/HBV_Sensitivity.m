clc;
clear;
global d_stream d_rainfall d_temp initial_c m_temp m_prec months

%--------
f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 1);
datee= data.Date;                             % Date
d_stream = data.Streamflow_m_3_s_ ;           % daily Streamflow
d_rainfall = data.Precipitation_mm_ ;         % daily Precipitation
d_temp=data.Temperature_mean_;
date_hist=data.Date; months=month(date_hist);

%--------
f = 'Final Dataset.xlsx';
data = readtable(f, 'Sheet', 6);
m_temp = data.Temp(1:12);  % monthly mean of temperature
m_prec = data.Prec(1:12);  % monthly mean of precipitation
%---------
initial_c = data.InitialConditions(1:5);
%---------
lb = data.LB(1:11); %lower bound
ub = data.UB(1:11); %upper bound
%---------
T_N=5;

for n=1:length(T_N)                           
 N=T_N(n);

 funPath = 'C:\Users\Lenovo\Desktop\Dr. Hassanzadeh\codes' ; 
 funFile = 'hbv_sobol';            
 smplMtd = 'LHS';                     
 seedNum=[];
 bootstrapFlag = 0;                    
 bootstrapSize = 500;                 
 confdLvl = 0.9;                       
 numGrp = 3;                           

%% Store Sobol Control Parameters and Specifications
 Sobol_inp.N = N; 
 Sobol_inp.seedNum = seedNum;
 Sobol_inp.funFile = funFile;
 Sobol_inp.funPath = funPath;
 Sobol_inp.funFile = funFile;
 Sobol_inp.smplMtd = smplMtd;
%% Randomization
 if isempty( seedNum ) == false
     rand('state',seedNum); 
 end

%% Generate the base sample from a unit hypercube

numDim = length(lb);  
SobolCost = N * (numDim + 2); 

switch smplMtd
    case 'RND'
        baseSample = rand(N, numDim * 2); 
    case 'LHS'
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
    case 'SymLHS'
        baseSample = SymLHS( N , numDim * 2,'maxmin', 10 );
    case 'PLHS'
        numSlices = ceil(N/sliceSize);
        if numSlices > 1
            smpSize = numSlices * sliceSize;
            p = PLHSdesign(smpSize, numDim * 2, numSlices, 10, 'maxmin');
            baseSample = p(1:N,:);
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        else
            baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
        end
    case 'Sobolseq'
        p = sobolset(numDim * 2,'Skip',1e3,'Leap',1e2); 
        p = scramble(p,'MatousekAffineOwen'); 
        baseSample = p(1:N,:);
    case 'Halton'
        p = haltonset(numDim * 2,'Skip',1e3,'Leap',1e2); 
        p = scramble(p,'RR2'); 
        baseSample = p(1:N,:);
    otherwise
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
end

%% Define sub-sample matrices A, B, and C, and scale them onto the actual factor space
A = baseSample(:, 1 : numDim);
B = baseSample(:, numDim + 1 : end);
for j = 1 : N
    A(j, :) = A(j, :) .* ( ub' - lb') + lb';
    B(j, :) = B(j, :) .* ( ub' - lb') + lb';
end
% A=Normalize(A);
% B=Normalize(B);
for i = 1 : numDim
    C{i} = B;
    C{i}(:, i) = A(:, i);
end
%% Run the function/model for all points in matrices A, B, and C and return the function/model response
[yA, yB, yC ] = funcEval(A, B, C, funPath, funFile);

%% Cacluate Sobol's First-Order (FO) and Total-Order (TO) sensitivity indices and rank factors based on TO
[ FO, TO, V, Mu ] = SobolCalc (yA, yB, yC);
rnkFO = factorRanking(FO);
rnkTO = factorRanking(TO);
Sobol_out.FO = FO;
Sobol_out.TO = TO;
Sobol_out.rnkFO = rnkFO;
Sobol_out.rnkTO = rnkTO;
Sobol_out.V = V;
Sobol_out.Mu = Mu;
T_rnkFO(n,:)=factorRanking(FO);
T_rnkTO(n,:)=factorRanking(TO);
T_FO(n,:)=FO;
T_TO(n,:)=TO;

%% Bootstrap
if bootstrapFlag == 1;
    [ Sobol_out.bootstrap ] = bootstrapSobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFO, rnkTO, numGrp);
end


%% Store Results
% if n==length(T_N)
save ('Results_Sobol', 'Sobol_out');
% end
end

function [ FO, TO, Vy, MUy ] = SobolCalc (yA, yB, yC)
[ N,  numDim ] = size(yC);
f0 = mean([ yA; yB ]);  
MUy = f0;
Vy = ( [ yA; yB ]' * [ yA; yB ] ) / ( 2 * N ) - f0 ^ 2; 
for i = 1 : numDim
    FO(i) = ( (( yA' * yC(:, i) ) / N) - mean(yA) * mean(yB) ) / Vy;     
    TO(i) = 1 - ( ( yB' * yC(:, i) ) / N - f0 ^ 2 ) / Vy; 
end
end



clearvars -except TO FO rnkFO rnkTO
%% ************************************************************************
function [ yA, yB, yC ] = funcEval (A, B, C, funPath, funFile)
currDir = pwd;
cd(funPath);
[ N , numDim ] = size(A);
SobolCost = N * (numDim + 2); 

yA = zeros(N, 1);
yB = zeros(N, 1);
yC = zeros(N, numDim);
fprintf('A Sobol experiment started: size of base sample = %g, total number of model runs = %g. \n', N, SobolCost );
for j = 1 : N
    fprintf('Group run (base sample) #%g started. Running model %s %g times...', j, funFile, numDim + 2 );
    tic;
    yA(j, 1) = feval(funFile, A(j, :) );
    yB(j, 1) = feval(funFile, B(j, :) );
    for i = 1 : numDim
        yC(j, i) = feval(funFile, C{i}(j, :) );
    end
    time = toc;
    fprintf(' Group run finished in %g seconds.\n', time);
end
cd (currDir);
end
%% ------------------------------------------------------------------------
function [ bootstrap ] = bootstrapSobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFOBnchmrk, rnkTOBnchmrk, numGrp)

[ N  numDim ] = size(yC);
[randstream1] = RandStream.create('mrg32k3a','NumStreams',1, 'Seed', 1234567);
for k = 1 : bootstrapSize
    rows = ceil ( rand(randstream1, N, 1) * N );
    yAboot = yA(rows);
    yBboot = yB(rows);
    yCboot = yC(rows, :);
    [ FO(k, :), TO(k, :), V(k, :), MU(k, :) ] = SobolCalc (yAboot, yBboot, yCboot);
     rnkFO(k, :) = factorRanking(FO(k, :));
     rnkTO(k, :) = factorRanking(TO(k, :));
end
FO_sorted = sort(FO, 1);
TO_sorted = sort(TO, 1);

bootstrap.FO_low = FO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.FO_upp = FO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

bootstrap.TO_low = TO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.TO_upp = TO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

for D = 1 : numDim
    bootstrap.rel_FO(1, D) = length ( find( rnkFO(:, D) == rnkFOBnchmrk(1, D) ) ) / bootstrapSize;
    bootstrap.rel_TO(1, D) = length ( find( rnkTO(:, D) == rnkTOBnchmrk(1, D) ) ) / bootstrapSize;
end
% ************ based on group ranking 
rank_benchmark_grp = groupRanking ( rnkTOBnchmrk , numGrp );
for iter = 1 : bootstrapSize
    rnkTO_grp( iter, :) = groupRanking ( rnkTO(iter, :) , numGrp );
end
for D = 1 : numDim
    relGrp(1, D) = length ( find( rnkTO_grp(:, D) == rank_benchmark_grp(1, D) ) ) / bootstrapSize;
end
 
end
%% ************************************************************************
function [ rank ] = factorRanking(SAindices)
[sorted, order] = sort(SAindices, 'descend');
temp = [ order; 1: length(SAindices) ]';
temp2 = sortrows(temp, 1)';
rank = temp2(2, :);
end
%% ************************************************************************
function rank_grp = groupRanking ( rank_indvl, numGrp )
numDim = length (rank_indvl);
grpSize = round ( numDim / numGrp );
grpNum = 1; temp = 0;
for rankNum = 1 : numDim
    if temp == grpSize
        temp = 0;
        grpNum = grpNum + 1;
    end
    temp = temp + 1;
    rank_grp ( 1,  rank_indvl == rankNum  ) = grpNum;
end
end

function Data=Normalize(data)
for i=1:size(data,2)
    Max(i)=max(data(:,i));
    Min(i)=min(data(:,i));
    Data(:,i)=(data(:,i)-Min(i))./(Max(i)-Min(i));
end
end

function [NSE]= hbv_sobol(par_values)
global d_stream
[ Q_sim ] = run_HBV ( par_values );
[NSE]  = error_model ( d_stream, Q_sim );
NSE=-NSE;
end

function [ Q_sim , ET_sim , SMS_sim , ini_values_valid, output_state] = run_HBV ( par_values )


global  d_rainfall d_temp initial_c m_temp m_prec months
watershed_area = initial_c(1); 
ini_values = [initial_c(1), initial_c(2), initial_c(3), initial_c(4)];
forcing = [d_rainfall, d_temp, months];
long_term=[m_temp, m_prec];

[ output_flux , output_state ] = HBV_SASK ( forcing , long_term , par_values , ini_values , watershed_area );
Q_sim = output_flux ( 1 : end, 1 );
ET_sim = output_flux ( 1 : end, 7 );
SMS_sim = output_state( 1 : end, 2 );
ini_values_valid= output_state(end,:) ;
end

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

function [NSE]  = error_model ( Q_obs, Q_sim )

mean_obs = mean ( Q_obs );
denominator = mean ( ( Q_obs - mean_obs * ones(length(Q_obs), 1 ) ) .^ 2 );
numerator = mean ( ( Q_obs - Q_sim ) .^ 2 );
errfun = 1 - numerator / denominator;
NSE = -1 * errfun; % used for optimization: to change problem to minimization

end

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

