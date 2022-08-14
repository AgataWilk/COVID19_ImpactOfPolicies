%% 
% Initialize values
p = [0.2605 0.1020];
countries = populations_europe.CountryName;
dates_train = datetime(2020,2,1): datetime(2020,11,30);
dates_test = datetime(2020,12,1): datetime(2021,1,31);
dates_total = datetime(2020,1,22): datetime(2021,1,31);
nCountries = 42;
nParams = 14;
% Length of training period and whole time horizon
skip = 10; % skip 10 days for restriction function estimation as the initial b_opt estimates are unstable
trainingLength = length(dates_train);
totalLength = length(dates_total);
% Construct training set
beta_train = beta_opt(ismember(beta_opt.Date, dates_train),:);
restrictions_train = restrictions{ismember(restrictions.Date, dates_train),3:end};

% Number of days for all countries 
daysTotal = size(beta_opt,1);
daysTraining = size(beta_train,1);
%%
% Common model 
bCountries_common = ones(daysTraining,1);
BCountries_common = ones(daysTotal,1);

% Restriction function
restrictionsMat_train = [ones(daysTraining,1),-restrictions_train];

% estimate parameters
x_common = lsqnonlin(@(x)nlinfun(x,restrictionsMat_train,bCountries_common,beta_train.beta,"common"),zeros(1,nParams));
beta_est_common = nlinbeta(x_common,[ones(daysTotal,1),-restrictions{:,3:end}],BCountries_common,"common");

%%
% Independent models
bCountries_independent = ones(trainingLength,1);
BCountries_independent = ones(totalLength,1);
x_independent = zeros(14,nCountries);
beta_est_independent = zeros(daysTotal,1);
for country = 1:nCountries
    beta_train_country = beta_train.beta(((country-1)*trainingLength+1):country*trainingLength);
    restrictions_train_country = restrictions_train(((country-1)*trainingLength+1):country*trainingLength,:);
    restrictionsMat_train_country = [ones(trainingLength,1),-restrictions_train_country];
    x_country = lsqnonlin(@(x)nlinfun(x,restrictionsMat_train_country,bCountries_independent,beta_train_country,"independent"),zeros(1,14));
    x_independent(:,country) = x_country;

    beta_est_independent(((country-1)*totalLength+1):country*totalLength) = nlinbeta(x_independent,[ones(totalLength,1),-restrictions{(((country-1)*totalLength+1):country*totalLength),3:end}],BCountries_independent,"independent");

end

%%
% Individualized models
bCountries_individualized = zeros(daysTraining,nCountries);
for k = 1 : nCountries
    bCountries_individualized(((k-1)*trainingLength+1):k*trainingLength,k) = 1;
end

BCountries_individualized = zeros(daysTotal,nCountries);
for k = 1 : nCountries
    BCountries_individualized(((k-1)*totalLength+1):k*totalLength,k) = 1;
end

x_individualized = lsqnonlin(@(x)nlinfun(x,restrictionsMat_train,bCountries_individualized,beta_train.beta,"individualized"),zeros(1,55));

beta_est_individualized = nlinbeta(x_individualized,[ones(daysTotal,1),-restrictions{:,3:end}],BCountries_individualized,"individualized");

%%
% Estimate daily infections using beta_opt for the training period and
% beta_est for the test period
yd_common = zeros(nCountries,totalLength);
yd_independent = zeros(nCountries,totalLength);
yd_individualized = zeros(nCountries,totalLength);

for k = 1:nCountries
N = populations_europe.Population(k);
Beta = zeros(totalLength,1);
Beta(1:(trainingLength+skip)) = beta_opt.beta(((k-1)*totalLength+1):((k-1)*totalLength+trainingLength+skip),:);
Beta((trainingLength+skip+1):totalLength) = beta_est_common(((k-1)*totalLength+trainingLength+skip+1):(k*totalLength),:);
[t,y] = ode15s(@(t,y)SEIR(t,y,Beta,N,p,(totalLength-1)),[0:1:(totalLength-1)],[N,0,1,0]);
yd_common(k,:) = p(1)*y(:,2);
Beta((trainingLength+skip+1):totalLength) = beta_est_independent(((k-1)*totalLength+trainingLength+skip+1):(k*totalLength),:);
[t,y] = ode15s(@(t,y)SEIR(t,y,Beta,N,p,(totalLength-1)),[0:1:(totalLength-1)],[N,0,1,0]);
yd_independent(k,:) = p(1)*y(:,2);
Beta((trainingLength+skip+1):totalLength) = beta_est_individualized(((k-1)*totalLength+trainingLength+skip+1):(k*totalLength),:);
[t,y] = ode15s(@(t,y)SEIR(t,y,Beta,N,p,(totalLength-1)),[0:1:(totalLength-1)],[N,0,1,0]);
yd_individualized(k,:) = p(1)*y(:,2);
end
%%
% Observed daily infections
infections_daily = reshape(pandemic_europe.ConfirmedDaily,totalLength,nCountries)';

%% Estimate errors
% RMSE
RMSE_common = sqrt(mean((infections_daily(:,(trainingLength+skip+1):end) - yd_common(:,(trainingLength+skip+1):end)).^2,2));
RMSE_independent = sqrt(mean((infections_daily(:,(trainingLength+skip+1):end) - yd_independent(:,(trainingLength+skip+1):end)).^2,2));
RMSE_individualized = sqrt(mean((infections_daily(:,(trainingLength+skip+1):end) - yd_individualized(:,(trainingLength+skip+1):end)).^2,2));

% Normalize by mean infections
NRMSE_common = RMSE_common./mean(infections_daily(:,(trainingLength+skip+1):end),2);
NRMSE_independent = RMSE_independent./mean(infections_daily(:,(trainingLength+skip+1):end),2);
NRMSE_individualized = RMSE_individualized./mean(infections_daily(:,(trainingLength+skip+1):end),2);

% Result
NRMSE = table(NRMSE_common,NRMSE_independent,NRMSE_individualized,'RowNames',countries,'VariableNames',{'Common','independent','individualized'});
NRMSE_mat = table2array(NRMSE);

