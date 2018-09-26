% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

addpath('utils')
addpath('models')
changeplot
clear all, close all, clc

filename = 'SEIR';
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

eps=  0.00025; % noise
numvalidation = 100; % number of crossvalidation experiments

%% generate Data

N  = 250; % number of time steps
nLags = 2; % number of lags
% Specify how many variables to use
nExo = 1 + nLags; % number of exogenous variables included in poolData
nFunc = 3; % number of equations to keep for model ID

% Exogenous variable
P = linspace(0,10,N+nLags)';
% Lag variables
P2 = P(1:(end-nLags));
P1 = P(nLags:(end-(nLags-1)));
P = P((nLags+1):end);
exo = [P P1 P2];

% Initial Conditions
Ntot = 1e4; % Total population
S(1) = 0.99*Ntot; % number of suceptibles in population
E(1) = 0.01*Ntot;
I(1) = 0;

plotTitle = 'Toy Model';
for ii = 2:N
% Make synthetic data
    
    % Toy Model: Brine tank cascade
    S(ii) = S(ii-1) - 0.5*S(ii-1);
    E(ii) = E(ii-1) + 0.5*S(ii-1) - 0.25*E(ii-1);
    I(ii) = I(ii-1) + 0.25*E(ii-1) - 0.7*I(ii-1) + ...
                        0.5*P(ii-1) + 1.0*P1(ii-1) + 4.0*P2(ii-1);

end

%% Put it in matrices

xFull = [S' E' I' exo];
x = xFull(1:end-1,:);
dx = xFull(2:end,:);

% add noise to state variables
rng(10);
x = x+eps*randn(size(x));

if plottag>=1
    figure(6)
    plot(x(:,1:4)/Ntot, 'o')
    xlabel('time step')
    ylabel(['% of population size: ' num2str(Ntot)])
    legend('S', 'E', 'I', 'P')
    legend('boxoff')
end

%% pool Data
polyorder = 1;  % search space up to how many order polynomials
usesine = 0;    % no trig functions

Theta = poolData(x,polyorder,usesine,0);
m = size(Theta,2);

% % % normalize the columns of Theta
% for k=1:m
%     normTheta(k) = norm(Theta(:,k));
%     Theta(:,k) = Theta(:,k)/normTheta(k);
% end

Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

%% compute Sparse regression
lambdavals.numlambda = 20;
lambdavals.lambdastart = -10;
lambdavals.lambdaend = 1;

% find coefficient vectors
[Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals, nFunc);

% display coefficients
for ii = 1:length(Xicomb)
    Xicomb{ii}
end

%% calculate validation data for new intial conditions.

 x0cross = 10.^(-1 + (4+1)*rand(nFunc,numvalidation));


for jj = 1:numvalidation
    % initialize
    S = zeros(N,1); % susceptibles
    E = zeros(N,1); % Latent/exposed
    I = zeros(N,1); % infected
    
    % Initial Conditions
    S(1) = x0cross(1, jj); % number of suceptibles in population
    E(1) = x0cross(2, jj);
    I(1) = x0cross(3, jj);
    
    for ii =2:N
    S(ii) = S(ii-1) - 0.5*S(ii-1);
    E(ii) = E(ii-1) + 0.5*S(ii-1) - 0.25*E(ii-1);
    I(ii) = I(ii-1) + 0.25*E(ii-1) - 0.7*I(ii-1) + ...
                        0.1*P(ii-1) + 0.05*P1(ii-1) + 0.025*P2(ii-1);

    end
    % create x and dx matrices with all function variables:
    x2= [S(1:end-1) E(1:end-1) I(1:end-1)];
    xA{jj} = x2 +eps*randn(size(x2));
    dxA{jj} = [S(2:end) E(2:end) I(2:end)]';
    
end
val.x0 = x0cross;
val.tA = N;
val.xA = xA;
val.options= 0;
val.exo = exo;

%% Perform validation and calculate aic for each model

clear abserror RMSE tB xB IC
for nn = 1:length(Xicomb)
    Xi = Xicomb{nn};
    clear error RMSE1 savetB savexB
    [error, RMSE1, savetB, savexB] = validateXi(Xi, Thetalib, val, plottag);
    ICtemp = ICcalculations(error', numcoeff(nn), numvalidation);
    abserror(:,nn) = error';
    RMSE(:,nn) = RMSE1;
    tB{nn} =savetB;
    xB{nn} = savexB;
    IC(nn) = ICtemp;
end

AIC_rel =cell2mat({IC.aic})-min(cell2mat({IC.aic}));
% plot number of terms vs AIC plot
AnalyzeOutput

rmpath('utils')
rmpath('models')