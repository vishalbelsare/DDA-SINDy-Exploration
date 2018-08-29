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
tic
%% generate Data
n = 4; % Number of equations

% all others are zero
% Transfer parameters
B_SE = 0.3; % Infectious rate
B_EI = 0.4; % Incubation rate
B_IR = 0.04; % Recovery rate
% Vital parameters
B_S = 0.02; % Birth rate
B_SEIR = 0.02; % Death rate
Ntot = 1e4; % Total (initial) population

N  = 250; % number of time steps

% Initial Conditions
S(1) = 0.99*Ntot; % number of suceptibles in population
E(1) = 0.01*Ntot;
I(1) = 0;
R(1) = 0;

plotTitle = 'SEIRvital';
% disease tranfer model
for ii =2:N
    %{
    % SEIR model, static pop.
    S(ii) = S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
    E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);
    % adding in the R data causes SINDy to fail.
    %}
    % SEIR model, vital dynamics
    S(ii) = S(ii-1) + B_S*Ntot - B_SEIR*S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
    E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1) - B_SEIR*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1) - B_SEIR*I(ii-1);
    R(ii) = R(ii-1) + B_IR*I(ii-1) - B_SEIR*R(ii-1);
    
end

% create x and dx matrices with all variables:
x = [S(1:end-1)' E(1:end-1)' I(1:end-1)' R(1:end-1)'];
dx = [S(2:end)' E(2:end)' I(2:end)' R(2:end)'];
% add noise to state variables
rng(10);
x = x+eps*randn(size(x));

% Plot specified data
if plottag>=1
    figure(6)
    plot(x/Ntot, 'o')
    xlabel('time step')
    ylabel(['% of population size: ' num2str(Ntot)])
    legend('S', 'E', 'I', 'R')
    legend('boxoff')
    title(plotTitle)
end

%% pool Data
polyorder = 2;  % search space up to 3rd order polynomials
usesine = 0;    % no trig functions
laurent = 0;
Theta = poolData(x,n,polyorder,usesine, laurent);
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
[Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals);

% display coefficients
for ii = 1:length(Xicomb)
    Xicomb{ii}
end

%% calculate validation data for new intial conditions.

 x0cross = 10.^(-1 + (4+1)*rand(n,numvalidation));


for jj = 1:numvalidation
    % initialize
    S = zeros(N,1); % Susceptible
    E = zeros(N,1); % Latent/exposed
    I = zeros(N,1); % Infected
    R = zeros(N,1); % Recovered
    
    % Initial Conditions
    S(1) = x0cross(1, jj); % number of suceptibles in population
    E(1) = x0cross(2, jj);
    I(1) = x0cross(3, jj);
    R(1) = x0cross(4, jj);
    
    for ii =2:N
    S(ii) = S(ii-1) + B_S*Ntot - B_SEIR*S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
    E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1) - B_SEIR*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1) - B_SEIR*I(ii-1);
    R(ii) = R(ii-1) + B_IR*I(ii-1) - B_SEIR*R(ii-1);
    end
    % create x and dx matrices with all variables:
    x2= [S(1:end-1) E(1:end-1) I(1:end-1) R(1:end-1)];
    xA{jj} = x2 +eps*randn(size(x2));
    dxA{jj} = [S(2:end) E(2:end) I(2:end) R(2:end)]';
    
end
val.x0 = x0cross;
val.tA = N;
val.xA = xA;
val.options= 0;

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

toc