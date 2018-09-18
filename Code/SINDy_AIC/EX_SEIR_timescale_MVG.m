% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

addpath('utils')
addpath('models')
changeplot
clear all, close all, clc
%clear all, clc

filename = 'SEIR';
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

eps=  0.00025; % noise
numvalidation = 100; % number of crossvalidation experiments
format shortE % Set notation to be able to read small numbers in output
tic
%% generate Data
n = 4; % Number of equations
N  = 250; % number of time steps



% All other parameters are zero
% Transfer parameters
B_SE = 0.3; % Infectious rate: between 0.3 and 0.5 doesn't break
B_EI = 0.4; % Incubation rate
B_IR = 0.04; % Recovery rate
% Vital parameters
B_S = 0.02; % Birth rate
B_SEIR = 0.02; % Death rate

% Total population
Ntot = 1e4; 

% Initial Conditions
S(1) = 0.99*Ntot; % Susceptible
E(1) = 0.01*Ntot; % Exposed
I(1) = 0;         % Infected
P(1) = 0;

%%
% Random walk precipitation
rng(6) % Set seed
N=250; % From SEI code
f = 4; % Number of high-frequency obs. per health model timestep 
g = 0; % Number of high-frequency timesteps of overlapping history you want

p = cumsum(randn(1,N*f+g)); 
figure; plot(p)
title('Precip random walk')
legend('seed = 6')
% Center and standardize to +/- 1
range = max(p) - min(p);
pCenter = p-max(p)+0.5*range; %plot(pCenter)
pStand = pCenter./max(pCenter); %plot(pStand)
%title('P (standardized to +/- 1)')
p = pStand; % Reassign p to pStand
pickInd = find(rem(1:N*f,f)==0); % Indices corresponding to health model timesteps 

% Make f timestep sum rolled P
p = p(g+1:N*4); % For this method, throw away leading overlapping history
pStack = zeros(4,N*4); % Initialize
pStack(:,4:end) = [p(4:end); % Based on f
                   p(3:end-1);
                   p(2:end-2);
                   p(1:end-3)];
% Sum over past f timesteps for only health model timesteps
pSumN = sum(pStack(:,pickInd),1);
figure; plot(pSumN); title('pSumN')
% Standardize pSumN to the range taken by B_SE = +/- 0.1
pSumRanged = 0.1.*(pSumN/max(pSumN)); % plot(pSumRanged) 
% Use the pSumRanged as input for our model
P = pSumRanged;

% % Make time delay p series for different lags
% % for f = 4 and g = 4:
% pickInd = pickInd + g; % Results in same picked values as in rolled version
% P = 0.1*(p./8); % Make maximum possible value of total P B_SE friendly 
% pt1 = P(pickInd);
% pt2 = P(pickInd - 1);
% pt3 = P(pickInd - 2);
% pt4 = P(pickInd - 3);
% % pt5 = P(pickInd - 4);
% % pt6 = P(pickInd - 5);
% % pt7 = P(pickInd - 6);
% % pt8 = P(pickInd - 7);


% % Not using derivative for now.
% % Calculate derivative time series (p = pStand)
% dp = [0,(p(2:end)-p(1:end-1))]; % Piece-wise slope
% figure; plot(dp)
% title('d(pStand)')

%%
plotTitle = 'SEIRvital';
% disease tranfer model
for ii =2:N
    
%     % SEIR model, static pop., B_SE = fn(P_tlags)
%     S(ii) = S(ii-1) - (B_SE+...
%                         0.01*pt1(ii-1)+pt2(ii-1)+0.05*pt3(ii-1)+0.08*pt4(ii-1))*...
%                         S(ii-1)*I(ii-1)/Ntot;
%     E(ii) = E(ii-1) + (B_SE+...
%                         0.4*pt1(ii-1)+0.015*pt2(ii-1)+pt3(ii-1)+0.07*pt4(ii-1))*...
%                         S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1);
%     I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);
%     % adding in the R data causes SINDy to fail.
    
        % SEIR model, static pop., B_SE = fn(sum(P)), 
    %  where sum is taken from last four timesteps
    S(ii) = S(ii-1) - (B_SE*S(ii-1)*I(ii-1)/Ntot + S(ii-1)*I(ii-1)*P(ii-1)/Ntot);
    E(ii) = E(ii-1) + (B_SE*S(ii-1)*I(ii-1)/Ntot + S(ii-1)*I(ii-1)*P(ii-1)/Ntot)...
                        - B_EI*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);    
    
%     % SEIR model, static pop.
%     S(ii) = S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
%     E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1);
%     I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);
%     % adding in the R data causes SINDy to fail.
    
%     % SEIR model, vital dynamics
%     S(ii) = S(ii-1) + B_S*Ntot - B_SEIR*S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
%     E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1) - B_SEIR*E(ii-1);
%     I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1) - B_SEIR*I(ii-1);
% %    R(ii) = R(ii-1) + B_IR*I(ii-1) - B_SEIR*R(ii-1);
    
end

% create x and dx matrices with all variables:
% %
% x =  [S(1:end-1)' E(1:end-1)' I(1:end-1)'...
%       pt1(1:end-1)' pt2(1:end-1)' pt3(1:end-1)' pt4(1:end-1)'];
%    %   pt5(1:end-1)' pt6(1:end-1)' pt7(1:end-1)' pt8(1:end-1)'];
% dx = [S(2:end)'   E(2:end)'   I(2:end)'...
%       pt1(2:end)'   pt2(2:end)'   pt3(2:end)'   pt4(2:end)'];
%  %     pt5(2:end)'   pt6(2:end)'   pt7(2:end)'   pt8(2:end)'];
%
x = [S(1:end-1)' E(1:end-1)' I(1:end-1)' P(1:end-1)'];
dx = [S(2:end)' E(2:end)' I(2:end)' P(2:end)']; % ??? NOT THE DERIVATIVE, NO?

% x = [S(1:end-1)' E(1:end-1)' I(1:end-1)'];
% dx = [S(2:end)' E(2:end)' I(2:end)'];

% add noise to state variables
rng(10);
x = x+eps*randn(size(x));

% Plot specified data
if plottag>=1
    figure(6)
    plot(x(:,1:3)/Ntot, 'o')
    xlabel('time step')
    ylabel(['% of population size: ' num2str(Ntot)])
    legend('S', 'E', 'I')
%    legend('S', 'E', 'I', 'P')
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
lambdavals.numlambda = 3;
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
%    R = zeros(N,1); % Recovered
    
    % Initial Conditions
    S(1) = x0cross(1, jj); % number of suceptibles in population
    E(1) = x0cross(2, jj);
    I(1) = x0cross(3, jj);
%    R(1) = x0cross(4, jj);
    
    for ii =2:N
        

        % SEIR model, static pop., B_SE = fn(sum(P)), 
    %  where sum is taken from last four timesteps
    S(ii) = S(ii-1) - (B_SE*S(ii-1)*I(ii-1)/Ntot + S(ii-1)*I(ii-1)*P(ii-1)/Ntot);
    E(ii) = E(ii-1) + (B_SE*S(ii-1)*I(ii-1)/Ntot + S(ii-1)*I(ii-1)*P(ii-1)/Ntot)...
                        - B_EI*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);    
    
        
        
%     S(ii) = S(ii-1) + B_S*Ntot - B_SEIR*S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
%     E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1) - B_SEIR*E(ii-1);
%     I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1) - B_SEIR*I(ii-1);
    end
    % create x and dx matrices with all variables:
    x2= [S(1:end-1) E(1:end-1) I(1:end-1)];
%    x2= [S(1:end-1) E(1:end-1) I(1:end-1) R(1:end-1)];
    xA{jj} = x2 +eps*randn(size(x2));
    dxA{jj} = [S(2:end) E(2:end) I(2:end)]';
%    dxA{jj} = [S(2:end) E(2:end) I(2:end) R(2:end)]';
    
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