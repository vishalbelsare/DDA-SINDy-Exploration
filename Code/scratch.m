%Making data

clear all

t = 0:500;

p = 2 + 0.01*t + cos(0:.1:(.1*500)); % Oscillator + line
d = 4.*(p(2:end)-p(1:end-1)); % Piece-wise slope
s = 2+0.6.*(d.*p(1:end-1)); % Slope .* p

figure
hold on
plot(p)
plot(d)
plot(s)
legend('p','d','s')


rng(1); 
r = rand(1,5);
R.sum = sum(r);
R.mean = mean(r);
R.min = min(r);
R.max = max(r);

% Make and plot a random walk:
figure
hold on
for i = [    5 6 7  9 ]
    rng(i)
    p = cumsum(randn(1,1000));
    plot(p)
    i
    max(p) - min(p)
end
legend('5', '6', '7', '9')
%%    
    
% Random walk precipitation
rng(6) % Set seed
N=250; % From SEI code
p = cumsum(randn(1,N*4));
figure; plot(p)
title('Precip random walk')
legend('rng = 6')
% Center and standardize to +/- 1
range = max(p) - min(p);
pCenter = p-max(p)+0.5*range; %plot(pCenter)
pStand = pCenter./max(pCenter); %plot(pStand)
%title('P (standardized to +/- 1)')
p = pStand; % Reassign p to pStand

% Make 4 timestep sum rolled P
pStack = zeros(4,N*4); % Initialize
pStack(:,4:end) = [p(4:end);
                   p(3:end-1);
                   p(2:end-2);
                   p(1:end-3)];
% Pick out only the indices corresponding to health model timesteps               
pickInd = find(rem((1:1000),4)==0);
pSumN = sum(pStack(:,pickInd),1);
figure; plot(pSumN); title('pSumN')
% Standardize pSumN to the range taken by B_SE = +/- 0.1
pSumRanged = 0.1.*(pSumN/max(pSumN)); % plot(pSumRanged) 

% % Not using derivative for now.
% % Calculate derivative time series (p = pStand)
% dp = [0,(p(2:end)-p(1:end-1))]; % Piece-wise slope
% figure; plot(dp)
% title('d(pStand)')

%%
% Probing AIC inequalities


figure; hold on; 
plot(cell2mat({IC.aic_c}))
plot(cell2mat({IC.aic}),'g')
plot(cell2mat({IC.bic}),'r')
legend('aic_c','aic','bic')

