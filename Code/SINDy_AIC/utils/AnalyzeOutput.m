%% Analyze results
minind = find(min(cell2mat({IC.aic_c})) == cell2mat({IC.aic_c}), 1, 'first')
mincoeff = numcoeff(minind);
x1coeff = Xicomb{minind}
lambdavec(minind,:)

% % TESTING
% minind2 = find(min(AIC_rel) == AIC_rel, 1, 'first')
% mincoeff2 = numcoeff(minind2)
% x1coeff2 = Xicomb{minind2}

% find models with index below or equal to relAICc = 7;
indless7 = find(AIC_rel<=7);
[AIC_rel_uni, ia, ic] = unique(AIC_rel);
for kk = 1:max(numcoeff)
    kk;
    % find indices of terms with kk coefficients
    ind = find(numcoeff==kk);
    % find unique indices in AIC_rel
    
    % find indices of unique AIC_rel valuse that have kk terms and
    % relAICc<=7
    indtog = intersect(intersect(indless7, ind), ia);
    
    if length(indtog)>0
        %     AIC_rel(indtog)
        %  figure(5)
        %     plot(kk, AIC_rel(indtog), 'og')
        Xiless7= Xicomb{indtog};
    end
end


%% plots
% plot everything on 2-D plot with the total number of terms on the x axis,
% and the aic for all models of that size on the y-axis.
    %{
    if exist('t')
        figure(111)
        plot(t,x, 'o')
        xlabel('t')
        ylabel('x')
        title('time series measurements')
        
    else
        figure(111)
        plot(x/Ntot, 'o')
        xlabel('t')
        ylabel('x')
        title('time series measurements')

    end
%}

%%
% Note: 
% minind == best model picked by aic_c metric above.
    
% Plot relative AICc
AICcRel =cell2mat({IC.aic_c})-min(cell2mat({IC.aic_c}));
   figure; hold on
    plot(numcoeff, AICcRel, 'o')
    plot(numcoeff(minind),AICcRel(minind),'om') % Plot best model by AIC_c
    xlabel('number of terms')
    ylabel('Relative AICc')
    title(['\lambda =',...
                    num2str(lambdavals.lambdastart),':',...
                    num2str(lambdavals.numlambda),':',...
                    num2str(lambdavals.lambdaend)])

%{
% Plot relative AIC
% AIC_Rel defined in EX_SEIR.m
    figure; hold on
    plot(numcoeff, AIC_rel, 'o')
     axis([min(numcoeff) max(numcoeff) 0 max(AIC_rel)])
%    axis([min(numcoeff) 20 0 max(AIC_rel)])
    xlabel('number of terms')    
    ylabel('Relative AIC')
%     patch([0 max(numcoeff)+1 max(numcoeff)+1 0], ...
%         [10 10  max(AIC_rel) max(AIC_rel)], ...
%         [0.2 0.2 0.2], 'FaceAlpha', 0.6)
%     patch([0 max(numcoeff)+1 max(numcoeff)+1 0], [4 4 7 7], [0.2 0.2 0.2], 'FaceAlpha', 0.3)
%     patch([0 max(numcoeff)+1 max(numcoeff)+1 0], [0 0 2 2], [0.2 0.2 0.2], 'FaceAlpha', 0.1)
    plot(numcoeff(minind),AIC_rel(minind),'om') % Highlight best model
    title(plotTitle)
    
% Plot relative BIC
BIC_Rel =cell2mat({IC.bic})-min(cell2mat({IC.bic}));
   figure; hold on
    plot(numcoeff, BIC_Rel, 'o')
    plot(numcoeff(minind),BIC_Rel(minind),'om') % Plot best model by AIC_c
    xlabel('number of terms')
    ylabel('Relative BIC')
%}
  
%% integrate identified systems
% this takes too long. may or may not require more debugging.
% execution in poolData line 28 when halted
%{
t = (1:(N-1))';
[tint,xint]=ode45(@(t,x)sparseGalerkin(t,x,Xicomb{minind},...
                                    polyorder,usesine),t',x(1,:));  % approximate ...set ODE options?
figure(6)
hold on
plot(tB,xint(:,1),'c--','LineWidth',1.2)
plot(tB,xint(:,2),'m--','LineWidth',1.2)
plot(tB,xint(:,3),'g--','LineWidth',1.2)
legend('S','E','I','Szeta','Ezeta','Izeta')
%}


