function [yout] = poolData(yin,polyorder,usesine, laurentorder)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data:
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1); % timesteps
nVars = size(yin,2); % number of variables
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
funstr{ind,1}  = '1';
ind = ind+1;


% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
    
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
            
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

thresh = 2e-7;

if(laurentorder>=1)
    for i=1:nVars
        yout(yin(:,i)>=thresh,ind) = 1./yin(yin(:,i)>=thresh, i);
        yout(yin(:,i)<thresh, ind) = 0;
        ind = ind+1
    end
end

if(laurentorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            ijthresh = yin(:,i)>=thresh & yin(:,j)>=thresh;
            yout(ijthresh,ind) = 1./yin(ijthresh,i)./yin(ijthresh,j);
            yout(~ijthresh,ind) = 0;
            ind = ind+1;
            
        end
    end
end
if(laurentorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                ijkthresh = yin(:,i)>=thresh & yin(:,j)>=thresh & yin(:,k)>=thresh;
                yout(ijkthresh,ind) = 1./yin(ijkthresh,i)./yin(ijkthresh,j)./yin(ijkthresh,k);
                yout(~ijkthresh, ind)= 0;
                ind = ind+1;
            end
        end
    end
end

if(usesine)
    for k=1:10;
        yout = [yout sin(k*yin) cos(k*yin)];
    end
end