
function [x,L1] = mergeStep(x,L0)
global pT nrmin usepenalty N T useExp D J p plotit coords Wav Ywav X
global sigma_prop  %A_g B_g A_gu B_gu B_beta A_beta pia pib

Dcen = D(x.center, x.center); Dcen = triu(Dcen); Dcen2 = Dcen(Dcen>0); indmin = find(Dcen2==min(Dcen2));
ind = randsample(1:length(indmin),1); ind = indmin(ind);
dicy = (2:x.d).*(1:(x.d-1))/2; j = find(dicy>=ind,1,'first'); i = j - (dicy(j)-ind); j = j+1;
DelCenter = randsample(x.center([i,j]),1);
% DelCenter = randsample(x.center,1);
ind = find(x.center==DelCenter);
ind0 = [i,j]; ind0(ind0 == ind) = [];

L_null = L0;
x1 = x; 
x1.center(ind) = []; 

[x1.labs, x1.bds] = ClusterGen(x1.center, []); 
x1.d = length(x1.center);
x1.beta(:,:,ind) = [];
x1.h(:,ind) = [];
x1.phi(:,ind) = [];
x1.sigma2(ind) = [];
lambdastar = x1.lambda(:,:,ind); x1.lambda(:,:,ind) = [];
pistar = x1.pis(2:end,:,ind); x1.pis(:,:,ind) = [];
gstar = x1.g(:,ind); x1.g(:,ind) = [];

mylabs = x.labs; mylabs(mylabs==ind) = -1;
mylabs(mylabs>ind) = mylabs(mylabs>ind) - 1; % x1 labels, but matched to labels in x1
mylabs(mylabs== -1) = x.d; 

PRatio = 0;
% correct boundary with the changing cluster under x1
ratio1 = 0;
labs = zeros(1,N); %those with boundary correction
for i = find(~cellfun(@isempty, x1.bds)) %calculate the part of proposal ratio involved for x
    if isempty(x.bds{i}) || x.labs(i) == ind % new boundary created in x1, do boundary correction
        % plus: old boundary but assigned to the
        % deleted cluster. But how about reversibility?
        C = length(x1.bds{i}); % size of choice set
        lik = zeros(1,C);
        for j = 1:C
            r = x1.bds{i}(j);
            err = Ywav(:,i) - ( Wav*( sum( X(:,:,i).*(Wav'*x1.beta(:,:,r)), 2) )  + x1.u(:,i) );
            lik(j) = -0.5*sum(log(x1.sigma2(r)*x1.g(:,r))) - 0.5/x1.sigma2(r)*sum(err.^2./x1.g(:,r));
        end
        lik = exp(lik-max(lik)); lik = lik/sum(lik);
        i1 = sum(cumsum([0, lik(1:(C-1))]) <= rand(1));
        x1.labs(i) = x1.bds{i}(i1);
        ratio1 = ratio1 + log(1/C) - log(lik(i1));
        labs(i) = 1;
    else % old boundary assigned to other clusters
        x1.labs(i) = mylabs(i); 
    end
end


ratio2 = 0;
labs0 = zeros(1,N); %"new" boundary in x
for i = find(~cellfun(@isempty, x.bds))
    if isempty(x1.bds{i}) 
        labs0(i) = 1;
    end
end

PRatio = PRatio + ratio1 - ratio2;


for myplotit = 1:plotit
    cmat = jet(x.d);
    figure(2), sz = 15*ones(1,N); sz(x.center) = 50;
    scatter(coords(:,1),coords(:,2),sz,cmat(mylabs,:),'filled');
    indr = x.center(ind);
    hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+90,[0,0,0],'d'); hold off
    indr = x.center(ind0);
    hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+60,[0,0,0],'ks'); hold off
    tmpind = find(x1.labs ~= mylabs);
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind),cmat(x1.labs(tmpind),:),'filled'); hold off
    tmpind = find(labs0 ==1);
    hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+60,cmat(mylabs(tmpind),:),'s'); hold off
    tmpind = find(labs==1);
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+60,cmat(x1.labs(tmpind),:),'d'); hold off
end


% propose g, lambda, and pi for the new cluster
indx = find(x.labs==ind);
weights = histc(x1.labs(indx), 1:x1.d); weights = weights/sum(weights);
ghat = sum(repmat(weights, [T,1]).*x1.g, 2);
% gnull = B_g./(A_g+1); gnull(1) = 1; %prior mode
% U_g = (gstar-gnull)./(ghat-gnull); %betarnd(A_gu, B_gu, [T,1]); % for all T g-compoenents
% Jacob = sum(log(abs(ghat(2:end) - gnull(2:end)))) - sum(log(betapdf(U_g(2:end), A_gu, B_gu)));
U_g = log(gstar(2:end)./ghat(2:end)); %normrnd(0, sigma_prop, [T,1]);
Jacob = sum(log(abs( gstar(2:end) ))) + 0.5*(T-1)*log(2*pi) + (T-1)*log(sigma_prop) + 0.5*sum(U_g.^2)/sigma_prop^2;

ghat = 0; ghat2 = 1;
for r = 1:x1.d
    ghat = ghat + weights(r)*x1.lambda(:,:,r);
    ghat2 = ghat2 .* (1./x1.pis(2:end,:,r) - 1).^weights(r);
end
U_g = log( lambdastar./ghat );  %normrnd(0, sigma_prop, [2*J+1,p]);
eq = 1./pistar-1;
U_g = [U_g; log(eq) - log(ghat2) ];
J1 = numel(U_g);
Jacob = Jacob + 0.5*J1*log(2*pi) + J1*log(sigma_prop) + 0.5*sum(sum(U_g.^2))/sigma_prop^2; % - pi(u)
Jacob = Jacob + sum(sum(log(abs( lambdastar )))) ;
Jacob = Jacob - 2*sum(sum(log(abs( 1+eq )))) + sum(sum(log(abs( eq )))) ;

ghat = 0; ghat2 = 0;
for r = 1:x1.d
    ghat = ghat + weights(r)*x1.lambda(:,:,r);
    ghat2 = ghat2 + weights(r)*x1.pis(:,:,r);
end
% gnull = B_beta./(A_beta+1); %IG prior mode
% gnull2 = ones(J+1, p); gnull2(2:end,:) = pia./(pia+pib); %beta prior mean
% U_g = (lambdastar-gnull)./(ghat-gnull); %betarnd(A_gu, B_gu, [J+1, p]); % for all T g-compoenents
% U_g2 = (pistar-gnull2)./(ghat2-gnull2); %betarnd(A_gu, B_gu, [J+1, p]); % for all T g-compoenents
% Jacob = Jacob + sum(sum(log(abs(ghat - gnull)))) - sum(sum(log(betapdf(U_g, A_gu, B_gu))))...
%     + sum(sum(log(abs(ghat2(2:end,:) - gnull2(2:end,:))))) - sum(sum(log(betapdf(U_g2(2:end,:), A_gu, B_gu)))) ;


% for myplotit = 1:plotit
%     cmat = jet(x.d);
%     cmat(ind,:) = [0,0,0];
%     figure(2), sz = 15*ones(1,N); sz(x1.center) = 50;
%     scatter(coords(:,1),coords(:,2),sz,cmat(x1.labs,:),'filled');
% end

nr = histc(x1.labs, 1:x1.d);
if(min(nr) < nrmin)
    x.alpha = -99; x.merge = 0; L1 = L0;
else
    % [L1,x1] = getLoglike(x1, 1);
    L1 = getLoglike0_r(x1, indx); % full likelihood
    L0 = getLoglike_r(x, ind, 0);  % marginal likelihood
    
    if usepenalty == 1
        logratio = L1 - L0 + 0.5*pT*log(N*T);
    elseif useExp == 1
        logratio = L1 - L0 + tan(pi/2*x1.kappa);
    else
        logratio = L1 - L0 - log(1-x1.kappa);
    end
    logratio = logratio + PRatio; % - Jacob;
    alpha = min(logratio,0);
    u = log(rand(1)); merge = 0;
    if u <= alpha
        x = x1; merge = 1;
        L1 = getLoglike0_r(x1, 1:N); % full likelihood
    else
        L1 = L_null;
    end
    x.alpha = alpha;
    x.merge = merge;
end
end