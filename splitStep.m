function [x,L1] = splitStep(x,L0)
global N nrmin p pT usepenalty T useExp D W J coords updateu
global plotit Ywav Wav X sigma_prop Gd_prior %A_g B_g A_gu B_gu B_beta A_beta pia pib

x1 = x; x1.d = x.d + 1; L_null = L0;

% newCenter = randsample(setdiff(1:N, x.center),1); % site ID

% [~, ties] = ClusterGen(x.center); ties = find(ties); % sample from ties
% newCenter = ties(randsample(length(ties),1));

% note x1.center is no longer ordered, so just add the new cluster to the
% end, i.e., ind = x1.d = x.d+1, so that older labels match in x and x1
x1.center = nan(1, x1.d);
ind = x1.d;
ind0 = randsample(x.d, 1); % the cluster to be splitted, must has min distance to ind
%ind = randsample(x1.d, 1);

% myd = D(x.center(ind0), :); 
% myind = find(myd == min(myd(myd>0))); 
myind = find(D(x.center(ind0), :)==1 & x.labs==ind0);
newCenter = myind(randsample(length(myind), 1));

x1.center(ind) = newCenter; i0 = setdiff(1:x1.d, ind);
x1.center(i0) = x.center; 
[x1.labs, x1.bds] = ClusterGen4Split(x1.center, x.labs, ind); 

PRatio = 0;
% need to compute the proposal ratio for the proposal above

% correct boundary with the changing cluster under x1
ratio1 = 0;
labs = zeros(1,N); %indicate if it is a boundary in x1
for i = find(~cellfun(@isempty, x1.bds))
    if any(x1.bds{i}==ind) || any(x1.bds{i}==ind0) % boundary that involves the two splitted clusters
        %         if isempty(x.bds{i}) % this is new boundary, still keep old labels
        %             if x.labs(i) ~= ind
        %                 x1.labs(i) = x.labs(i);
        %                 %randsample(x1.bds{i}, 1); % this is safe when x1.bds{i} contains at least 2 entries (should be)
        %                 labs(i) = 1;
        %             end
        if isempty(x.bds{i}) % this is new boundary, assign to the nearest cluster or new cluster if it has ties
            DistMi = D(i, x1.center);
            minind = find(abs(DistMi - min(DistMi)) < 1e-10);
            if length(minind)==1
                x1.labs(i) = minind;
            else
                x1.labs(i) = ind;
                if all(x1.bds{i} ~= ind)
                    warning('new boundary has tie and does not involve the new cluster, suspicious')
                end
            end
        elseif x.labs(i) == ind0 && any(x1.bds{i} == ind) % old boundary only involves ind, not ind0 && all(x1.bds{i} ~= ind0) 
            x1.labs(i) = ind; % ind inherits its father (older ind0 before its split) boundaary assignment
            % note x1.bds{i} may also include ind0. 
        else % belong to old boundary, keep the label
            x1.labs(i) = x.labs(i);
        end
    else % boundary value does not involve two splitted clusters, usually from x that is far away from the splitted center
        x1.labs(i) = x.labs(i); % note: i may not be a boundary in x
    end
end

ratio2 = 0;
labs0 = zeros(1,N); %indicate if it is a new boundary in x
for i = find(~cellfun(@isempty, x.bds)) %calculate the part of proposal ratio involved for x
    if isempty(x1.bds{i}) || x1.labs(i) == ind  % "new" boundary in x, will do boundary correction as in merge step
        C = length(x.bds{i}); % size of choice set
        lik = zeros(1,C);
        for j = 1:C
            r = x.bds{i}(j);
            err = Ywav(:,i) - ( Wav*( sum( X(:,:,i).*(Wav'*x.beta(:,:,r)), 2) )  + x.u(:,i) );
            lik(j) = -0.5*sum(log(x.sigma2(r)*x.g(:,r))) - 0.5/x.sigma2(r)*sum(err.^2./x.g(:,r));
        end
        lik = exp(lik-max(lik)); lik = lik/sum(lik);
        ratio2 = ratio2 - (Gd_prior==2)*log(C) - log(lik(x.bds{i}==x.labs(i)));
        labs0(i) = 1;
    end
end

PRatio = PRatio + ratio1 - ratio2;



for myplotit = 1:plotit
    cmat = jet(x1.d);
    figure(2), sz = 15*ones(1,N); sz(x.center) = 50;
    scatter(coords(:,1),coords(:,2),sz,cmat(x.labs,:),'filled');
    indr = x1.center(ind);
    hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+90,cmat(x1.d,:),'*'); hold off
    indr = x1.center(ind0);
    hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+90,cmat(ind0,:),'s'); hold off
    tmpind = find(x1.labs ~= x.labs);
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind),cmat(x1.labs(tmpind),:),'filled'); hold off
    tmpind = find(labs0 ==1); % old boundary (in x)
    hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+60,cmat(x.labs(tmpind),:),'s'); hold off
    tmpind = find(labs==1); % newly added boundary (in x1, not in x)
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+60,cmat(x1.labs(tmpind),:),'d'); hold off
end


nr = histc(x1.labs, 1:x1.d);
x1.beta = zeros([T, p, x1.d]); x1.beta(:,:,i0) = x.beta;
x1.lambda = zeros(J+1,p,x1.d); x1.lambda(:,:,i0) = x.lambda;
x1.pis = zeros(J+1,p,x1.d); x1.pis(:,:,i0) = x.pis;
x1.g = zeros(T,x1.d); x1.g(:,i0) = x.g;
x1.h = zeros(T,x1.d); x1.h(:,i0) = x.h;
x1.phi = zeros(T,x1.d); x1.phi(:,i0) = x.phi;
x1.sigma2 = zeros(1,x1.d); x1.sigma2(i0) = x.sigma2;

% % propose g for the new cluster
indx = find(x1.labs==ind);
weights = histc(x.labs(indx), 1:x.d); weights = weights/sum(weights);
ghat = sum(repmat(weights, [T,1]).*x.g, 2);
% % this proposal is not reversible
% gnull = B_g./(A_g+1); gnull(1) = 1; %prior mode
% U_g = betarnd(A_gu, B_gu, [T,1]); % for all T g-compoenents
% x1.g(:,ind) = ghat.*U_g + gnull.*(1-U_g);
% Jacob = sum(log(abs(ghat(2:end) - gnull(2:end)))) - sum(log(betapdf(U_g(2:end), A_gu, B_gu)));
U_g = normrnd(0, sigma_prop, [T-1,1]);
x1.g(1,ind) = 1; x1.g(2:end,ind) = exp(U_g).*ghat(2:end,:);
Jacob = sum(log(abs(x1.g(:,ind)))) + 0.5*(T-1)*log(2*pi) + (T-1)*log(sigma_prop) + 0.5*sum(U_g.^2)/sigma_prop^2;

%propose lambda and pis for wavelets
ghat = 0; ghat2 = 1;
for r = 1:x.d
    ghat = ghat + weights(r)*x.lambda(:,:,r);
    ghat2 = ghat2 .* (1./x.pis(2:end,:,r) - 1).^weights(r);
end
U_g = normrnd(0, sigma_prop, [2*J+1,p]); J1 = numel(U_g);
Jacob = Jacob + 0.5*J1*log(2*pi) + J1*log(sigma_prop) + 0.5*sum(sum(U_g.^2))/sigma_prop^2; % - pi(u)
U_g = exp(U_g);
x1.lambda(:,:,ind) = ghat.*U_g(1:(J+1),:);
Jacob = Jacob + sum(sum(log(abs(x1.lambda(:,ind))))) ;
ghat2 = ghat2.*U_g((J+2):end,:);
Jacob = Jacob - 2*sum(sum(log(abs( 1+ghat2 )))) + sum(sum(log(abs( ghat2 )))) ;
x1.pis(1,:,ind) = 1;
x1.pis(2:end,:,ind) = 1./( 1+ghat2 );

if updateu == 1 %update random effects
    % propose h
    ghat = sum(repmat(weights, [T,1]).*x.h, 2);
    U_g = normrnd(0, sigma_prop, [T,1]);
    x1.h(:,ind) = exp(U_g).*ghat;
    Jacob = Jacob + sum(log(abs(x1.h(:,ind)))) + 0.5*T*log(2*pi) + T*log(sigma_prop) + 0.5*sum(U_g.^2)/sigma_prop^2;
    
    % propose phi
    ghat = prod(( (1-x.phi)./(x.phi+1) ).^repmat(weights, [T,1]), 2); % b_phi = 1, a_phi = -1
    U_g = normrnd(0, sigma_prop, [T,1]);
    ghat = ghat.*exp(U_g); %now ghatexp(u)
    x1.phi(:,ind) = 2./(1+ghat) - 1; %(b-a)/(1+ghatexp(u)) + a
    Jacob = Jacob + sum(log( 2*ghat./(1+ghat).^2 )) + ... % ghatexp(u)(b-a)/(1+ghatxep(u))^2
        0.5*T*log(2*pi) + T*log(sigma_prop) + 0.5*sum(U_g.^2)/sigma_prop^2;
end

% for myplotit = 1:plotit
%     cmat = jet(x1.d);
%     figure(2), sz = 15*ones(1,N); sz(x1.center) = 50;
%     scatter(coords(:,1),coords(:,2),sz,cmat(x1.labs,:),'filled');
%     sz(x1.center(ind)) = 50;
%     hold on; scatter(coords(indx,1),coords(indx,2),sz(indx)+50,[0,0,0],'d'); hold off
% end

if(min(nr) < nrmin)
    x.alpha = -99; x.growth = 0; L1 = L0;
else
    [L1, x1] = getLoglike(x1, 1, ind); % marginal likelihood
    L0 = getLoglike0_r(x, indx);  % full likelihood
    
    for myplotit = 1:plotit
        cmat = jet(x1.d);
        figure(3); plot(Wav'*Ywav, 'k.');
        hold on
        for r = 1:x1.d
            indr = find(x1.labs==r); nr=length(indr);
            tmp = Wav*( squeeze(sum( X(:,:,indr).*repmat(Wav'*x1.beta(:,:,r), [1,1,nr]), 2 )) ) + x1.u(:,indr);
            plot(Wav'*tmp, 'LineStyle','-','Color',cmat(r,:)); %cols{r}
        end
        hold off
    end
    
    if usepenalty == 1
        logratio = L1 - L0 - 0.5*pT*log(N*T);
    elseif useExp == 1
        logratio = L1 - L0 - tan(pi/2*x1.kappa);
    else
        logratio = L1 - L0 + log(1-x1.kappa);
    end
    logratio = logratio + PRatio + Jacob;
    alpha = min(logratio,0);  u = log(rand(1)); growth = 0;
    if u <= alpha
        x = x1; growth = 1;
        L1 = getLoglike0_r(x1, 1:N); % full likelihood
    else
        L1 = L_null;
    end
    x.alpha = alpha;  x.split = growth;
end
end
