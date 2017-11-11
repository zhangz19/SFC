function [x, w, L1] = driftStep(x, w, u, L0, W, plotit, X, Wav, Ywav, T, Gd_prior)
% no need to propose beta under new cluster
% 2015-8-28: note now Ywav and u are T by N matrix, not structured
global N nrmin coords  D 
x1 = x; w1 = w; 
subCenter = w.center(arrayfun(@(idx) ~all(ismember(find(W(idx,:)), w.center)), w.center));
nsize = length(subCenter); % number of centers with >=1 non-center neighbors
target = subCenter(randsample(nsize,1));
ind = find(w.center==target);
group = find(W(target, :)); % neighbors of selected center
group(ismember(group, w.center) |  w.labs(group)~=ind) = []; % move to non-center neighbor with the same label
gsize = length(group); % number of non-center neighbors with the same labels
if gsize == 0
    warning('something is wrong')
end

% mylen = nan(1,gsize);
% for i0 = 1:gsize
%     % i0 = randsample(gsize,1);
%     newCenter = group(i0);
%     PRatio = log(nsize) + log(gsize);
%     ind = find(x.center==target);
%     x1.center(ind) = newCenter;  % replace
%     [x1.labs, x1.bds] = ClusterGen(x1.center, []); %membership may change
%     labs2 = 0; 
%     for i = find(~cellfun(@isempty, x1.bds))
%         if any(x1.bds{i}==ind) 
%             if isempty(x.bds{i})
%                 labs2 = labs2 + 1;
%             end
%         end
%     end
%     mylen(i0) = labs2;
% end
% i0 = find(mylen == max(mylen));
i0 = randsample(gsize,1);
newCenter = group(i0);
PRatio = log(nsize) + log(gsize);
w1.center(ind) = newCenter;  % replace
[w1.labs, w1.bds] = ClusterGen(w1.center, []); %membership may change   x.labs

% correct boundary with the changing cluster under x1
ratio1 = 0;
labs = zeros(1,N); % =1 indicates corrected old boundary (match in x and x1)
labs2 = zeros(1,N); % =1 indicates new boundary
labs3 = zeros(1,N); % =1 indicates corrected old boundary (dismatch)
for i = find(~cellfun(@isempty, w1.bds))
    if any(w1.bds{i}==ind) % boundary correction iff changing center is involved
        C = length(w1.bds{i}); % size of choice set
        % j0 = w.labs(i); % the current label of i
        lik = zeros(1,C);
        for j = 1:C
            r = w1.bds{i}(j);
            % err = Ywav{j0}(:,i0) - Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  - x1{j0}.u(:,i0) ;
            err = Ywav(:,i) - Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  - u(:,i) ;
            lik(j) = -0.5*sum(log(x1{r}.sigma2*x1{r}.g)) - 0.5/x1{r}.sigma2*sum(err.^2./x1{r}.g);
        end
        loglik = lik; 
        lik = exp(lik-max(lik)); lik = lik/sum(lik);
        i1 = sum(cumsum([0, lik(1:(C-1))]) <= rand(1));
        w1.labs(i) = w1.bds{i}(i1);
        if ~isempty(w.bds{i}) % boundary under both x and x1
            if isempty(setdiff(w.bds{i},w1.bds{i})) && isempty(setdiff(w1.bds{i},w.bds{i}))
                % choice set matches. In this case, only need to consider
                % those with labels corrected, since no change yields ratio 1
                if w1.labs(i) ~= w.labs(i)
                    % note this will be canceled in L1-L0
                    ratio1 = ratio1  + loglik(w.bds{i}==w.labs(i)) - loglik(i1) + ...
                        (Gd_prior==2)*(log(length(w.bds{i})) - log(C)); 
                    labs(i) = 1;
                end
            else  % warning('boundaries do not match')
                ratio1 = ratio1  - log(lik(i1)) - (Gd_prior==2)*log(C);
                labs3(i) = 1;
            end
        elseif isempty(w.bds{i}) % it is a new boundary
            ratio1 = ratio1 - log(lik(i1)) - (Gd_prior==2)*log(C);
            labs2(i) = 1;
        end
    else % keep old label when boundary value does not involve ind, shall not contribute to ratio1
        w1.labs(i) = w.labs(i); % note: i may not be a boundary in x
    end
end

ratio2 = 0; % proposal density for LB under x
labs0 = zeros(1,N); %indicate if it is a boundary in x
for i = find(~cellfun(@isempty, w.bds)) %calculate the part of proposal ratio involved for x
    if isempty(w1.bds{i}) || labs3(i)==1 % compute "new" or dismatched boundary only! difference in boundary correction in ratio1 already
        C = length(w.bds{i}); % size of choice set
        % j0 = w.labs(i); % the current label of i
        % i0 = find( indr{j0} == i); % find the relative index of i in cluster j0
        lik = zeros(1,C);
        for j = 1:C
            r = w.bds{i}(j);
            % err = Ywav{j0}(:,i0) - ( Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  + x1{j0}.u(:,i0) );
            err = Ywav(:,i) - ( Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  + u(:,i) );
            lik(j) = -0.5*sum(log(x1{r}.sigma2*x1{r}.g)) - 0.5/x1{r}.sigma2*sum(err.^2./x1{r}.g);
        end
        lik = exp(lik-max(lik)); lik = lik/sum(lik);
        ratio2 = ratio2 + log(lik(w.bds{i}==w.labs(i))) + (Gd_prior==2)*log(C);
        labs0(i) = 1;
    end
end

PRatio = PRatio + ratio1 + ratio2;

for myplotit = 1:plotit
    cmat = jet(w1.d);
    figure(2), sz = 15*ones(1,N); sz(w.center) = 50;
    scatter(coords(:,1),coords(:,2),sz,cmat(w.labs,:),'filled');
    inds = w1.center(ind);
    hold on; scatter(coords(inds,1),coords(inds,2),sz(inds)+50,[0,0,0],'filled'); hold off
    tmpind = find(w1.labs ~= w.labs);
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind),cmat(w1.labs(tmpind),:),'filled'); hold off
    tmpind = find(labs0 ==1); % old boundary (in x)
    hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+60,cmat(w.labs(tmpind),:),'s'); hold off
    tmpind = find(labs==1); % corrected boundary
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+60,cmat(w1.labs(tmpind),:),'d'); hold off
    tmpind = find(labs3==1); % new boundary
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+60,cmat(w1.labs(tmpind),:),'k+'); hold off
    tmpind = find(labs2==1); % new boundary
    % tmpind = find(~cellfun(@isempty, x1.bds));
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+60,'ks'); hold off
end

% x1.labs = labs;



% indr1 = cell(1,w1.d);
% for r = 1:w1.d
%     indr1{r} = find(w1.labs == r);
% end


nr = histc(w1.labs, 1:w1.d);
if(min(nr) < nrmin)
    x.alpha = -99;  x.shift = 0; L1 = L0;
else
    % [L1,x1] = getLoglike(x1, 0);
    L1 = getLoglike0_r(x1, w1, u, 1:N, Ywav, X, T, Wav, plotit); 
    subCenter = w1.center(arrayfun(@(idx) ~all(ismember(find(W(idx,:)), w1.center)), w1.center));
    Revnsize = length(subCenter);
    group = find(W(newCenter, :));
    group(ismember(group, w1.center)) = [];
    Revgsize = length(group);
    PRatio = PRatio - log(Revnsize) - log(Revgsize);
    % L1 - L0 should match following:
    % tmpind = find(x1.labs ~= x.labs); getLoglike0_r(x1, tmpind) - getLoglike0_r(x, tmpind)
    logratio = L1 - L0 + PRatio; alpha = min(0,logratio);
    u = log(rand(1)); drift = 0;
    if u <= alpha
        x = x1; w = w1; drift = 1;
    else
        L1 = L0;
    end
    w.alpha = alpha; w.drift = drift;
end


% % post-hoc check
% for i = 1:N
%     labs = unique(x.labs(W(i, :)~=0));
%     if all( labs ~= x.labs(i)) % i is surrounded by neighbors from other clusters, suspicious
%         C = length(labs); % size of choice set
%         lik = zeros(1,C);
%         for j = 1:C
%             r = labs(j);
%             err = Ywav(:,i) - ( Wav*( sum( X(:,:,i).*(Wav'*x.beta(:,:,r)), 2) )  + x.u(:,i) );
%             lik(j) = -0.5*sum(log(x.sigma2(r)*x.g(:,r))) - 0.5/x.sigma2(r)*sum(err.^2./x.g(:,r));
%         end
%         lik = exp(lik-max(lik)); lik = lik/sum(lik);
%         i1 = sum(cumsum([0, lik(1:(C-1))]) <= rand(1));
%         x.labs(i) = labs(i1);
%     end
% end


nr = histc(w.labs, 1:w.d);
if any(nr ==0)% this only happens when a new cluster has all points as boundary, and corrected to other clusters
    warning('some cluster vanishes in drift step')
    ind = find(nr == 0);
    w.center(ind) = []; [w.labs, w.bds] = ClusterGen(w.center, []); w.d = length(w.center);
    x(ind) = [];
end
end
