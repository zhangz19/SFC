function [x, w, E, W0] = initPara(Y, X, W)
% input: Y is T by N, X is T by p by N
% output: x parameters, w partition, E environment variables

E.interceptonly = 0; 
E.simu = 0; 
E.useflat = 0;
E.computeFL = 1; 
E.computeML = 0;
E.updateb = 1;  
E.updatee = 1;  
E.updateu = 0;

[T, p, N] = size(X);

p = 1; %for test run with first p covariates

% include an intercept
X = X(:,[1,1:p],:);  X(:,1,:) = 1; p = p+1;
if E.interceptonly == 1; X = X(:,1,:); p = 1; end

pT = p*T; J = floor(log2(T)); indJ = zeros(1,T); t = 1; for j = 0:J; for k = 1:round(2^(j-1)); indJ(t) = j; t = t+1; end; end;  indJ = indJ + 1;
% construct wavelet transformation matrix
tmp = eye(T); Wav = tmp; for i1 = 1:T; Wav(:,i1) = wavedec(tmp(:,i1),J,'db1'); end
Ywav = Wav*Y;

%% precalculation
% tic
Xalls = zeros(T, N, pT);
for i = 1:p
    %     for j = 1:T
    %         % %----------------- same calculation but the second is faster
    %         %         tic
    %         %         X_rij_0 = zeros(T, N);
    %         %         for s = 1:N
    %         %             tmp = (Wav*diag(X(:,i,s))*Wav');
    %         %             X_rij_0(:,s) = tmp(:,j);
    %         %         end
    %         %         toc
    %         %         tic
    %         X_rij = zeros(T, N);
    %         for t = 1:T
    %             X_rij = X_rij + ...
    %                 repmat(squeeze(X(t,i,:))', [T,1]) .*repmat(Wav(:,t)*Wav(j,t), [1,N]);
    %         end
    %         %         toc
    %         %         norm(X_rij_0-X_rij)
    %         Xalls(:,:, (i-1)*T+j) = X_rij;
    %         % figure(1); plot(X_rij(:,1),'b-'); tmp = (Wav*diag(X(:,i,1))*Wav'); hold on; plot(tmp(:,j),'r-'); hold off
    %     end
end
% toc  %Elapsed time is 25.433665 seconds.
% save('Xalls.mat', 'Xalls')
load('Xalls.mat', 'Xalls')
Xalls = Xalls(:,:,1:pT); 
if E.interceptonly == 1; Xalls = Xalls(:,:,1:pT); end

%% set hyperparameters
mean_IG = 0.01; var_IG = 10^4;
alpha_IG = 2+mean_IG^2/var_IG;
invbeta_IG = mean_IG*(alpha_IG-1);
E.A_beta = alpha_IG; %*ones(J+1,p);
E.B_beta = invbeta_IG; %*ones(J+1,p); % better not Jefferys
if E.useflat == 1;  alpha_IG = 0; invbeta_IG = 0;  end
E.A_e = alpha_IG; E.B_e = invbeta_IG;
E.A_h = alpha_IG;% *ones(T,1);
E.B_h = invbeta_IG; %*ones(T,1);
E.A_g = alpha_IG;% *ones(T,1);
E.B_g = invbeta_IG; %*ones(T,1);
% A_gu = 100; B_gu = 1; %A_gu large: closer to 1
% alpha_delta = alpha_IG; invbeta_delta = invbeta_IG;
% alpha_phi = alpha_IG; invbeta_phi = invbeta_IG;
% alphasig = 0; invbetasig = 0; %Jeffery
% alphaeta=0; invbetaeta = 0;
E.pia = 1; E.pib = 1;
E.indJ = indJ; E.p = p; E.T = T; E.J = J;

% store the nearest neighbors based on W
W0 = (N+1)*ones(N, max(sum(W==1,2)));  % pad "N+1" in the end to make same length
for i = 1:N; inds  = find(W(i,:)==1);  W0(i, 1:length(inds)) = inds;  end

%% initialize parameters
load('lab_kmeans7.mat', 'lab1');
labs = formSPclust(lab1, Y, W, W0);
w.labs = labs;  w.d = max(labs);  nr = histc(w.labs, 1:w.d);
[~, I]  = sort(w.labs);
Ywav = mat2cell(Ywav(:, I), T, nr);
% tic
Xalls = Xalls(:,I, :);   Xalls = mat2cell(Xalls, T, nr, pT);
% toc

u = zeros(T, N);
x = cell(1, w.d); % for parallelization
for r = 1:w.d
    indr = find(w.labs==r);
    if numel(indr)==0, indr = w.center(r); end
    % x{r}.indr = indr;
    x{r}.beta = zeros(T, p);
    for t = 1:T % this is in original domain
        Xt = squeeze(X(t,:,indr));
        if E.interceptonly == 0
            Xt = Xt'; %Xt is now nr by p
        end
        x{r}.beta(t,:) = (Xt'*Xt + 1e-4*eye(size(Xt,2)))\(Xt'*Y(t,indr)');
    end
    Yhat = Wav*( squeeze(sum( X(:,:,indr).*repmat(x{r}.beta, [1,1,nr(r)]), 2 )) ); %response surface
    x{r}.beta = Wav*x{r}.beta; %transform beta into wavelet domain as well
    if numel(indr) == 1; x{r}.h = 1e-4*ones(T,1); else; x{r}.h = var( Ywav{r} - Yhat, [], 2); end
    x{r}.phi = 0.9*ones(T,1);
    x{r}.g = x{r}.h*.9; % x.u = 0
    if E.updateu == 0;   x{r}.g = x{r}.h;  end% without considering the random effects
    x{r}.h = x{r}.h - x{r}.g;
    
    x{r}.sigma2 = x{r}.g(1);
    x{r}.g = x{r}.g/x{r}.sigma2;
    x{r}.h = x{r}.h/x{r}.sigma2;
    if E.updateu == 1
        for t = 1:T
            if E.simu == 1
                Lo = diag(M(indr)) - x{r}.phi(t)*W(indr,indr); Lo = chol(Lo, 'lower');
                u(t, indr) = (Lo'\normrnd(0,sqrt(x{r}.h(t)*x{r}.sigma2), [nr(r),1]))';
            else
                u(t, indr) = Ywav{r}(t) - Yhat(t,:) - ...
                    normrnd(0,sqrt(x{r}.g(t)*x{r}.sigma2),[1,nr(r)]);
            end
        end
    end
    x{r}.u = u(:,indr);
    x{r}.Y = Ywav{r};
    x{r}.err = Ywav{r} - Wav*( squeeze(sum( X(:,:,indr).*repmat(Wav'*x{r}.beta, [1,1,nr(r)]), 2 )) )  - x{r}.u;
    x{r}.X = Xalls{r};
    x{r}.n = nr(r);
    x{r}.lambda = 1e4*ones(J+1,p);    %var(Ywav,[],2);
    x{r}.pis = 0.0001*ones(J+1, p);   %betarnd(pia, pib, [p J+1]);
    x{r}.pis(1,:) = 1;
end

end

% not run
