function [x] = gibbsStep(x)
global E

n = x.n; T = E.T; p = E.p; J = E.J; indJ = E.indJ;

%================== UPDATE beta components
for myloop = 1:E.updateb
    
    %----------------------- update beta
    for i = 1:p %p %1:p
        for j = 1:T
            ij = (i-1)*T+j; indj = indJ(j);
            tmp = x.X(:,:,ij);
            x.err = x.err +  tmp*x.beta(j,i); % exclude old contribution of beta(j,i,r)
            Lo = sum(sum( tmp.^2./repmat(x.g, [1,n]) ));
            Mu = sum(sum( tmp.*x.err./repmat(x.g, [1,n]) ));
            logBF = - 0.5*log( 1+x.lambda(indj,i)*Lo );
            Lo = Lo + 1/x.lambda(indj,i);  Mu = Mu/Lo;   Lo = 1/Lo*x.sigma2;
            logBF = logBF + 0.5*Mu^2/Lo;
            Odds = exp(logBF)*x.pis(indj,i)/(1-x.pis(indj,i));
            alpha = 1;   if ~isinf(Odds);  alpha = Odds/(1+Odds);  end
            u = rand(1);
            if u <= alpha
                %if x.beta(j,i) ~=0
                x.beta(j,i) = normrnd(Mu, sqrt(Lo));
                x.err = x.err - tmp*x.beta(j,i); % include new contribution pf beta(j,i,r)
                %end
            else
                x.beta(j,i) = 0; % the contribution is 0
            end
        end
    end
    
    %----------------------- update scaling
    for i = 1:p
        for j = 0:J
            vec = x.beta(indJ==j+1,i);
            x.lambda(j+1,i) = 1./gamrnd(0.5*sum(vec~=0)+E.A_beta, ... %A_beta(j+1,i)
                (0.5*sum(vec.^2)/x.sigma2+E.B_beta)^(-1));
            if j > 0 %fix pi at j==0 to be 1
                vec = (vec~=0);
                x.pis(j+1,i) = betarnd(E.pia+sum(vec), E.pib+sum(1-vec)); %pia(j,i)
            end
        end
    end
    
end

%================== UPDATE random effects components
for myloop = 1:E.updateu
    x.err = x.err + x.u; % the error = Y - Xbeta - u now becomes Y-Xbeta
    for t = 1:T
        % tic
        Lo = (diag(x.M) - x.phi(t)*x.W)/x.h(t) + sparse(eye(n))/x.g(t); % precision
        Lo = chol(Lo, 'lower');
        % toc
        
        Mu = Lo\( x.err(t,:)'/x.g(t) );
        Mu = Mu + normrnd(0, sqrt(x.sigma2), size(Mu));
        u0 = Lo'\Mu;
        x.u(t,:) = u0';
        
        % update dependence
        u0 = x.u(t,:)';
        wDw = u0'*x.W*u0;
        loglike = lphi0 + (0.5/(x.h(t)*x.sigma2))*wDw*phis;
        MaxLogLike = max(loglike);
        P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
        % plot(phis, P)
        U0 = rand(1);
        cump = cumsum([0, P(1:(end-1))]);
        i0 = sum(U0 > cump);
        x{r1}.phi(t) = phis(1);
        if i0 > 1
            x{r1}.phi(t) = phis(i0-1) + gap_phi/P(i0)*(U0-cump(i0));
        end
        
        % update scaling
        wDw = sum(M{r1}.*u0.^2) - x{r1}.phi(t)*wDw;
        x{r1}.h(t) = 1./gamrnd(A_h + 0.5*nr(r1), (B_h + 0.5*wDw/x{r1}.sigma2)^-1 );
        
        
    end
    
    % figure(1), plot(x0.u,'r-'); hold on, plot(x.u,'b--'); hold off
    % figure(1), plot(x0.f,'r-'); hold on, plot(x.h,'b--'); hold off
    
    x.err = x.err - x.u; % Y-Xbeta becomes Y - Xbeta - u
end

%================== UPDATE error components
for myloop = 1:E.updatee
    err2 = x.err.^2; %squared residual
    
    %----------------------- update scaling
    for t = 2:T	% for t = 1, x.g(1,r) = 1 fixed
        x.g(t) = 1./gamrnd( E.A_g + 0.5*n,  ( E.B_g + 0.5*sum(err2(t,:))/x.sigma2 )^-1 );
    end
    
    %----------------------- update variation
    A = 0.5*( n*T + sum(sum(x.beta~=0))  ); %+ nr*T
    B = 0.5*(  sum(sum( err2, 2)./x.g) + ...
        sum(sum( x.beta.^2./x.lambda(indJ,:) ))  ); %mynorm2( x.u(:,indr)./repmat(x.h(:,r), [1,nr]) )
    
    if E.updateu == 1 % if there is random effect, need add this term
        A = A + 0.5*n*T;
        B = B + 0.5*sum(sum( x.u.^2./repmat(x.h, [1,n]) ));
    end
    
    x.sigma2 = 1./gamrnd( E.A_e + A,  ( E.B_e + B )^-1 );
end

%================== COMPUTE full likelihood
for myloop = 1:E.computeFL
    x.L0 = - 0.5*T*log(2*pi*x.sigma2) - 0.5*sum(log(x.g)) - 0.5*sum( err2./repmat(x.sigma2*x.g, [1, n]), 1);
    x.RMSE = sum(sum(err2));
end

%================== COMPUTE marginal likelihood
for myloop = 1:E.computeML
    myL1 = 0; myL2 = 0;
    ind1 = find(x.beta~=0);
    XX = 0; XY = 0;
    
    for t = 1:T
        Xt = squeeze(x.X(t,:,:)); %nr by pT
        Xt = Xt(:, ind1);
        if E.updateu == 1 %with random effect u
            Lo = eye(n)*x.g(t)  + ((diag(x.M) - x.phi(t)*x.W)/x.h(t))\eye(n); % covariance
            Lo = chol(Lo, 'lower'); myL2 = myL2 - sum(log(diag(Lo)));
            tmp = Lo\Xt; XX = XX + tmp'*tmp;
            Yt = Lo\x.Y(t,:)'; XY = XY + tmp'*Yt;
            myL1 = myL1 + sum(Yt.^2);
        else  %without random effect u
            XX = XX + Xt'*Xt/x.g(t);
            XY = XY + Xt'*x.Y(t,:)'/x.g(t);
            myL1 = myL1 + sum(x.Y(t,:).^2/x.g(t));
            myL2 = myL2 - 0.5*n*log(x.g(t));
        end
    end
    
    lambda = repmat(x.lambda, [1, E.p]); lambda = lambda(ind1);
    Lo = XX + diag(1./lambda);
    Lo = chol(Lo, 'lower');
    Mu = Lo\XY;
    n1 = 0.5*n*T + E.A_e;
    myL1 = 0.5*( myL1 - sum(Mu.^2) );
    myL2 = myL2 + sum(log(diag(Lo))) - 0.5*sum(log(lambda));
    if E.useflat == 0
        myL1 = -n1*log( 1 + myL1/E.B_e );
        myL2 = myL2 + gammaln(n1) - gammaln(E.A_e) - 0.5*n*T*log(2*pi*E.B_e);
    else
        myL1 = -n1*log( myL1 );
        myL2 = myL2 + gammaln(n1) -0.5*n*T*log(2*pi);
    end
    x.L1 = myL1 + myL2; % sum over clusters
end

end

% not run
