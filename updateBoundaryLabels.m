function [w1] = updateBoundaryLabels(x1, w1, u, Ywav, plotit, X, Wav, T)
% 2015-7-30: simply update the boundary labels
global verbose coords N W D
w = w1; 
useBIC = 0;

mylabs = zeros(1,N); 
for i = find(~cellfun(@isempty, w1.bds))
    C = length(w1.bds{i}); % size of choice set
    lik = zeros(1,C);
    
    %     % check a weird thing
    %     inds = find(W(i,:) == 1);
    %     if isempty(setdiff(inds, find(~cellfun(@isempty, x1.bds)) )) && sum(x.labs(inds) ~= x.labs(i))==3
    %         mylabs(i) = 1;
    %     end
    
    err = zeros(T, C); 
    for j = 1:C
        r = w1.bds{i}(j);
        err(:,j) = Ywav(:,i) - ( Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  + u(:,i) );
        lik(j) = -0.5*sum(log(x1{r}.sigma2*x1{r}.g)) - 0.5/x1{r}.sigma2*sum(err(:,j).^2./x1{r}.g);
        if useBIC == 1
            lik(j) = lik(j) - 0.5*log(T)*sum(reshape(x1{r}.beta, [1, numel(x1{r}.beta)])~=0); 
        end
    end
    
    
    
    %     % 2015-7-31: check the updating boundary values
    %     cmat = jet(w1.d);
    %     figure(2), sz = 15*ones(1,N); sz(w.center) = 50;
    %     scatter(coords(:,1),coords(:,2),sz,cmat(w.labs,:),'filled');
    %     indr = w.center;
    %     hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+50,cmat,'filled'); hold off
    %     tmpind = i;
    %     hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+40,cmat(w1.labs(tmpind),:),'ks'); hold off
    %
    %     figure(4), plot( Wav'*Ywav(:,i),'k.')
    %     myps = zeros(1,C);
    %     err1 = zeros(T, C);
    %     for j = 1:C
    %         r = w1.bds{i}(j);
    %         err1(:,j) = Wav*( sum( X(:,:,i).*(Wav'*x1{r}.beta), 2) )  + u(:,i);
    %         hold on, plot(Wav'*err1(:,j),'Color',cmat(r,:)); hold off
    %         myps(j) = sum(reshape(x1{r}.beta, [1, numel(x1{r}.beta)])~=0);
    %     end
    
    
    
    
    
    lik = exp(lik-max(lik)); lik = lik/sum(lik);
    i1 = sum(cumsum([0, lik(1:(C-1))]) <= rand(1));
    w1.labs(i) = w1.bds{i}(i1);
    
    %     if x1.labs(i) ~= x.labs(i)
    %         tmpind = i;
    %         hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+40,cmat(x1.labs(tmpind),:),'kd'); hold off
    %     end
    
end
% changed = any(w.labs ~= w1.labs); 

% if verbose ==1
%     fprintf('L[%1d]  ', changed)
% end

for myplotit = 1:plotit
    cmat = jet(w1.d);
    figure(2), sz = 15*ones(1,N); sz(w.center) = 50;
    scatter(coords(:,1),coords(:,2),sz,cmat(w.labs,:),'filled');
    indr = w.center;
    hold on; scatter(coords(indr,1),coords(indr,2),sz(indr)+50,[0,0,0],'filled'); hold off
    tmpind = find(~cellfun(@isempty, w.bds));
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+40,cmat(w.labs(tmpind),:),'kd'); hold off
    tmpind = find(w1.labs ~= w.labs);
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind),cmat(w1.labs(tmpind),:),'filled'); hold off
    hold on; scatter(coords(tmpind,1)+.4,coords(tmpind,2),sz(tmpind)+40,cmat(w1.labs(tmpind),:),'ks'); hold off
    tmpind = find(mylabs==1);
    hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+60,cmat(w1.labs(tmpind),:),'k+'); hold off
    tmpind = find(~cellfun(@isempty, w1.bds));
    hold on; scatter(coords(tmpind,1),coords(tmpind,2),sz(tmpind)+20,cmat(w1.labs(tmpind),:),'k*'); hold off
end

end