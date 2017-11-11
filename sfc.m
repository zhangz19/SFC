function [out] = sfc(Y, X, W)
% Y is T by N, X is T by p by N, lab is N-vector of labels
global E

%% prepare data
[x, w, E, W0] = initPara(Y, X, W); 

%% MCMC
verbose = 1; tot = 100;	burn = 50;  thin = 1; nsample = (tot-burn)/thin;
filename = 'out.mat'; 
out = cell(1, nsample);

if verbose == 1; fprintf('MCMC for SFC is now running...\n'); end
tic
iter_save = 1; 
for iter = 1:tot
    if verbose == 1; fprintf('%4d ', iter); if(~mod(iter,20)); fprintf('\n'); end; end
    
    x = cellfun(@gibbsStep, x, 'UniformOutput', false);
    
    if iter > burn &&  ~mod(iter-burn, thin)
        a = cell2mat(x);
        out{iter_save}.beta = reshape([a.beta], [size(x{1}.beta), w.d]); %x.beta;
        out{iter_save}.labs = w.labs; 
        iter_save = iter_save + 1; 
    end
    
end
toc

save(filename, 'out')

end



