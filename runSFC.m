function [] = runSFC()

load('datLM.mat', 'Y','X','W', 'W0','coords')

% if isempty(gcp('nocreate')); parpool('local'); end

sfc(Y, X, W); 

% if ~isempty(gcp('nocreate'));    poolobj = gcp('nocreate');    delete(poolobj); end

end

% not run

