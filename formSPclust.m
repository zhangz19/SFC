function [labs] = formSPclust(lab1, Y, W, W0)
% find spatial clusters based on a partition and neighborhood matrix W

% load('lab_kmeans7.mat', 'lab1');

%============== find subgroups that are spatially contiguous
labs = nan(1, numel(lab1));
k = 1;
for r = 1:max(lab1)
    indr = find(lab1==r);
    A = W(indr, indr);
    [v_per, freq] = lrcm(A);
    freq1 = [0, freq];
    % A(v_per, v_per)  %block-diagonal after sorting
    for i = 1: (numel(freq1)-1)
        labs(indr(  v_per( (freq1(i)+1):freq1(i+1) )	)) = k;
        k = k +1;
    end
end

%============== reduce number of clusters
d = max(labs);
nr =  histc(labs, 1:d); 
[~, I]  = sort(labs);  T = size(Y, 1); 
ys = mat2cell(Y(:, I), T, nr);
ym = cellfun(@(x) mean(x,2), ys, 'UniformOutput', false); 
ym = cell2mat(ym);   %plot(ym) the cluster mean
% a = dist(ym);  %distance matrix at cluster level

%----------------- absorb minor groups to its nearest cluster neighbor
nrmin = 3;  %forcely delete clusters with <= nrmin members. Caution with possible singleton
% note: a spatial cluster with 3 members have all boundary points (8NN)
inds = find(nr <= nrmin); 
labs0 = [labs, 0];   % note W0 is padded with N+1, so pad labs with 0 and later exclude 0 when counting
for r = inds  % for each minor group
    indr = find(labs==r);
    blabs = labs0(W0(indr,:));
    binds = arrayfun(@(i) unique(blabs(i,blabs(i,:)~=0)), 1:size(blabs,1),'UniformOutput',false); 
    % note, binds will include the neighbors in the same clusters
    [~,I] = sort(cellfun(@length, binds), 'descend'); % start with the "boundary" points with less purity (more l)
    
    for i = I % for each site, assign it to its
        s = indr(i);
        blabs = labs(W(s,:)==1);    uc = unique(blabs);     uc(uc==r) = []; %exclude cluster it currently belong to
        ds = dist([Y(:,s), ym(:, uc)]); ds = ds(1, 2:end); labs(s) = uc(ds==min(ds));
        ds = dist([Y(:,s), ym(:, uc)]); ds = ds(1, 2:end); labs(s) = uc(ds==min(ds));
    end
end

%----------------- finally, update labels by remonving empty clusters
nr =  histc(labs, 1:d);  k = 1; for r = 1:d; if nr(r) ~= 0; labs(labs==r) = k; k = k+1; end; end
% d = max(labs); nr =  histc(labs, 1:d);  % useful check


% save('checkFindSP.mat', 'labs')

% not run
