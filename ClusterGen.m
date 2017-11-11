
function [labs, bds] = ClusterGen( CenterIndex, labs0 )
% 2015-7-23: 
% return cluster labels and boundary set with label assigned 0
% bds has non-empty cell indicating the neighboring labels of each boundary
% point
global N D K
labs = NaN(1, N);
labs(CenterIndex) = 1:length(CenterIndex);
bds = cell(1, N); 
for i = setdiff(1:N,CenterIndex)
    DistMi = D(i, CenterIndex);
    ind = find(abs(DistMi - min(DistMi)) < 1e-10);
    if length(ind) > 1 % distance tie
        bds{i} = ind;
        labs(i) = 0;
    else
        if isempty(labs0)
            labs(i) = ind;
        else
            labs(i) = labs0(i); % inherit the info from labs0
            % warning: in this case, it could be no longer a voronoi tessellation
            % boundary set is an independent set from voronoi tessellation
        end
    end
end
labs2 = ones(1,N);
% note the following step can rely on the provided labs0
for i = 1:N % center could be boundary point
    vec = D(i,:); 
    % ind = find(vec<=K & vec>0); %K-neighborhood
    vec = unique(labs(vec<=K & vec>0));
    vec = vec(vec~=0); % new rule: append boundary if its K-neighors, except tied points, have different labels
    C = length(vec);
    if C > 1 % append boundary set
        %         vec2 = []; tmp = vec(vec~=0); 
        %         if length(tmp)>1
        %             vec2 = tmp;
        %         end
        %         if any(vec==0)
        %             for j = find(labs(ind)==0)
        %                 vec2 = [vec2, bds{ind(j)}]; % neighboring labels of its neighbor with tie
        %             end
        %         end
        %         bds{i} = unique(vec2);
        bds{i} = unique(vec);
        labs2(i) = 0;
    end
end
labs(labs2==0) = 0;
end
