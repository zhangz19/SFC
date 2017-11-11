function [v_per, freq] = lrcm( A )
% from https://www.mathworks.com/examples/matlab/community/20074-l-rcm-matrix-block-detection-algorithm
%Calculation of the Laplacian
L = sparse(diag( sum(A,2) ) - A);
%Applies RCM to identify the blocks
v_per = symrcm( L );
%Permute the Laplacian to separate the blocks
Lp = L(v_per,v_per);
%Sum the permuted laplacian to identify the blocks
value = sum( triu(Lp) );
freq = find( value == 0 );
end