function [A,B,C]=cp3_init(X,R,init_type)
%CP3_INIT Initialization of the loading matrices for the CP3 decomposition of X
%extracted from tensorlab
% INPUTS: 
% - X: 3rd order tensor of size (IxJxK)
% - R: Number of rank-1 components in the PARAFAC model
% - init_type='dtld' to initialize by Direct Trilinear Decomposition if possible,
%   otherwise randomly,
% - init_type='random' to enforce random initialization.
% OUTPUTS:
% A(IxR) B(JxR) and C(KxR): Loading Matrices

[I,J,K]=size(X);
% Check if 2 dimensions are greater than the rank to see if DTLD can be used
[size_sort,perm_vec]=sort([I J K],'descend');
if ((size_sort(2))<R) 
    init_type='random';    % DTLD can not be used so use random init
end
        
if strcmp(init_type,'dtld')==1
    [A,B,C] = cp3_dtld(X,R);    
elseif strcmp(init_type,'random')==1
    if isreal(X); A=rand(I,R);B=rand(J,R);C=rand(K,R);
    else;A=rand(I,R)+j*rand(I,R);B=rand(J,R)+j*rand(J,R);C=rand(K,R)+j*rand(K,R);
    end
    % Orthogonal matrices if possible else well conditioned matrices
    if I>=R;A=orth(A);else;A=orth(A')';end
    if J>=R;B=orth(B);else;B=orth(B')';end
    if K>=R;C=orth(C);else;C=orth(C')';end
end
end