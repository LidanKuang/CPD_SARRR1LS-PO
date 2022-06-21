function A = krprod(B,C);
% Khatri Rao product
sb = size(B,1);
sc = size(C,1);
A = zeros(sb*sc,size(B,2));
for k = 1:size(B,2)
    A(:,k) = reshape(C(:,k)*B(:,k).',sb*sc,1);
end