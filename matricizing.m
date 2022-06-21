function Y = matricizing(X,n)
% Turn tensor to matrix along n'th mode
if isempty(n)
   Y = X;
else
   sX = size(X);
   N = ndims(X);
   n2 = setdiff(1:N,n);
   Y = reshape(permute(X, [n n2]),prod(sX(n)),prod(sX(n2)));
end