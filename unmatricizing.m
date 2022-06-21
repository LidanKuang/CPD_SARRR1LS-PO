% -------------------------------------------------------------------------
function X = unmatricizing(X,n,D)
% Inverse function of matricizing
ind = 1:length(D);
ind(n) = [];
if n == 1
    perm = [1:length(D)];
else
    perm = [2:n(1) 1 n(end)+1:length(D)];
end
  
X = permute(reshape(X,D([n 1:n(1)-1 n(end)+1:length(D)])),perm);
