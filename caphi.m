function [phiz] = caphi(A,Sf,C,X3,T,f)
for k=1:size(A,1)
    Sft = Sf.*exp(T(k,:)'*f);
    S = ifft(Sft,[],2); 
    phi(k) = norm(X3(:,k) - kr(C,S.')*A(k,:).','fro')^2;
end
phiz = sum(phi);