function [S,B,C,T] = SARRR1LSPO(X,R)
% X: the complex-valued multi-subject fMRI data
% R: the number of components
% S: shared spatial maps (SMs)
% B: shared time courses (TCs)
% C: subject-specific intensities
% T: time delays
[S,B,C]=cp3_init(X,R,'random'); % randomly initilize the loading matrices
[I J K] = size(X);T = zeros(K,R);
maxiter = 200; % the maximum number of iterations
X1 = matricizing(X,1).';  % mode-1 unfloding 
vox = 1/3; sigma = 4; lambda = 2.5; % parameters for phase sparsity constraint
for i=1:maxiter
    for n=1:R
     FACT{1} = C(:,n);          % subject intensity 
     FACT{2} = B(:,n);          % shared TC
     [m] = reconstructShiftCPc(FACT,T(:,n)); m=m';
     M(:,n) = m(:);            % new aggregating mixing matrix
    end
    S = (pinv(M)*X1).'; % update shared SMs by rank-R least-square fit
    %%%%%%% add phase sparsity constriant %%%%%%%%%%%%%%
   for n=1:R
      [cS,cC(:,n),~,~,~,Rs22] = SM_correctionnew((S(:,n).'),abs(S(:,n).'));% peform phase de-ambiguity on shared SM S 
      sa = sort(abs(angle(cC(:,n))));
      if vox~=1
        al = sa(round(length(S)*(1-vox)));
        Cs(:,n) = max(eps+pi-al,(pi-abs(angle(cC(:,n)))));%+sigm
        in1(n,:) = find(Cs(:,n)==eps+pi-al);
      else in1(n,:) = 1:length(S); end
   end
   delta = S.*exp(-abs(S).^2/(2*sigma^2))/sigma^2;
   S1 = zeros(size(delta));
   for n=1:R
       S1(in1(n,:),n) = lambda*delta(in1(n,:),n);
   end
   S = S - S1*inv(M'*M); 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [U,Sigm,V]=svd(S,'econ'); S=U*V'; % add orthonormality constriaint
   M = X1*pinv(S.'); % update aggregating mixing matrix by rank-R least-square fit
   %%%%%%% shift-invariant rank-1 approximation %%%%%%%%%%%%%%
  for n=1:R
     mm = M(:,n);
     mm1 = reshape(mm,[size(X,2) size(X,3)]);
     xx = squeeze(mm1).';
     [C(:,n),B(:,n),~,T(:,n),phiz] = complexShiftedCP(xx,1);
     FACT{1} = C(:,n);          % subject intensity 
     FACT{2} = B(:,n);          % shared TC
     [m] = reconstructShiftCPc(FACT,T(:,n)); m=m';
     M(:,n) = m(:);            % new aggregating mixing matrix
  end
  sigma = sigma*0.995;%
  phi(i)=norm(X1-kr(C,B)*S.','fro')^2;  % error
  if phi(i)<1e-10; break; end
end