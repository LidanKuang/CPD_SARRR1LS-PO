function [A,B,C,T,phiz,err] = complexShiftedCP(X,R)
Tol1=1e-8;  MaxIt1=500; Ninit=1; 
% Matrix Unfoldings of X
[I,J,K]=size(X); nrmodes = ndims(X);
X1=tens2mat(X,1);  % X1 is IKxJ (or smaller dimensions if compression was done)
X2=tens2mat(X,2);  % X2 is JIxK (idem)
if nrmodes==3, X3=tens2mat(X,3); end % X3 is KJxI (idem)
noc = R; N = [I J K];  Ns = N(2); 
SST=norm(X(:))^2;
sigma_sq=SST/((1+10^(0/10))*(prod(N)));
Lambda=ones(1,noc)*eps*sigma_sq;
f = j*2*pi*[0:Ns-1]/Ns;
xx2 = matricizing(X,2).';
Xf = fft(xx2,[],2); Nf = [I J K];
Xf = unmatricizing(Xf.',2,[N(1), Nf(2), N(3:end)]);
T = zeros(I,R);
TauW = ones(noc,1)*[-10 10];
TauWMatrix = generateTauWMatrix(TauW,N(2));
for ninit=1:Ninit   % Try Ninit different initializations
             
    % INITIALIZATION
    % Case 1: none of the input matrices A,B,C is given
    if K>1
      [A,B,C]=cp3_init(X,R,'random');  % force random init 
    else
      [A,B]=cp3_init(X,R,'random');  % force random init 
    end

    A1=A;A2=A; if K>1, C1=C;C2=C; end  % useful for Line Search

    % LOOP for alternating updates
    Sf = fft(B.',[],2); Sf1 = Sf; Sf2 = Sf;
    krprt = ones(1,noc);krkrt = ones(noc);
    if K>1
        krprt = krprod(C,krprt); krkrt = krkrt.*(C'*C);
    end
    krSf = krprod(krprt,Sf.').'; krf = krprod(ones(size(krprt,1),1),f.').';
    for k = 1:I   
       Resf(k,:) = Xf(k,:)-A(k,:)*(krSf.*exp(T(k,:)'*krf));
    end
    phi = norm(Resf,'fro')^2;
    stop=0;
    phiz=[];
    it1=0;
    
% if strcmp(lsearch,'elsc')==1
    Pa2=fliplr(pascal(4));
    Pa=zeros(4);
    for n=1:4
        Pa(n,n:end)=diag(Pa2,n-1).';
    end     
    Ja=repmat([1i^3,1i^2,1i,1],4,1);

    Pb2=fliplr(pascal(4));
    Pb=zeros(4);
    for n=1:4
        Pb(n,n:end)=diag(Pb2,n-1).';
    end 
    Jb=toeplitz([1,0,0,0],[1,1i,1i^2,1i^3]);
    
  while stop==0 
     it1=it1+1;
    %% Line Search  %%
     dA=A1-A2;  % search direction for A
     dSf=Sf1-Sf2;  % search direction for B
      if K==1, 
          dC = zeros(1,noc); C2 = ones(1,noc);
          [A,Sf,C] = cp3_lsearch(A2,Sf2,C2,dA,dSf,dC,X1.','elsc',it1,Xf,T,f,Pa,Ja,Pb,Jb); 
          phi_old = caphi(A,Sf,C,X1.',T,f);
      else
          dC=C1-C2;
          [A,Sf,C] = cp3_lsearch1(A2,Sf2,C2,dA,dSf,dC,X3,'elsc',it1,Xf,T,f,Pa,Ja,Pb,Jb); 
          phi_old = caphi(A,Sf,C,X3,T,f);
      end
      B = (ifft(Sf,[],2)).';
   %% update B %%
   krprt = ones(1,noc);
   krkrt = ones(noc);
   if nrmodes>2
       krprt = krprod(C,krprt);
       krkrt = krkrt.*(C'*C);
   end
   if nrmodes>2
       Xfp = krprt'*matricizing(Xf, [3:nrmodes]);
       Xfp = unmatricizing(Xfp,3,[Nf(1:2), noc]);
       Xtp = krprt'*matricizing(X, [3:nrmodes]);
       Xtp = unmatricizing(Xtp,3,[N(1:2), noc]);       
       krprtt = sparse(1:noc,1:noc,[ones(1,noc)]);   
   else
       Xfp = Xf;
       Xtp = X;
       krprtt = ones(1,noc);   
   end
   Sf = fft(B.',[],2);
   [Sf]=updateS(Sf,Xfp,A,T,f,krprtt,krkrt,Nf,sigma_sq*Lambda);
   B = (ifft(Sf,[],2)).';
    
   %% update T %%
   P.Sf = krprod(krprt,Sf.').';
   P.f = krprod(ones(size(krprt,1),1),f.').';
   P.w = ones(1,length(f));
   P.w(2:end-1) = 2;
   P.w = krprod(ones(size(krprt,1),1),P.w.').';
   P.sizeX2 = size(X,2);       
   P.At = A;
   P.Xft = matricizing(Xf,1);
   for m = 1:size(Xf,1)
       if K>1
          [T(m,:)] = estTimeAutCorc(P.Xft(m,:),X3,A(m,:),B,C,P.Sf,P.f,T(m,:),Nf,TauWMatrix);
       else
          [T(m,:),A(m,:)] = estTimeAutCorc1(P.Xft(m,:),A(m,:),B,ones(1,noc),P.Sf,P.f,T(m,:),Nf,TauWMatrix);
       end
   end
   
    %%  update A %%
      ind = 1:N(1):prod(N(1:2));
       for k=1:I
           Sft = Sf.*exp(T(k,:)'*f);
           S = ifft(Sft,[],2); 
           if K>1
               A(k,:) = (inv((C'*C).*(conj(S)*S.'))*(kr(C,S.')'*X3(:,k))).';
               Ft(ind+k-1,:) = krprod(S.',A(k,:));
           else
               A(k,:) = (inv((krprt'*krprt).*(conj(S)*S.'))*(kr(krprt,S.')'*X1(k,:).')).'; 
           end
       end
     if K>1, krkrt = Ft'*Ft; end
   %%  update C %%
   if K>1
      krpr = Ft;
      krkr = krkrt;
      Xt = matricizing(X,3);
      Xtkrpr=krpr'*X2;
      C = (inv(krkr)*Xtkrpr).';
   end
    %%%%%%%%%%%%%%%%%%%%%%
      % Normalization to avoid overflow
       normA=sqrt(sum(A.*conj(A)));  % Frobenius norm of each column of A
       normB=sqrt(sum(B.*conj(B)));
       if K>1, normC=sqrt(sum(C.*conj(C))); else normC = 1; end
       prod_norm=normA.*normB.*normC;
       if K>1, Scale_mat=diag(prod_norm.^(1/3)); else  Scale_mat=diag(prod_norm.^(1/2)); end
       % equal repartition of power of each rank-1 tensor over the 3 vectors:
       A=A*diag(1./normA)*Scale_mat;
       B=B*diag(1./normB)*Scale_mat; 
       if K>1, C=C*diag(1./normC)*Scale_mat; end
       % Calculate the new fit to the model and decide to stop or not
       Sf = fft(B.',[],2);
       if K>1, phi = caphi(A,Sf,C,X3,T,f); else phi = caphi(A,Sf,krprt,X1.',T,f); end
       % Define a stop criterion
       err(it1) = ((phi-phi_old)/phi_old);
       if  (it1==MaxIt1) || (phi<Tol1) || (abs(err(it1))<Tol1)
          stop=1;
       end
       phiz(it1) = phi;
       % Store matrices to prepare next Line Search step
        A2=A1;Sf2=Sf1;
        A1=A;Sf1=Sf;
        if K>1, C1=C;C2=C1; end
        if (phi_old-phi)/phi_old<0
            A = A1; Sf = Sf1; phiz(it1) = phi_old;
        end
    end    % Algorithm has converged for this initialization
    
% figure,plot(phiz(1:it1-1))
end   % END of ALL init