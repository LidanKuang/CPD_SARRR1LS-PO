function [A,Sf,C] = cp3_lsearch(A,Sf,C,dA,dSf,dC,X,lsearch,it,Xf,T,f,Pa,Ja,Pb,Jb) 
%CP3_LSEARCH Line search for CANDECOMP/PARAFAC order 3.
%   New loading matrices A (IxR), B (JxR), C(KxR) of cp3 are computed
%   from their previous values and from the search directions dA, dB, dC,
%   as follows:
%  
%      A <- A + mu * dA
%      B <- B + mu * dB
%      C <- C + mu * dC
%
%   Line Search for cp3 can for instance be used in gradient-based
%   algorithms or with Alternating Least Squares (ALS), in order to speed 
%   up convergence.
%
%   For instance, if Line Search is used with ALS, the search directions
%   may be defined outside this function by dA=A1-A, dB=B1-B, dC=C1-C and
%   the input matrices by A=A1, B=B1 and C=C1 where (A1,B1,C1) denote
%   conditional least squares estimates at the current ALS iteration it,
%   whereas (A,B,C) denote old estimates, such that
% 
%      A <- A + mu*(A1-A)
%      B <- B + mu*(B1-B)
%      C <- C + mu*(C1-C)
%
%   This means that, in the context of ALS, setting mu=1 annihilates the
%   line search, since A, B and C will be set to their current values A1,
%   B1 and C1.
%
%   For instance, if Line Search is used with gradient-based algorithms,
%   the search directions dA, dB, dC, could be (- gradient) of the cost 
%   function wrt A, B and C, respectively, in a gradient descent 
%   algorithm, or other more sophisticated search directions like in a
%   conjugate gradient algorithm.
%
%   Several choices of mu are possible:
%
%      - if lsearch = 'none' then no interpolation is done (in the sense of 
%        lsearch for als), i.e. we enforce mu=1.
%      - if lsearch = 'lsh' (Line Search proposed by Harshman), then mu is 
%        fixed to 1.25 and if the interpolated matrices do not decrease 
%        the cost function, then mu=1 is enforced.
%      - if lsearch = 'lsb' (Line Search proposed by Bro), then mu is fixed
%        to it^(1/3) and if the interpolated matrices do not decrease 
%        the cost function, then mu=1 is enforced.
%      - if lsearch = 'elsr' (Exact Line Search with real step), then we 
%        seek for the optimal real-valued step mu that minimizes the cost 
%        function. This amounts to minimizing a polynomial of degree 6 in
%        mu. It can also be used for complex-valued data.
%      - if lsearch = 'elsc' (Exact Line Search with complex step), then we
%        seek for the optimal complex-valued step mu that minimizes the
%        cost function. If the data are real-valued, the step has to be
%        real-valued, and in this case elsc is replaced by elsr.
%
%   INPUTS: 
%   
%      - A,B,C: estimates of the loading matrices A,B,C at it-1
%      - dA,dB,dC: search directions
%      - X: KJxI matrix unfolding of the observed tensor X 
%      - lsearch: = 'none', 'lsh','lsb', 'elsr' or 'elsc'
%      - it is the iteration step number       
%      - Pa,Ja,Pb,Jb: matrices used only by elsc in order to alternate
%        between updates of real(mu) and imag(mu) in an efficient way.
%       These matrices are generated outside this function as follows:
%           Pa2=fliplr(pascal(4));Pa=zeros(4);
%    		for n=1:4;Pa(n,n:end)=diag(Pa2,n-1).';end     
%    		Ja=repmat([1i^3,1i^2,1i,1],4,1);
%    		Pb2=fliplr(pascal(4));Pb=zeros(4);
%           for n=1:4;Pb(n,n:end)=diag(Pb2,n-1).';end 
%           Jb=toeplitz([1,0,0,0],[1,1i,1i^2,1i^3]);
%
%   OUTPUTS: 
%
%      - Updated loading matrices A,B,C

%   Copyright 2010
%   Version: 09-07-10
%   Authors: Dimitri Nion (dimitri.nion@gmail.com),
%
%   References:
%   [1] R.A. Harshman, "Foundations of the PARAFAC procedure: Model and 
%       Conditions for an explanatory Multi-mode Factor Analysis", UCLA 
%       Working Papers in Phonetics, vol.16, pp 1-84, 1970.
%   [2] R. Bro, "Multi-Way Analysis in the Food Industry: Models, 
%       Algorithms, and Applications", PhD. dissertation, University of 
%       Amsterdam, 1998.
%   [3] M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A 
%       Novel Method to Accelerate  PARAFAC", SIAM Journal Matrix Anal. and
%       Appl. (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008.
%   [4] D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for 
%       Complex-Valued Tensor Decompositions. Application in DS-CDMA", 
%       Signal Processing, vol. 88, issue 3, pp. 749-755, March 2008.

it_start=3; % Line Search will start after it_start iterations (at least 2)
if isreal(X) && strcmp(lsearch,'elsc')==1; lsearch='elsr';end

%   The cp3 least-squares cost function is  
%   phi = norm (X - kr(C,B)*A.' , 'fro')^2, 
%   where kr is the Khatri-Rao product and X is the mode-1 KJ by I matrix 
%   unfolding of the IxJxK tensor that is decomposed. Note that, in the 
%   row indexing of X, the index k=1,..,K has to vary more slowly than 
%   j=1,..,K, since the cost function involves kr(C,B) and not kr(B,C).
%
%   The crucial point is to choose a good step size mu such that
%      phinew = norm (X - kr((C+mu*dC),(B+mu*dC))*(A+mu*dA).' , 'fro')^2,   
%   is smaller than phi.

% lsearch='none', i.e., standard ALS            
if strcmp(lsearch,'none')==1
     mu=1; 
     A=A+mu*dA;
     B=B+mu*dB;
     C=C+mu*dC;

% lsearch='lsh'
elseif strcmp(lsearch,'lsh')==1
    if it<it_start
        mu=1;
        A=A+mu*dA;
        B=B+mu*dB;
        C=C+mu*dC;   
    else
        % Compute phi with mu=1
        A=A+dA;
        B=B+dB;
        C=C+dC;
%         phi=norm(X-kr(C,B)*A.','fro');
        phi = sqrt(caphi(Xf,A,B,C,T,f));
        % Compute phi with mu=1.25
        mu=1.25;
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
%         phim=norm(X-kr(Cm,Bm)*Am.','fro');
        phim = sqrt(caphi(Xf,Am,Bm,Cm,T,f));
        % Accept or Reject Interpolation                   
        if phim < phi     % accept
            A=Am;
            B=Bm;
            C=Cm;              
        end
    end   
% lsearch='lsb'         
elseif strcmp(lsearch,'lsb')==1
    if it<it_start
        mu=1;
        A=A+mu*dA;
        B=B+mu*dB;
        C=C+mu*dC;  
    else
        % Compute phi with mu=1
        A=A+dA;
        B=B+dB;
        C=C+dC;
%         phi=norm(X-kr(C,B)*A.','fro');
        phi = sqrt(caphi(Xf,A,B,C,T,f));
        % Compute phi with mu=it^(1/3)
        mu=it^(1/3);
        Am=A+mu*dA;
        Bm=B+mu*dB;
        Cm=C+mu*dC;
%         phim=norm(X-kr(Cm,Bm)*Am.','fro');
        phim = sqrt(caphi(Xf,Am,Bm,Cm,T,f));
        % Accept or Reject Interpolation                   
        if phim < phi     % accept
            A=Am;
            B=Bm;
            C=Cm;              
        end
    end
% lsearch='elsr', compute optimal real-valued mu
elseif strcmp(lsearch,'elsr')==1
    if it<it_start
        mu=1;  
    else
        KdCdB=kr(dC,dB);
        KdCB=kr(dC,B);
        KCdB=kr(C,dB);
        KCB=kr(C,B);
        Mat3=KdCdB*dA.';
        Mat2=KdCdB*A.' + (KdCB+KCdB)*dA.';
        Mat1=KCB*dA.' + (KdCB+KCdB)*A.';
        Mat0=KCB*A.'-X;
        M=[Mat3(:) Mat2(:) Mat1(:) Mat0(:)]; 
        H_mat=real(M'*M);  
        % Now we define the coefficients of the 6th order polynomial
        d6=H_mat(1,1);
        d5=2*H_mat(1,2);
        d4=2*H_mat(1,3)+H_mat(2,2);
        d3=2*(H_mat(1,4)+H_mat(2,3));
        d2=2*H_mat(2,4)+H_mat(3,3);
        d1=2*H_mat(3,4);
        d0=H_mat(4,4);
        pol=[d6 d5 d4 d3 d2 d1 d0];
        pol_der=[6*d6 5*d5 4*d4 3*d3 2*d2 d1];
        sqrts=roots(pol_der);
        sqrts=sqrts(imag(sqrts)==0); % real roots
        sqrts=[sqrts;1];   
        % Choice of optimal mu
        extremum=polyval(pol,sqrts);
        mu=sqrts(find(extremum==min(extremum),1));
        mu=mu(1);
    end
    A=A+mu*dA;
    B=B+mu*dB;
    C=C+mu*dC;   
% lsearch='elsc'
% Alternate between updates of real and imaginary part of mu
elseif strcmp(lsearch,'elsc')==1
    if it<it_start
        mu=1; 
    else
        Tolelsc=1e-4;
        Niterelsc=50;
%         Sf = fft(B.',[],2); dSf = fft(dB.',[],2); 
        for k=1:size(A,1)
           Sft = Sf.*exp(T(k,:)'*f); dSft = dSf.*exp(T(k,:)'*f);
           S = ifft(Sft,[],2); dS = ifft(dSft,[],2); 
           KdCdB = kr(dC,dS.'); %KdCdB=kr(dC,dB);
           KdCB = kr(dC,S.'); %KdCB=kr(dC,B);
           KCdB = kr(C,dS.'); %KCdB=kr(C,dB);
           KCB = kr(C,S.');   %KCB=kr(C,B);
           Mat3(:,k)=KdCdB*dA(k,:).';         %Mat3=KdCdB*dA.';
           Mat2(:,k)=KdCdB*A(k,:).' + (KdCB+KCdB)*dA(k,:).';%Mat2=KdCdB*A.' + (KdCB+KCdB)*dA.';
           Mat1(:,k)=KCB*dA(k,:).' + (KdCB+KCdB)*A(k,:).';  %Mat1=KCB*dA.' + (KdCB+KCdB)*A.';
           Mat0(:,k)=KCB*A(k,:).'-X(:,k);                   %Mat0=KCB*A.'-X;
           M(:,:,k)=[Mat3(:,k) Mat2(:,k) Mat1(:,k) Mat0(:,k)]; 
           H_mat(:,:,k)=M(:,:,k)'*M(:,:,k);  
        % Initialization
        mu=1;
        a=1;  % real part of step
        u=[mu^3 mu^2 mu 1].';
        phi_new(k)=abs(u'*H_mat(:,:,k)*u);
        end
%         phi_new=caphi1(X2,A+mu*dA,B+mu*dB,C+mu*dC,T,f);
        phi0=sum(phi_new);   % initial value of the cost function (with mu=1)
        phi_diff=sum(phi_new);
        it_in=0;
        % Alternate between updates of real part a and imag part b of mu 
        while (phi_diff > Tolelsc) && (it_in < Niterelsc)
            it_in=it_in+1;
            phi_old=phi_new;
            % Update imag part b      
            Da=Pa.*Ja.*toeplitz([1,0,0,0],[1 a a^2 a^3]); 
            for k=1:size(A,1)
            Ha_mat=real(Da'*H_mat(:,:,k)*Da);
            % Now we define the coefficients of the 6th order
            % polynomial in b
            d6(k)=Ha_mat(1,1);
            d5(k)=Ha_mat(1,2)+Ha_mat(2,1);
            d4(k)=Ha_mat(1,3)+Ha_mat(2,2)+Ha_mat(3,1);
            d3(k)=Ha_mat(1,4)+Ha_mat(2,3)+Ha_mat(4,1)+Ha_mat(3,2);
            d2(k)=Ha_mat(2,4)+Ha_mat(3,3)+Ha_mat(4,2);
            d1(k)=Ha_mat(3,4)+Ha_mat(4,3);
            d0(k)=Ha_mat(4,4);
            end
            pol=[sum(d6) sum(d5) sum(d4) sum(d3) sum(d2) sum(d1) sum(d0)];
            pol_der=[6*sum(d6) 5*sum(d5) 4*sum(d4) 3*sum(d3) 2*sum(d2) sum(d1)];
            sqrts=roots(pol_der);
            sqrts=sqrts(imag(sqrts)==0);
            sqrts=[sqrts;1];      
            % Choice of optimal b
            extremum=polyval(pol,sqrts);
            b=sqrts(find(extremum==min(extremum),1));
            b=b(1);

            % Update real part a                   
            Db=Pb.*Jb.*toeplitz([1,0,0,0],[1 b b^2 b^3]);
            for k=1:size(A,1)
            Hb_mat=real(Db'*H_mat(:,:,k)*Db);
            % Now we define the coefficients of the 6th order
            % polynomial in a
            d6(k)=Hb_mat(1,1);
            d5(k)=Hb_mat(1,2)+Hb_mat(2,1);
            d4(k)=Hb_mat(1,3)+Hb_mat(2,2)+Hb_mat(3,1);
            d3(k)=Hb_mat(1,4)+Hb_mat(2,3)+Hb_mat(4,1)+Hb_mat(3,2);
            d2(k)=Hb_mat(2,4)+Hb_mat(3,3)+Hb_mat(4,2);
            d1(k)=Hb_mat(3,4)+Hb_mat(4,3);
            d0(k)=Hb_mat(4,4);
            end
            pol=[sum(d6) sum(d5) sum(d4) sum(d3) sum(d2) sum(d1) sum(d0)];
            pol_der=[6*sum(d6) 5*sum(d5) 4*sum(d4) 3*sum(d3) 2*sum(d2) sum(d1)];
            sqrts=roots(pol_der);
            sqrts=sqrts(imag(sqrts)==0);
            sqrts=[sqrts;1];   % we add 1 in the set of possible values    
            % Choice of optimal a
            extremum=polyval(pol,sqrts);
            a=sqrts(find(extremum==min(extremum),1));
            a=a(1);

            % update mu
            mu=a+1i*b;
            u=[mu^3 mu^2 mu 1].';
            for k=1:size(A,1)
               phi_new(k)=abs(u'*H_mat(:,:,k)*u);
            end
%             phi_new=caphi1(X,A+mu*dA,Sf+mu*dSf,C+mu*dC,T,f);
            phi_diff=abs(sum(phi_new)-sum(phi_old));         
        end
        if phi_new>phi0; mu=1;end
    end
    % Interpolate
    A=A+mu*dA;
    Sf=Sf+mu*dSf;
    C=C+mu*dC;
end
end