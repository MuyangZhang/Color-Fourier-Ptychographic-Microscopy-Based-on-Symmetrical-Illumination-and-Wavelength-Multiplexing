function [ O,P ] = GDUpdate_Multiplication_rank_rm(O,P,dpsi,Omax,cen,Ps,alpha,beta)
%GDUPDATE_MULTIPLICATION update estimate of O and P according to gradient
%descent method, where psi = O*P
%   Inputs:
%   O0: object estimate, n1xn2
%   P0: pupil function estimate: m1xm2
%   psi: update estimate field estimate
%   psi0: previous field estimate
%   cen: location of pupil function
%   alpha: gradient descent step size for O
%   betta: gradient descent step size for P
%   Ps: support constraint for P0, e.g. spatially confined probe or
%   objective with known NA
%   iters: # of iterations to run on updates
%
% last modified by Lei Tian, lei_tian@alum.mit.edu, 3/1/2014

% size of O
[mO, nO, ~] = size(O);
[mP, nP, ~] = size(P);
No = [mO, nO]; No = No(:);
% size of P, Np<=No
Np = [mP, nP]; Np = Np(:);
% cen0 = round((No+1)/2); %%%
r = size(dpsi,3);

dO = zeros(No(1),No(2),3);
dP = 0;
sumP = zeros(No(1),No(2),3);
sumO = 0;
for m = 1:r
    % operator to put P at proper location at the O plane
    n1 = cen(:,m)-floor(Np/2);
    n2 = n1+Np-1;
    
    % operator to crop region of O from proper location at the O plane
    downsamp = @(x) x(n1(1):n2(1),n1(2):n2(2));
    
    dO0 = abs(P(:,:,m)).*conj(P(:,:,m)).*dpsi(:,:,m);
    dO(n1(1):n2(1),n1(2):n2(2),m) = dO0;%dO(n1(1):n2(1),n1(2):n2(2),m)+
    
    O1 = downsamp(O(:,:,m));
    dP = abs(O1).*conj(O1).*dpsi(:,:,m);%dP + 
    
    sumP(n1(1):n2(1),n1(2):n2(2),m) = abs(P(:,:,m)).^2;%sumP(n1(1):n2(1),n1(2):n2(2),m)+
    sumO = abs(O1).^2;%sumO+
    P(:,:,m) = P(:,:,m)+1/Omax*dP./(sumO+beta).*Ps(:,:,m);%conj(O1).*dpsi(:,:,m)/max(max(abs(O1).^2));%
    O(:,:,m) = O(:,:,m)+1/max(max(abs(P(:,:,m))))*dO(:,:,m)./(sumP(:,:,m)+alpha);
%     O(n1(1):n2(1),n1(2):n2(2),m) = O(n1(1):n2(1),n1(2):n2(2),m)+conj(P(:,:,m)).*dpsi(:,:,m)/max(max(abs(P(:,:,m)).^2));%

end



end

