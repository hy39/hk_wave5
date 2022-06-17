% calculate effective R for SEIR model
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk) 
function [Re] = calRt_fullpar_age(beta, tau, Tqr, q, Tc, Thos, ratr, Si, arrh)

N = 7480000; 

Tqr = Tqr;
    
if exist('Undet_CR')
else
    Undet_CR = 0;
end    
CR = 1;

%Transmission Matrix
ItoE = sum(beta.*(Si).*arrh(1,:,1,1));
QitoE = sum(q*beta.*(Si).*arrh(1,:,1,1));

T =   [    0         ItoE        0             QitoE;
           0         0           0             0;
           0         0           0             0; 
           0         0           0             0];


gamma2 = 1/Thos;

S = [-(1/tau+1/Tqr)   0             		0               0;
     1/tau            -(1/Tc+1/Tqr+ratr)  	0               0;
     1/Tqr            0              		-1/tau          0;
     0                1/Tqr+ratr     		1/tau           -gamma2];

K = -T*inv(S);
ek = eig(K);
Re = ek(1);
end




