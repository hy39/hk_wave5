% Return the matrix of Force of Infection including the effects of temperature
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)

function [ matFOI ] = getMOFQTemp( pars,y )

maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
beta = pars.beta;   % transmission rate for unquarantined
betaq = pars.betaq; % transmission rate for quarantined
arrf = pars.arrf; %infectivity array
arrIlu = pars.arrIlu; %I lookup
arrQilu = pars.arrQilu; %Q lookup
R = pars.R;
Tg = pars.Tg;
q = pars.q;           % infectivity of quarantined people
matM = pars.matM; %M[a b], b infects a
lambda = 0; %force of infection
matFOI = zeros(1,maxa);
%calculate for of infection
    for a=1:maxa %a represents the age groups in susceptible individuals
        tmpval_transmisibility = 0;
        for b=1:maxa  %b represents the age groups in infected individuals
            tmpval = 0;
            for i=1:maxi
                for j=1:maxj
                    for k=1:maxk
                        tmpval = tmpval + beta*arrf(b,i,j,k)*y(arrIlu(b,i,j,k));
                        tmpval = tmpval + betaq*y(arrQilu(b,i,j,k)); % few isolated patients can still infect others
                    end
                end
            end
            tmpval_transmisibility = tmpval_transmisibility + tmpval*matM(a,b);  
            
        end
        matFOI(a) = tmpval_transmisibility;
    end
return;
end

