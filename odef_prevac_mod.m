% odef_prevac_mod
% ODE model to obain disease dynamics
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ ydot ] = odef_prevac_mod( t,y,pars )
ydot = zeros(1,length(y));
maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
beta = pars.beta;
arrf = pars.arrf; %infectivity array
arrh = pars.arrh; %susceptibility array
%alpha = pars.alpha; %immunity decay rate
Tg = pars.Tg; %infectious period
%lookup tables
arrSlu = pars.arrSlu; 
arrIlu = pars.arrIlu; 
arrCIlu = pars.arrCIlu;
low_titre = 1;

%[InBoost OutBoost ImmDecayS ImmDecayI ] = getImmBoost( pars,y ); %Immune boosting and decay
[InVBoost OutVBoost] = getVacBoost( pars,y ); %Immune boosting and decay

if t > (90-76) && t <= (90-26) %76 days befefore January 25 (a week before February 1)
    pars.vaccination_rate3 = 0.001*0.6; % Check the current rate
    pars.vaccination_rate3c = 0.001*0.4; % Check the current rate
    [InVBoost3 OutVBoost3] = getVacBoost3( pars,y ); %Immune boosting and decay
    [InVBoost3c OutVBoost3c] = getVacBoost3c( pars,y ); %Immune boosting and decay
end
if t > 90-26
    [InVBoost3 OutVBoost3] = getVacBoost3( pars,y ); %Immune boosting and decay
    [InVBoost3c OutVBoost3c] = getVacBoost3c( pars,y ); %Immune boosting and decay
end


InVBoost3 = zeros(size(InVBoost));
OutVBoost3 = zeros(size(OutVBoost));
InVBoost3c = zeros(size(InVBoost));
OutVBoost3c = zeros(size(OutVBoost));


%getMOF should return a two by two matrix.
    Total_NewInfected = 0;
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk                    


                    X = 1;
                    ydot(arrIlu(a,i,j,k)) = 0;
                    ydot(arrCIlu(a,i,j,k)) = 0;
                    ydot(arrSlu(a,i,j,k)) = InVBoost(a,i,j,k) - OutVBoost(a,i,j,k) + InVBoost3(a,i,j,k) - OutVBoost3(a,i,j,k) + InVBoost3c(a,i,j,k) - OutVBoost3c(a,i,j,k);

                end
                
            end
        end
    end
  


ydot = ydot';
return;
end

