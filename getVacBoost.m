% Antibody response after full immunisation
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ InVBoost OutVBoost ] = getVacBoost( pars,y )

maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
arrg = pars.arrg; %immune boosting array
vac_rate = pars.vaccination_rate*pars.bntproportion; %daily vaccination rate for BNT 2nd dose
age_arr = pars.age_arr; %age group population size
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
InBoostPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); % goes into S compartment, per strain
InVBoost = zeros(maxa,maxi,maxj,maxk); % goes into S compartment
OutVBoost = zeros(maxa,maxi,maxj,maxk); % leaves from I compartment
ImmDecaySPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); %


%calculate for of ImmBoosting into S and out from I 

    for a=1:maxa
        for i=1:maxi %titre i is the original titre
            for j=1:maxj
                for k=1:maxk
                    % Immune boosting
                    % InBoost(a,i,j,k) = 0;
                    if i == 1
                    %OutVBoost(a,i,j,k) = OutVBoost(a,i,j,k) + vac_rate*y(arrSlu(a,i,j,k)); % Third dose
                    OutVBoost(a,i,j,k) = OutVBoost(a,i,j,k) + vac_rate*age_arr(a); % Prime dose 
                    end
                        %OutBoost(X,a,i,j,k) = 0;
                        for l=1:maxi  %titre l is the target titre
                            for m=1:maxj
                                for n=1:maxk
                                  % boosting from l to i 
                                  if l == 1
                                    InVBoost(a,i,j,k) = InVBoost(a,i,j,k) + vac_rate*age_arr(a)*arrg(a,l,m,n,i,j,k);
                                  end
                                end
                            end
                        end
                        %InBoost(a,i,j,k);
                            
                end
            end    
        end
    end
    
    
end

