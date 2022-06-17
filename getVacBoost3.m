% Antibody response after booster vaccination from BNT
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ InVBoost OutVBoost ] = getVacBoost3( pars,y )

maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
arrg3 = pars.arrg3; %immune boosting array
arrg3m = arrg3;
mprot = 0;
if isfield(pars, 'mixbntproportion')
   mprot = pars.mixbntproportion; %proportion of 3rd BNT that is mixed 
   arrg3m = pars.arrg3m; %immune boosting array mixed doses
end
vac_rate = pars.vaccination_rate3; %daily vaccination rate
age_arr = pars.age_arr; %age group population size
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
InBoostPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); % goes into S compartment, per strain
InVBoost = zeros(maxa,maxi,maxj,maxk); % goes into S compartment
OutVBoost = zeros(maxa,maxi,maxj,maxk); % leaves from I compartment
ImmDecaySPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); %


%calculate for of ImmBoosting into S and out from I 
%%%Vac 應該只要針對titre 1,2,3,4
%%%



    for a=1:maxa
        for i=1:maxi %titre i is the original titre
            for j=1:maxj
                for k=1:maxk
                    %要計算Sai的注射率 = (Na x r) x Sai/(Sai=1 + Sai=2 + ...)
                    %只有titre level 1, 2的人注射第三劑
                    if i <= 2
                        if maxi == 1
                          largest_titre_tobe_vaccinated = 1;
                        else
                          largest_titre_tobe_vaccinated = 2;
                        end
                        totalSforAforVac = sum(y(arrSlu(a,1:largest_titre_tobe_vaccinated,j,k))); %計算group當中titre level 1,2的Susceptible人數
                        %totalSforA = sum(y(arrSlu(a,1:end,j,k))); %計算group當中的Susceptible人數
                        %應該計算每個age group人數
                        totalNforA = pars.age_arr(a);
                        ratioforVac = totalNforA/totalSforAforVac; %Na/(Sai=1 + Sai=2 + ...)
                        if y(arrSlu(a,largest_titre_tobe_vaccinated,j,k)) - vac_rate*ratioforVac*y(arrSlu(a,largest_titre_tobe_vaccinated,j,k)) > 0
                            OutVBoost(a,i,j,k) = vac_rate*ratioforVac*y(arrSlu(a,i,j,k)); % Third dose
                        else
                            if i == 1 && largest_titre_tobe_vaccinated == 2
                               S1forAforVac = y(arrSlu(a,1,j,k));
                               totalNforA = pars.age_arr(a);
                               ratioforVac = totalNforA/S1forAforVac;
                               OutVBoost(a,i,j,k) = vac_rate*ratioforVac*y(arrSlu(a,i,j,k));
                               %disp("no one can be vaccinated");
                               InVBoost = zeros(maxa,maxi,maxj,maxk); % goes into S compartment
                               OutVBoost = zeros(maxa,maxi,maxj,maxk); % leaves from I compartment
                               return; 
                            end  
                        end
                    end 
                        for l=1:largest_titre_tobe_vaccinated  %titre l is the source titre
                            for m=1:maxj
                                for n=1:maxk
                                  % boosting from l to i 
                                    if mprot == 0 
                                        InVBoost(a,i,j,k) = InVBoost(a,i,j,k) + vac_rate*ratioforVac*y(arrSlu(a,l,m,n))*arrg3(a,l,m,n,i,j,k); % Third dose
                                    else
                                        InVBoost(a,i,j,k) = InVBoost(a,i,j,k) + mprot*(vac_rate*ratioforVac*y(arrSlu(a,l,m,n))*arrg3(a,l,m,n,i,j,k)) + (1-mprot)*(vac_rate*ratioforVac*y(arrSlu(a,l,m,n))*arrg3m(a,l,m,n,i,j,k)); % Third dose
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

