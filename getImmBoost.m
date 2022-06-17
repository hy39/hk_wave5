% getImmBoost
% Calculate antibody responses after natural infection
% Return the ImmBoosting 1) into S and 2) out from I
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk) 
function [ InBoost OutBoost ImmDecayS ImmDecayI] = getImmBoost( pars,y )

maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
arrng = pars.arrng; %immune boosting array
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
InBoostPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); % goes into S compartment, per strain
InBoost = zeros(maxa,maxi,maxj,maxk); % goes into S compartment
OutBoost = zeros(maxa,maxi,maxj,maxk); % leaves from I compartment
ImmDecaySPerStrain = zeros(maxX,maxa,maxi,maxj,maxk); %
ImmDecayS = zeros(maxa,maxi,maxj,maxk); %
ImmDecayI = zeros(maxX,maxa,maxi,maxj,maxk); %

%calculate for of ImmBoosting into S and out from I 

    for a=1:maxa
        for i=1:maxi %titre i is the original titre
            for j=1:maxj
                for k=1:maxk
                    % Immune boosting
                     InBoost(a,i,j,k) = 0;
                     OutBoost(a,i,j,k) = 0;
                        for l=1:maxi  %titre l is the target titre
                            for m=1:maxj
                                for n=1:maxk
                                  InBoost(a,i,j,k) = InBoost(a,i,j,k) + y(arrIlu(a,l,m,n))*arrng(a,l,m,n,i,j,k);
                                  OutBoost(a,i,j,k) = OutBoost(a,i,j,k) + y(arrIlu(a,i,j,k))*arrng(a,i,j,k,l,m,n); %total number of cases leaving from I(a,i)   
                                end
                            end
                        end   
                end
            end    
        end
    end
    
    
end

