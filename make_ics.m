% make_ics
% Setep the initial conditions with everyone susceptible
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ yini age_arr pars] = make_ics( pars, arrSlu, arrIlu, arrCIlu)
%initial values are set to be 0
yini = zeros(1,pars.novars);
%age_arr = get_popweights_18_hk_mid_2019(pars); %get_popweights_5_uk_2003(pars);
age_arr = get_unipopwights(pars)*pars.N;
for a1=length(age_arr):pars.maxa-1
  age_arr = [age_arr age_arr(end)];
end
age_arr = age_arr*(pars.N/sum(age_arr));
pars.age_arr = age_arr;
seed_arr = (age_arr.*pars.seed)./pars.N; %seed uniformly
%seed_arr = [0.5 0.5 zeros(1,length(age_arr)-2)].*pars.seed;

for a=1:pars.maxa
    for i=1:pars.maxi
        for j=1:pars.maxj
            for k=1:pars.maxk
                yini(arrSlu(a,i,j,k)) = (age_arr(a) - seed_arr(a))*initial_prev_allAges(pars,a,i,j,k); %initial_prev_johnson_2009(pars,a,i,j,k);
            end
        end
    end
end

for a=1:pars.maxa
    for X=1:pars.maxX % this will cause total number not 1
        %X = 1;
        yini(arrIlu(a,1,1,1)) = seed_arr(a);
    end
end


function popweights = get_popweights_5_uk_2003(pa)
    popweights = zeros(1,pa.maxa);
    pop = [3009303, 8731641, 19889635, 12707659, 8453916]; 
    if pa.maxa~=5
      %disp(['max number of age classes is ' num2str(pa.maxa)]);%
      pop = pop(1:pa.maxa);
      popweights = pop./sum(pop);
    else
      pop = [3009303, 8731641, 19889635, 12707659, 8453916]; 
      popweights = pop./sum(pop);
    end
end


function popweights = get_popweights_18_hk_mid_2019(pa)
    popweights = zeros(1,pa.maxa);
    pop = [277600, 307400, 288600, 282500, 402100, 492200, 561400, 614700, 569800, 584600, 583000, 651100, 576700, 446600, 308200, 191500, 174900, 211200]; % refer to https://www.censtatd.gov.hk/hkstat/sub/sp150.jsp?tableID=002&ID=0&productType=8
    popweights = pop./sum(pop);
end

function popweights = get_unipopwights(pa)
    popweights = zeros(1,pa.maxa);
    w = 1/pa.maxa;
    pop = w*ones(1,pa.maxa); % can also refer to https://www.censtatd.gov.hk/hkstat/sub/sp150.jsp?tableID=002&ID=0&productType=8
    popweights = pop./sum(pop);
end


function init_prev = initial_prev_johnson_2009(pa,a,i,j,k)
        init_prev = [];
        if(a==1)
            prev = 25/58;
            if i==2 
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==2)
            prev = 0.25;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==3)
            prev = 0.18;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==4)
            prev = 0.2;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        elseif(a==5)
            prev = 0.30;
            if i==2
                init_prev = prev;
            elseif i==1
                init_prev = 1-prev;
            else
                init_prev = 0;
            end
        else
            error('Problem in initial_prev_johnson_2009');
        end
end


function init_prev = initial_prev_allAges(pa,a,i,j,k) %same initial prevalence for all age groups
    init_prev = [];
    prev = 0;
    if i==2
        init_prev = prev;
    elseif i==1
        init_prev = 1-prev;
    else
        init_prev = 0;
    end
end


end