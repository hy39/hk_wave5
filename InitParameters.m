% InitParameters
% Read and initialize Parameters
% Function to return all required parameters and lookup arrays in the model
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [Params] = InitParameters()

%pa = readparams(filename);
%pa.R0 = 1.3;
pa.Tg = 6;          % Tg = infectious period until isolation
pa.Tgu = 7;          % Tgu = infectious period for undetected cases
pa.R = 5.08;
%pa.beta = (pa.R*(1-0.37))/pa.Tg;
%pa.betam = pa.beta;
%pa.betalate = pa.beta;
[pa.mob pa.dailymob]= get_dailymobility();
pa.beta = (pa.R*(1-0.254))/pa.Tg;      % 0.25 default, around Chinese new year
pa.betam = (pa.R*(1-0.278))/pa.Tg;      % 0.25 default, around Chinese new year
pa.betalate = (pa.R*(1-0.37))/pa.Tg;      % 0.37 default, transmission rate February 21 
pa.beta_mob = 1; % coefficient
%pa.betalate = (pa.R*(1-0.5))/pa.Tg;      % 0.5 default, transmission rate February 21 
pa.IFD50 = 200000;
pa.IFD50Boost = 1; %improvement of reporting rate
pa.det_rate_matrix = get_detection_rate(pa.IFD50);
pa.late_det_rate_matrix = get_detection_rate(pa.IFD50*pa.IFD50Boost);
pa.early_det_rate = 0.8;
pa.det_rate = 0.8;    %assumping 80% of cases are detected
%https://www.nytimes.com/interactive/2022/01/22/science/charting-omicron-infection.html
%https://www.cdc.gov/mmwr/volumes/70/wr/mm705152e3.htm

pa.latent = 3;    % 3 days latent period
pa.alpha = 0/365;   % immune decay
pa.wan = 0/30;      % period of transient immunity
pa.qr = 0/10;       % 1/10 allows 50% of infected cases (E and I) to be quarantined or isolated 
%pa.ages = [0 4; 5 9; 10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44; 45 49; 50 54; 55 59; 60 64; 65 69; 70 74; 75 79; 80 84; 85 120]; % [0 18; 19 39; 40 65; 65 100];
%pa.ages = [1 2; 2 3; 4 5; 6 7; 8 9; 10 11; 12 13]; % each age group define a lower to upper bound
pa.ages = [1 2];
pa.maxi = 7;        % antibody titre levels to strain A  
pa.maxj = 1;        % antibody titre levels to strain B (Not used)
pa.maxk = 1;        % antibody titre levels to strain C (Not used)
pa.maxa = size(pa.ages,1);        % number of age groups
pa.maxX = 1;        % max number of strains (Not used)
pa.N = 7.48E6;        % total population size
pa.seed = 2200;      % initial infected seeding
pa.trickle = 0;     % (Not used)
pa.age_flag = false; % True will create heterogeneous contact matrix M
pa.semi_mechanistic = 0; % (Not used)
pa.matM = make_M_simple(pa.age_flag);
%pa.matM = make_M_default(); % make_M_china();   % make contact mixture 
%pa.matM = make_M_china();
pa.AbB = 1.8;         % antibody boosting rate in each age group
%pa.AbB3 = 3.65;         % antibody 3rd dose mean value
pa.AbB2mu = 0.263;         % 2nd dose antibody boosting mean
pa.AbB2sigma = 0.475;      % 2nd dose antibody boosting standard deviation

%pa.AbB3mu = 1.6367;         % antibody 3rd dose mean value
%pa.AbB3sigma = 0.1889;    % antibody boosting standard deviation of lognormal
pa.AbB3mu = 1.566;         % antibody 3rd dose mean value for BNT
pa.AbB3sigma = 0.288;    % antibody boosting standard deviation of lognormal
%pa.AbB3cmu = 0.5663;         % antibody 3rd dose mean value for coronavac
%pa.AbB3csigma = 0.4036;    % antibody boosting standard deviation of lognormal
pa.AbB3cmu = 0.516;         % antibody 3rd dose mean value for coronavac
pa.AbB3csigma = 0.435;    % antibody boosting standard deviation of lognormal
pa.AbB3mmu = 1.48;         % antibody 3rd dose mean value for 2 coronavac mixed with BNT
pa.AbB3msigma = 0.281;    % antibody boosting standard deviation of lognormal
pa.AbBn = 50;         % antibody boosting rate in each age group
pa.AbB_mat = make_AbB_simple(pa.maxa); % assign antibody boosting for primary dose for each group
pa.AbB3_mat = make_AbB3_simple(pa.maxa); % assign antibody boosting for third dose for each group
pa.AbB3c_mat = make_AbB3c_simple(pa.maxa); % assign antibody boosting for third dose for each group
pa.AbBn_mat = make_AbBn_simple(pa.maxa); % assign antibody boosting for natural infection each age group
pa.arrf = make_f_simple();  % make infectivity matrix
pa.arrg = make_g_simple();  % make antibody boosting matrix
pa.arrg3 = make_g3_simple();  % make antibody boosting for third dose matrix
pa.arrg3c = make_g3c_simple();  % make antibody boosting for third dose (coronavac) matrix
pa.arrg3m = make_g3m_simple();  % make antibody boosting for a mixed (2coronavac + 1BNT) third dose matrix
pa.arrng = make_ng_simple();  % make antibody boosting for natural infection matrix
pa.ka = 1.77;   % old: 1.25 To define susceptibility
pa.kb = 3.38;  % old: 3.3561 To define susceptibility
pa.arrh = make_h_simple();  % make susceptibility matrix h(age,i titre ,j titre,k titre)
pa.vaccination_rate = 0.0005; % daily pre vaccination rate for the past 3 months before February 1  
pa.bntproportion = 0.61; %current scenario
v1 = 0.004; %current scenario 0.00285
v2 = 0.007;
pa.vaccination_rate3 = v1*pa.bntproportion; % daily vaccination rate for BNT 3rd Dose
pa.vaccination_rate3c = v1*(1-pa.bntproportion); % daily vaccination rate for Coronavac 3rd Dose 
[bnt sino tot bntprop] =  get_dailybooster(); %get daily booster rate
pa.dailybooster_bnt = bnt; %get daily booster rate
pa.dailybooster_sino = sino;
pa.dailybooster_tot = tot;
pa.bntproportion = bntprop;
pa.mixbntproportion = 0; %proportion of people take BNT as a mixed third dose

pa.temperature = get_dailytemperature(); %get daily temperature
pa.beta_temperature = 0; %set daily temperature coefficient
pa.relhumid = get_dailyrelhumidity(); %get daily relative humidity
pa.beta_relhumid = 0; %set daily relative humidity coefficient
pa.tqr = 82.525; %%10% patient can be traced and quarantiend
pa.ratr = 0; %% rate of rapid antigen test
pa.alpha1 = 0; % graduate change of rapid antigen test
pa.q = 0.114; %% 11.4% quarantiend people can contact others
[Slu Ilu CIlu ] = make_S_I_lookups();
pa.arrSlu = Slu;
pa.arrIlu = Ilu;
pa.arrCIlu = CIlu;
[ Slu Elu Ilu Rlu Qelu Qilu Hlu HNlu HRlu Vlu Vclu CIlu CHlu] = make_SEIRQH_lookups();
pa.arrSlu = Slu;
pa.arrElu = Elu;
pa.arrIlu = Ilu;
pa.arrRlu = Rlu;
pa.arrQelu = Qelu;   %Qe
pa.arrQilu = Qilu;  %Qi
pa.arrHlu = Hlu;
pa.arrHNlu = HNlu;  %From Necleic testing 
pa.arrHRlu = HRlu;  %From Rapid antigen testing 
pa.arrVlu = Vlu;
pa.arrVclu = Vclu;
pa.arrCIlu = CIlu;
pa.arrCHlu = CHlu;
%pa.beta = pa.R0./pa.Tg;
%pa.novars = max(pa.arrCIlu);
pa.novars = max(max(max(pa.arrCHlu)));
Params = pa;


function [ B ] = make_AbB_simple ( numAgeGroups )
    B = zeros(1,numAgeGroups);
    B(1:end) = pa.AbB;
end

function [ B ] = make_AbB3_simple ( numAgeGroups )
    B = zeros(1,numAgeGroups);
    B(1:end) = pa.AbB3mu;
end

function [ B ] = make_AbB3c_simple ( numAgeGroups )
    B = zeros(1,numAgeGroups);
    B(1:end) = pa.AbB3cmu;
end


function [ Bn ] = make_AbBn_simple ( numAgeGroups )
    Bn = zeros(1,numAgeGroups);
    Bn(1:end) = pa.AbBn;
end

function [ M ] = make_M_simple( age_flag )
%make_M_simple 
% Simple function to return a null mixing matrix
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
M = zeros(maxa,maxa);
for a=1:maxa
    for b=1:maxa
        if (age_flag)   %with age structure
            if a==b
                M(a,b)=1;
                %if a == 1
                % M(a,b)=5;
                %end
            else
                M(a,b)=0.2;
            end
        else            %homogeneous
            M(a,b)=1;
        end
        
    end
end
end


function [ M ] = make_M_default()
    maxa = pa.maxa;
    M = ones(maxa);
end


function [ M ] = make_M_china()
%make_M_simple 
% Simple function to return a null mixing matrix
% Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
M = ones(maxa,maxa);
%column is the infected source
%row is the susceptible target
top = 59364; % as the referece
tar = [237456 206068 196144 189552 149324 142848 141924 105712 101900 95412 94632 91244 86416 85104 82176 77252 77116 62636 62220 52920]; 

tar_k = tar(3)/top; %relative transportation frequency regarding to Toppest one
%M = [1 0;
%     0.001*tar_k 0];
%M = [
%10.92 0.6 0.8 0.4; 
%0.5 1.6 1 0.91;
%0.6 0.79 1.13 0.94;
%0.1 0.56 0.7 2.33;   
%    ];
end


function [ f ] = make_f_simple( input_args )
%make_f_simple Summary of this function goes here
%   Simple function to return a null infectivity array of the correct dimensons
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
f = zeros(maxa,maxi,maxj,maxk);
for a=1:maxa
    for i=1:maxi
        for j=1:maxj
            for k=1:maxk
                f(a,i,j,k)=1;
            end    
        end
    end
end
end


function [ f_] = make_ftable_simple( input_args )
%make_ftable_simple Summary of this function goes here
%   Simple function to return a null infectivity array of the correct dimensons
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C

f_ = zeros(maxa*maxi*maxj*maxk, 4);
cnt=0;
for a=1:maxa
    for i=1:maxi
        for j=1:maxj
            for k=1:maxk
                f_(cnt+1,1)=a;
                f_(cnt+1,2)=i;
                f_(cnt+1,3)=j;
                f_(cnt+1,4)=k;
                f_(cnt+1,5)=1;
                cnt = cnt + 1;
            end    
        end
    end
end
end


function [ n2 cnt] = make_g_simple( input_args )
%make_g_simple Summary of this function goes here
% Function to return the most simple immune state boost
% Written by Sean Yuan (hyuan@imperial.ac.uk) 

maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
n2 = zeros(maxa,maxl,maxm,maxn,maxl,maxm,maxn);

if pa.age_flag == 0
  pa.AbB_mat(1:end) = pa.AbB2mu; % 3rd dose immunization
end

%----------------- Based on Log normal --------------------%

%Add antibody boosting in each (age) group 
for i=1:pa.maxi
   for j=1:pa.maxa
       mu = pa.AbB2mu;
       boostpr_1 = lognpdf(i:pa.maxi,mu,pa.AbB2sigma);
       boostpr = zeros(1,length(boostpr_1));
       if sum(boostpr_1) > 0
            boostpr = boostpr_1./sum(boostpr_1);
       end
       age(j).Boosting(i).abl = boostpr; % Antibody boosting among each age group
   end
end
cnt = 0;
for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                         n2(a,l,m,n,l:pa.maxi,m,n)=age(a).Boosting(l).abl;
                    end
                end
            end
        end
    end %--end of age loop
end
end


function [ n3 cnt] = make_g3_simple( input_args )
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
n3 = zeros(maxa,maxl,maxm,maxn,maxl,maxm,maxn);

if pa.age_flag == 0
  pa.AbB3_mat(1:end) = pa.AbB3mu; % 3rd dose immunization
end

%----------------- Based on Log normal --------------------%

%Add antibody boosting in each (age) group 
for i=1:pa.maxi
   for j=1:pa.maxa
       %mu = log(pa.AbB3)-pa.AbB3sigma^2/2;
       mu = pa.AbB3mu;
       boostpr_1 = lognpdf(i:pa.maxi,mu,pa.AbB3sigma);
       boostpr = zeros(1,length(boostpr_1));
       if sum(boostpr_1) > 0
            boostpr = boostpr_1./sum(boostpr_1);
       end
       %boostpr(end) = 1-sum(boostpr(1:end-1)); % Add remaining to the last element
       age(j).Boosting(i).abl3 = boostpr; % Antibody boosting among each age group
   end
end
cnt = 0;
for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                         %1 loop:1 2 3 4 5 ...
                         %2 loop:  2 3 4 5 ...
                         %3 loop:    3 4 5 ...
                         %n3(a,l,m,n,l+[0:pa.maxi-l],m,n)=age(a).Boosting(l).abl3;
                         n3(a,l,m,n,l:pa.maxi,m,n)=age(a).Boosting(l).abl3;
                    end
                end
            end
        end
    end %--end of age loop
end
end


%% 3rd dose coronavac

function [ n3c cnt] = make_g3c_simple( input_args ) %Boosting for coronavac 3rd dose
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
n3c = zeros(maxa,maxl,maxm,maxn,maxl,maxm,maxn);

if pa.age_flag == 0
  pa.AbB3c_mat(1:end) = pa.AbB3cmu;
end

%----------------- Based on Log normal --------------------%

%Add antibody boosting in each (age) group 
for i=1:pa.maxi
   for j=1:pa.maxa
       %mu = log(pa.AbB3c)-pa.AbB3csigma^2/2; 
       mu = pa.AbB3cmu;
       boostpr = lognpdf(i:pa.maxi,mu,pa.AbB3csigma);
       boostpr = boostpr./sum(boostpr);
       %boostpr(end) = 1-sum(boostpr(1:end-1)); % Add remaining to the last element
       age(j).Boosting(i).abl3c = boostpr; % Antibody boosting among each age group
   end
end
cnt = 0;
for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                         %1 loop:1 2 3 4 5 ...
                         %2 loop:  2 3 4 5 ...
                         %3 loop:    3 4 5 ...
                         n3c(a,l,m,n,l+[0:pa.maxi-l],m,n)=age(a).Boosting(l).abl3c;
                         n3c(a,l,m,n,l:pa.maxi,m,n)=age(a).Boosting(l).abl3c;
                    end
                end
            end
        end
    end %--end of age loop
end
end


%% 3rd dose mixed BNT

function [ n3c cnt] = make_g3m_simple( input_args ) %Boosting for coronavac 3rd dose
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
n3c = zeros(maxa,maxl,maxm,maxn,maxl,maxm,maxn);

if pa.age_flag == 0
  pa.AbB3c_mat(1:end) = pa.AbB3cmu;
end

%----------------- Based on Log normal --------------------%

%Add antibody boosting in each (age) group 
for i=1:pa.maxi
   for j=1:pa.maxa
       mu = pa.AbB3mmu;
       boostpr = lognpdf(i:pa.maxi,mu,pa.AbB3msigma);
       boostpr = boostpr./sum(boostpr);
       age(j).Boosting(i).abl3c = boostpr; % Antibody boosting among each age group
   end
end
cnt = 0;
for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                         n3c(a,l,m,n,l+[0:pa.maxi-l],m,n)=age(a).Boosting(l).abl3c;
                         n3c(a,l,m,n,l:pa.maxi,m,n)=age(a).Boosting(l).abl3c;
                    end
                end
            end
        end
    end %--end of age loop
end
end

function [ ng cnt] = make_ng_simple( input_args )
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxl = pa.maxi; %Ab level to strain A
maxm = pa.maxj; %Ab level to strain B
maxn = pa.maxk; %Ab level to strain C
ng = zeros(maxa,maxl,maxm,maxn,maxl,maxm,maxn);

if pa.age_flag == 0
  pa.AbBn_mat(1:end) = pa.AbBn; % Natual Infection
end
%----------------- Based on Poisson --------------------%

%Add antibody boosting in each (age) group 
for i=1:pa.maxi
   for j=1:pa.maxa
       boostpr = ztpoisspdf(0:pa.maxi-i,pa.AbBn_mat(j));
       boostpr(end) = 1-sum(boostpr(1:end-1)); % Add remaining to the last element
       age(j).Boosting(i).abln = boostpr; % Antibody boosting among each age group
   end
end
cnt = 0;

for X=1:maxX
    for a=1:maxa
        for l=1:maxl
            for m=1:maxm
                for n=1:maxn
                    if X==1
                         ng(a,l,m,n,l+[0:pa.maxi-l],m,n)=age(a).Boosting(l).abln;
                    end
                    if X==2
                        if(m<=maxm-abl)
                            ng(a,l,m,n,l,m+abl,n)=1;
                            cnt = cnt + 1;
                        else
                            ng(a,l,m,n,l,m,n)=1;
                            cnt = cnt + 1;
                        end
                    end
                    if X==3
                        if(n<=maxn-abl)
                            ng(a,l,m,n,l,m,n+abl)=1;
                            cnt = cnt + 1;
                        else
                            ng(a,l,m,n,l,m,n)=1;
                            cnt = cnt + 1;
                        end
                    end
                end
            end
        end
    end %--end of age loop
end
end


function [ h ] = make_h_simple( input_args )
%make_h_simple Summary of this function goes here
% Takes a scalar set of parameters and returns an array with values for h.
% h(X,a,i,j,k) is the susceptibility of individuals to strain X dependent on their
% age and their antibody levels. Values are referenced to h(X,1,1,1,1)=1, i.e.
% those in the youngest age group with no detectable titres always have a relative
% suceptibility of 1
%   Written by Sean Yuan (hyuan@imperial.ac.uk) 
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
%h = zeros(maxX,maxa,maxi,maxj,maxk);
h = zeros(maxa,maxi,maxj,maxk);
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                        % h is linearly decreasing

                        %    h(a,i,j,k) = exp(-0.29*(i-1));
                        
                        %Change from 1.3 to 1.25
                        h(a,i,j,k)  = 1./(1+exp(pa.ka*[i-pa.kb])); %titre 1:25.6 represent titre level 3.3561
                        
                        %
                        if maxi==2
                           if i==1
                               h(a,i,j,k) = 1;
                           else
                               h(a,i,j,k) = 0;
                           end
                           
                        end
                end
            end
        end
    end

end

function [bnt sino total bntprop] = get_dailybooster()
  dat = load('HK_virus');
  bnt = dat.dailybnt/pa.N;
  sino = dat.dailysinovac/pa.N;
  total = dat.dailytotal/pa.N;
  bntprop = sum(dat.dailybnt(end-13:end))/(sum(dat.dailybnt(end-13:end))+sum(dat.dailysinovac(end-13:end)));
end

function [mob dailymob] = get_dailymobility()
  dat = load('HK_virus');
  mob = dat.mobility_7d;
  dailymob = dat.dailymobility;
end

function [tem] = get_dailytemperature()
  dat = load('HK_virus');
  tem = dat.temperature;
end

function [rh] = get_dailyrelhumidity()
  dat = load('HK_virus');
  rh = dat.relhumid;
end

function [ S I CI ] = make_S_I_lookups( input_args )
% make a lookup array for the linearized index of the isl model
% Can also be used to make list of column names for the solution table

maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
S = zeros(maxa,maxi,maxj,maxk);
I = zeros(maxa,maxi,maxj,maxk);
CI = zeros(maxa,maxi,maxj,maxk);
counter = 1;

%assign the number for susceptible
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    S(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

%assign the number for infected
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    I(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

%assign the number for accumulated infected
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    CI(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

end

%% Make SEIR lookup table
function [ S E I R Qe Qi H HN HR V Vc CI CH] = make_SEIRQH_lookups( input_args )
% make a lookup array for the linearized index of the isl model
% Can also be used to make list of column names for the solution table
maxX = pa.maxX; %number of strains
maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
S = zeros(maxa,maxi,maxj,maxk);
E = zeros(maxa,maxi,maxj,maxk);
I = zeros(maxa,maxi,maxj,maxk);
R = zeros(maxa,maxi,maxj,maxk);
Qe = zeros(maxa,maxi,maxj,maxk); %Quarantine latent 
HE = zeros(maxa,maxi,maxj,maxk); %Reporting latent
Qi = zeros(maxa,maxi,maxj,maxk); %Quarantine infectious
HN = zeros(maxa,maxi,maxj,maxk); %Reporting infectious after nucleic testing
HR = zeros(maxa,maxi,maxj,maxk); %Reporting infectious after rapid antigen testing
HB = zeros(maxa,maxi,maxj,maxk);
H = zeros(maxa,maxi,maxj,maxk);
Iu = zeros(maxa,maxi,maxj,maxk);
V = zeros(maxa,maxi,maxj,maxk);
Vc = zeros(maxa,maxi,maxj,maxk);
CI = zeros(maxa,maxi,maxj,maxk);
CH = zeros(maxa,maxi,maxj,maxk);
counter = 1;

%assign the number for susceptible
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    S(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for exposed (during latent period)
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    E(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for infected
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    I(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for infectious
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    R(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for quarantine or self isolation
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    Qe(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

 %assign the number for self isolation
        for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    Qi(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end


    %assign the number for hospital isolation or reported cases
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    H(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

    %assign the number for reported cases by necleic test 
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    HN(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
    %assign the number for reported cases by rapid antigenic test
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    HR(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    

%assign the number for vaccinated titres by BNT
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    V(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for vaccinated titres by coronavac
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    Vc(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
    
%assign the number for cumulated infected
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    CI(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end

%assign the number for accumulated infected
    for a=1:maxa
        for i=1:maxi
            for j=1:maxj
                for k=1:maxk
                    CH(a,i,j,k) = counter;
                    counter = counter + 1;
                end
            end
        end
    end
end

    


end