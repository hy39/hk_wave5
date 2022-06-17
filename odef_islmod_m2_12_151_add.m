% odef_islmod 
% ODE model for ploting the impact of each intervention
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ ydot ] = odef_islmod_m2_12_151_add( t,y,pars )
ydot = zeros(1,length(y));
maxX = pars.maxX; %number of strains
maxa = pars.maxa; %number of age groups
maxi = pars.maxi; %Ab level to strain A
maxj = pars.maxj; %Ab level to strain B
maxk = pars.maxk; %Ab level to strain C
beta = pars.beta;
arrf = pars.arrf'; %infectivity array
arrh = pars.arrh'; %susceptibility array
%alpha = pars.alpha;  %immunity decay rate
R = pars.R;          %intrinsic R0
Tg = pars.Tg;        %infectious period
latent = pars.latent; %latent period 
%lookup tables
arrSlu = pars.arrSlu;
arrElu = pars.arrElu; 
arrIlu = pars.arrIlu;
arrRlu = pars.arrRlu;
arrQelu = pars.arrQelu;
arrQilu = pars.arrQilu;
arrHlu = pars.arrHlu;
arrHNlu = pars.arrHNlu;
arrHRlu = pars.arrHRlu;
arrVlu = pars.arrVlu;
arrVclu = pars.arrVclu;
arrCIlu = pars.arrCIlu;
dailybooster_bnt = pars.dailybooster_bnt; %get daily booster rate
dailybooster_sino = pars.dailybooster_sino;
dailybooster_tot = pars.dailybooster_tot;
bntproportion = pars.bntproportion;
temperature = pars.temperature;
beta_temperature = pars.beta_temperature; %coeeficient of temperature
relhumid = pars.relhumid;
beta_relhumid = pars.beta_relhumid; %coeeficient of RH
ratr = pars.ratr; % Rapid antigen testing rate
alpha1 = pars.alpha1; % Gradual change of rapid antigen testing rate
q = pars.q; % contact rate of quarantined
%mob = pars.dailymob % don't use daily mobility
mob = pars.mob.*pars.beta_mob;
qr = 1/pars.tqr; %quarantine rate
wan = pars.wan;
TestTime = 100;
MassTestFlag = false;
T_ntr = 0.1; % average waiting time for testing
ntr = 1/T_ntr;

%% Scenario 1: Basic Social Distancing + Rapid Antigenic Test
%% Mobility as the average between 01 - 09 February
%% No T2, T3, or RAT 
if pars.scenario == 1
ratr = 0; % no RAT  
% use daily mobility
if t < 9 % Until 9 February 
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mean(mob(1:9)));
end

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 2: Basic Social Distancing + Tightened Regulations
%% "No RAT"
if pars.scenario == 2
ratr = 0; % no RAT  
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 3: Current NPIs 
%% Considering weather effects
if pars.scenario == 3
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 4: Current NPIs 
%% Weather fixed at day 20
if pars.scenario == 4
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
if t < length(dailybooster_bnt)
  tid = floor(t) + 1;
  pars.vaccination_rate3 = dailybooster_bnt(tid);
  pars.vaccination_rate3c = dailybooster_sino(tid);
else
  pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
  pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < 20
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(20))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(20))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 5: Current NPIs 
%% With no vaccine booster after Feburary
if pars.scenario == 5
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
pars.vaccination_rate3 = 0;
pars.vaccination_rate3c = 0;
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 6: Current NPIs 
%% With 2x vaccine booster after Feburary
if pars.scenario == 6
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid)*2;
pars.vaccination_rate3c = dailybooster_sino(tid)*2;
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end))*2;
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end))*2;
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end


%% Scenario 7: Current NPIs 
%% Temperature without the cooldest two days 
if pars.scenario == 7
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3);
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
end
if t > 19 && t < 21
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(1:28))+5-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(1:28))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 8: Current NPIs 
%% Without considering weather effects
%% Without booster since February
%% Without social distancing tightening
if pars.scenario == 8
ratr = 0; % no RAT  
% use daily mobility
if t < 17
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(17));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
pars.vaccination_rate3 = 0;
pars.vaccination_rate3c = 0;

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end

end


%%==
%% Scenario 9: Current NPIs 
%% Temperature fixed after February
if pars.scenario == 9
% use daily mobility
if t < length(mob)
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mob(end));
end
pars.betaq = q*pars.beta;

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < 28
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(1:28))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(1:28))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end

%% Scenario 10: Basic Social Distancing + Rapid Antigenic Test
%% Mobility as the average between 10 - 23 February
%% No T3 or RAT
if pars.scenario == 10
ratr = 0; % no RAT  
% use daily mobility
if t < 23 % Until 23 February 
  tid = floor(t) + 1;
  pars.beta = (R/Tg)*(1-mob(tid));
else
  pars.beta = (R/Tg)*(1-mean(mob(10:23)));
end

%% Get daily booster rate
if t < length(dailybooster_bnt)
tid = floor(t) + 1;
pars.vaccination_rate3 = dailybooster_bnt(tid);
pars.vaccination_rate3c = dailybooster_sino(tid);
else
    pars.vaccination_rate3 = mean(dailybooster_bnt(end-13:end));
    pars.vaccination_rate3c = mean(dailybooster_sino(end-13:end));
end

%use daily temperature
%temperature = movmean(temperature,3); 
if t < length(temperature)
  tid = floor(t) + 1;
  pars.beta = pars.beta*exp(beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
  pars.beta = pars.beta*exp(beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
  pars.betaq = q*pars.beta;
else
  pars.beta = pars.beta*exp(beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28))));
  pars.beta = pars.beta*exp(beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28))));
  pars.betaq = q*pars.beta;
end
end


%First calculate the force of infection
matFOI = getMOFQTemp( pars,y ); %matrix of lambda for each age group
%Second setup the variables for immune boosting
%InBoost % goes into S compartment
%OutBoost % leaves from I compartment
[InBoost OutBoost ImmDecayS ImmDecayI ] = getImmBoost( pars,y ); %Immune boosting and decay

%% set detection rate
% https://service.hket.com/search/result?dis=basic&keyword=%E6%96%B0%E5%86%A0%E8%82%BA%E7%82%8E
%det_rate = pars.early_det_rate;
Iprev = round(sum(y(pars.arrIlu)));
maxIprev = 5000000;
if (Iprev<1)
    Iprev = 1;
end
if (Iprev < maxIprev)
    det_rate = pars.det_rate_matrix(Iprev);
else
    det_rate = pars.det_rate_matrix(maxIprev);
end

%% set reporting delay
if Iprev > 30000
    T_ntr2 = T_ntr + (1.2-T_ntr)*(1-exp(-(Iprev-30000)/300000));
    ntr = 1/T_ntr2;
end

if t > 26
    det_rate2 = pars.early_det_rate;
    t_point = 31;
    det_max = det_rate+(1-det_rate)*det_rate2;
    if t < t_point
         det_rate_tmp = det_max - (det_max-det_rate)*exp(-alpha1*(t-26));
    else
         det_rate_tmp = (det_max-det_rate)*exp(-alpha1*(t-t_point)) + det_rate;
    end
    det_rate = det_rate_tmp;
end

%% set RAT rate
t_rat = 25;
if t<=t_rat
    ratr = 0;
end

%% When 5 millions of people vaccinated by 3rd dose, stop it
totalvaccinecov = sum(y(arrVlu(1,:,1,1))) + sum(y(arrVclu(1,:,1,1)));
if totalvaccinecov > pars.N*0.65 
     pars.vaccination_rate3 = 0; % Check the current rate
     pars.vaccination_rate3c = 0; % Check the current rate
end

[InVBoost3 OutVBoost3] = getVacBoost3( pars,y ); %Immune boosting and decay
[InVBoost3c OutVBoost3c] = getVacBoost3c( pars,y ); %Immune boosting and decay

InBoost = InBoost';
OutBoost = OutBoost';
InVBoost3 = InVBoost3';
OutVBoost3 = OutVBoost3';
InVBoost3c = InVBoost3c';
OutVBoost3c = OutVBoost3c';

                    Total_NewInfected = 0;                  
                    %calculate for infection
                    NewInf = [];
                    DecayInf = [];
                    RecoverInf = [];
                    X = 1;
                    RecoverInf = 0;

                    %% Classic SEIR
                    NewInf = matFOI*(y(arrSlu)./pars.N).*arrh; %h represent susceptibility
                    Vaccinated = InVBoost3 - OutVBoost3;
                    VaccinatedC = InVBoost3c - OutVBoost3c; % Coronavac vaccinated people
                    % Test
                    %Vaccinated = 0;
                    %VaccinatedC = 0;
                    %Vaccinated2 = InVBoost(a,i,j,k) - OutVBoost(a,i,j,k);
                    NewVaccinated = InVBoost3;
                    NewVaccinatedc = InVBoost3c;

                    AfterLatent = 1/latent*y(arrElu);
                    ReSusc = wan*y(arrRlu);
                    RecoverInffromI = (1./Tg).*OutBoost;
                    RecoverInftoR = ((1-det_rate)/Tg).*InBoost; %undetection
                    RecoverInftoH = (det_rate/Tg).*InBoost;
                    
                    
                    %RecoverQtoH = (1./Tg)*y(arrQilu);
                    RecoverQtoR = (1./Tg)*y(arrQilu);
                    QuarantinefromE = qr*y(arrElu);
                    QuarantinefromI = qr*y(arrIlu);
                    QetoQi = (1./latent)*y(arrQelu);
                    ydot(arrSlu) = -NewInf + ReSusc + Vaccinated + VaccinatedC;
                    ydot(arrElu) = NewInf - AfterLatent - QuarantinefromE;
                    ydot(arrIlu) = AfterLatent - RecoverInffromI - QuarantinefromI;
                    ydot(arrRlu) = -ReSusc + RecoverInftoR + RecoverQtoR;
                    ydot(arrQelu) = QuarantinefromE - QetoQi;
                    ydot(arrQilu) = QuarantinefromI + QetoQi - RecoverQtoR;
                    %When infectious cases are quarantined, they are tested and reported immediately 
                    ydot(arrHlu) = RecoverInftoH + QuarantinefromI + (1./latent)*y(arrQelu) - ntr*y(arrHlu); % + RecoverQtoH;
                    ydot(arrHNlu) = ntr*y(arrHlu);
                    ydot(arrVlu) = NewVaccinated;
                    ydot(arrVclu) = NewVaccinatedc; 
                    ydot(arrCIlu) = NewInf;

                    if (t>t_rat)
                        Rate = ratr;
                        Reported = Rate*y(arrIlu);
                        ydot(arrIlu) = AfterLatent - RecoverInffromI - QuarantinefromI - Reported;
                        ydot(arrQilu) = QetoQi + QuarantinefromI - RecoverQtoR + Reported;
                        
                        % option 1 Contact-traced report late
                        ydot(arrHlu) = RecoverInftoH + QuarantinefromI + QetoQi - ntr*y(arrHlu);
                        ydot(arrHNlu) = ntr*y(arrHlu);
                        ydot(arrHRlu) = Reported;
                        %% option 2 Contact-traced report immediately
                        %ydot(arrHlu) = RecoverInftoH; 
                        %ydot(arrHRlu) = Reported + QuarantinefromI + QetoQi;
                    else
                        %Rate = ratr;
                        %Reported = Rate*y(arrIlu);
                        %ydot(arrIlu) = AfterLatent - RecoverInffromI - QuarantinefromI - Reported;
                        %ydot(arrQilu) = QetoQi + QuarantinefromI - RecoverQtoR + Reported;
                        %ydot(arrHlu) = RecoverInftoH + QuarantinefromI + QetoQi + Reported; 
                    end


ydot = ydot';
return;
end

