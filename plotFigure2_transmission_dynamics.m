% plotFigure2_transmission_dynamics
% Plot the transmission dynamics and cumulative incidence of the fifth wave
% Figure 2CD in the paper
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [] = plotFigure2_transmission_dynamics()

%figure('Renderer', 'painters', 'Position', [20 500 1080 400]); %for Mac 
figure('Renderer', 'painters', 'Position', [20 500 1200 450]); %for Windows

load('mcmc/fullmod_mcmc_output');
ps = PosteriorSamples; %posterior distributions
sys_par = sys_par;%McMC system parameters 
pars = Par_stat.par;%model parameters

%% model parameters
names = ps.Properties.VariableNames
%generate 100 random number 
no = height(ps);%the number of row in table ps
a = sys_par.burnIn;
b = no;
sample = 20;
r = round((b-a).*rand(sample,1)) + a;%size: 100*1

y_pre = pars.y0;
totaldays = 120;
times = 0:1:totaldays;
sample_time = 1; % orginal 40d
pars.sample_time = sample_time;
times_sim = 0:0.25:times(end);

R_list = ps.R(r);
tqr_list = ps.Tqr(r);
q_list = ps.q(r);
IFD50_list = ps.IFD50(r);
early_det_rate_list = ps.early_det_rate(r);
ka_list = ps.ka(r);
kb_list = ps.kb(r);
beta_temperature_list = ps.beta_temperature(r);
beta_relhumid_list = ps.beta_relhumid(r);
%pars.temperature = movmean(pars.temperature,3); 
ratr_list = ps.ratr(r);
alpha1_list = ps.alpha1(r);
beta_mob_list = ps.beta_mob(r);
%% Create prediction intervals
for  i=1:length(r)
    i
pars.R = R_list(i);
pars.tqr = tqr_list(i);
pars.q = q_list(i);
pars.IFD50 = IFD50_list(i);
pars.det_rate_matrix = get_detection_rate(pars.IFD50);
pars.early_det_rate = early_det_rate_list(i);
pars.ka = ka_list(i);
pars.kb = kb_list(i);
pars.beta_temperature = beta_temperature_list(i);
pars.beta_relhumid = beta_relhumid_list(i);
pars.ratr = ratr_list(i);
pars.alpha1 = alpha1_list(i);
pars.beta_mob = beta_mob_list(i);
%% Preseting
% set susceptibility
pars.arrh = make_h(pars);

[t y] = ode23(@(t,y)odef_islmod_m1_12_151(t,y,pars), times_sim, y_pre);
tt = t(1:4:end);
y = y(1:4:end,:);
local_cum(i,:) = sum(y(:,pars.arrCIlu(1,:,1,1)),2); %True infections
local_total(i,:) = diff(local_cum(i,:));
local_cum_H(i,:) = sum(y(:,pars.arrHNlu(1,:,1,1)),2); %Reported infections
local_total_H(i,:) = diff(local_cum_H(i,:));
local_cum_HR(i,:) = sum(y(:,pars.arrHRlu(1,:,1,1)),2); %Reported infections
local_total_HR(i,:) = diff(local_cum_HR(i,:));

local_cum_R(i,:) = sum(y(:,pars.arrRlu(1,:,1,1)),2); %Undetected
local_daily_R(i,:) = diff(local_cum_R(i,:)); %Daily undetected 
local_cum_CH(i,:) = sum(y(:,pars.arrCHlu(1,:,1,1)),2); %Detected
local_daily_CH(i,:) = diff(local_cum_CH(i,:)); %Daily Detected

    mob = pars.mob;
    temperature = pars.temperature;
    relhumid = pars.relhumid;
    q = pars.q;
    ratr = pars.ratr;   
    for tid = 1:tt(end)
    ratr = pars.ratr;    
    p.beta = 0;
    Thos = pars.Tg;
    if tid <= length(mob)
        p.beta = (pars.R/pars.Tg)*(1-mob(tid));
    else
        p.beta = (pars.R/pars.Tg)*(1-mob(end));
    end
    t_rat = 25;
    if tid <= t_rat
      ratr = 0;
    else
      ratr = pars.ratr;   
    end

    if tid < length(temperature)
        p.beta = p.beta*exp(pars.beta_temperature*(temperature(tid)-mean(temperature(1:28)))); % mean temperature in February
        p.beta = p.beta*exp(pars.beta_relhumid*(relhumid(tid)-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
        p.betaq = q*p.beta;
    else
        p.beta = p.beta*exp(pars.beta_temperature*(mean(temperature(end-6:end))-mean(temperature(1:28)))); % mean temperature in February
        p.beta = p.beta*exp(pars.beta_relhumid*(mean(relhumid(end-6:end))-mean(relhumid(1:28)))); % mean relative humidity (RH) in February is reference RH
        p.betaq = q*p.beta;
    end   
    Re = calRt_fullpar_age(p.beta, pars.latent, pars.tqr, pars.q, pars.Tg, Thos, ratr, y(tid,pars.arrSlu(1,:,1,1))./pars.N, pars.arrh);
    Re_p3(i,tid) = Re;
    end
end

tt2 = tt'
tt2(end) = []; % the last day of new cases is blank
n = length(tt2); 
t0 = 1;

subplot(1,2,1);
hold on;
local_total_mean = mean(local_total);
local_total_bound = prctile(local_total,[2.5 97.5]);
XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
YY = [local_total_bound(1,t0), local_total_bound(2,t0:end), local_total_bound(1,n), fliplr(local_total_bound(1,t0:end))];
h1 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
b1 = plot(tt2(t0:end),local_total_mean(t0:end),'b','linewidth',1.5);

dat = load('HK_virus');
b3 = bar(0:length(dat.cases)-1,[dat.cases dat.cases_rat] ,'stacked');

local_total_H_mean = mean(local_total_H+local_total_HR);
local_total_H_bound = prctile(local_total_H+local_total_HR,[2.5 97.5]);
XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
YY = [local_total_H_bound(1,t0), local_total_H_bound(2,t0:end), local_total_H_bound(1,n), fliplr(local_total_H_bound(1,t0:end))];
h2 = fill(XX,YY,[237/255 218/255 218/255],'Linestyle','none');
b2 = plot(tt2(t0:end),local_total_H_mean(t0:end),'r','linewidth',1.5);

local_total_H_PCR_mean = mean(local_total_H);
b14 = plot(tt2(t0:end),local_total_H_PCR_mean(t0:end),'r--','linewidth',1.5);

legend([b1 b2 b14 b3 b3],'Model predicted true infections','Model predicted reported cases','Model predicted reported cases by PCR', 'Reported cases by PCR', 'Reported cases by RAT')

set(h2,'facealpha',.5);
date0 = ([0 14 28 42 58 89]);
date = ({'01/02','15/02','01/03','15/03','31/03','01/05'});
set(gca,'xtick',date0,'XTickLabel',date);
xtickangle(90);
ylabel('Daily number of infections'); 
set(gca,'FontSize',14);
xlim([0 59]);
ylim([0 3.5E05]);

% Plot the cumulative dynamics
subplot(1,2,2);
hold on;
local_cum(:,end)=[];
local_cum = local_cum./pars.N*100;
local_cum_mean = mean(local_cum);
local_cum_bound = prctile(local_cum,[2.5 97.5]);
XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
h1 = fill(XX,YY,[218/255 218/255 227/255],'Linestyle','none');
b1 = plot(tt2(t0:end),local_cum_mean(t0:end),'b','linewidth',1.5);

cumcases = cumsum(dat.cases);
cumcases_rat = cumsum(dat.cases_rat);
b3 = bar(0:length(dat.cases)-1, [(cumcases./pars.N)*100 (cumcases_rat./pars.N)*100],'stacked');


local_cum_H(:,1) = [];
local_cum_HR(:,1) = [];
local_cum_Htot = (local_cum_H+local_cum_HR)./pars.N*100;
local_cum_H_mean = mean(local_cum_Htot);
local_cum_H_bound = prctile(local_cum_Htot,[2.5 97.5]);
XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
YY = [local_cum_H_bound(1,t0), local_cum_H_bound(2,t0:end), local_cum_H_bound(1,n), fliplr(local_cum_H_bound(1,t0:end))];
h2 = fill(XX,YY,[227/255 218/255 218/255],'Linestyle','none');
b2 = plot(tt2(t0:end),local_cum_H_mean(t0:end),'r','linewidth',1.5);
set(h2,'facealpha',.5);

local_cum_H_mean = mean(local_cum_H)./pars.N*100;
b4 = plot(tt2(t0:end),local_cum_H_mean(t0:end),'r--','linewidth',1.5);
legend([b1 b2 b4 b3 b3],'Model predicted true infections','Model predicted reported cases', 'Model predicted reported cases by PCR', 'Reported cases by PCR', 'Reported cases by RAT');

date0 = ([0 14 28 42 58 89]);
date = ({'01/02','15/02','01/03','15/03','31/03','01/05'});
set(gca,'xtick',date0,'XTickLabel',date);
xtickangle(90);
ylabel('Cumulative infections (% of the population)'); 
set(gca,'FontSize',14);
xlim([0 59]);
end

