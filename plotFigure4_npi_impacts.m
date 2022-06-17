% plotFigure4_npi_impacts
% Plot the impacts of individual NPIs and weather conditions
% Figure 4 in the paper
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [] = plotFigure4_npi_impacts()

load('mcmc/fullmod_mcmc_output');
ps = PosteriorSamples; %posterior distributions
sys_par = sys_par;%McMC system parameters 
pars = Par_stat.par;%model parameters

%figure('Renderer', 'painters', 'Position', [20 500 1050 420]); %for Mac
figure('Renderer', 'painters', 'Position', [20 500 1200 450]); %for Windows
hold on;

% draw samples
no = height(ps);%the number of row in table ps
a = sys_par.burnIn;
b = no;
sample = 20;
r = round((b-a).*rand(sample,1)) + a;%size: 100*1
ci = zeros(5,sample);

%% model parameters
for subfig = 1:2
if subfig == 1  % T1-3+RAT(3), with T1(1), T1,2(10), T1-3(2), T1-3+RAT-without large temperature drop or humidity increase(7)
    Sc = [3 1 10 2]; 
end
if subfig == 2  % T1-3+RAT(3), with T1(1), T1,2(10), T1-3(2), T1-3+RAT-without large temperature drop or humidity increase(7)
    Sc = [3 7 9]; 
end
names = ps.Properties.VariableNames
endaxis = 59;

y_pre = pars.y0;
totaldays = 120;
times = 0:1:totaldays;
sample_time = 1; % orginal 40d
pars.sample_time = sample_time;
times_sim = 0:0.25:times(end);
pars.scenario = 1;

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
%pars.temperature = movmean(pars.temperature,3); 
% 1 no interventions
% 2 with intervention but without the effect of weather
% 3 with intervention and with the effect of weather
b1 = [];
b2 = [];
b3 = [];
b4 = [];
b5 = [];
b6 = [];
b21 = [];
b22 = [];
b23 = [];
b24 = [];
b25 = [];
b26 = [];
%Sc = [1 2 3]; 
for sid=1:length(Sc)
%% 建立prediction intervals
s = Sc(sid);
local_cum = []; %True infections
local_total = [];
local_cum_H = []; %Reported infections
local_total_H = [];
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
%% Preseting
% set susceptibility
pars.arrh = make_h(pars);
pars.scenario = s;
pars2 = pars;
%[t y1] = ode23(@(t,y)odef_islmod_m1_12_14(t,y,pars1), times_sim, y_pre);
[t y2] = ode23(@(t,y)odef_islmod_m2_12_151_add(t,y,pars2), times_sim, y_pre);
tt = t(1:4:end);
y = y2(1:4:end,:);
%transform y into proportion
%y = y./pars.N;
%detected_local(1:end-1);
local_cum(i,:) = sum(y(:,pars.arrCIlu(1,:,1,1)),2); %True infections
local_total(i,:) = diff(local_cum(i,:));
local_cum_H(i,:) = sum(y(:,pars.arrHlu(1,:,1,1)),2); %Reported infections
local_total_H(i,:) = diff(local_cum_H(i,:));
end

tt2 = tt';
tt2(end) = []; % the last day of new cases is blank
n = length(tt2); 
t0 = 1;


% Plot the cumulative dynamics
local_cum = local_cum./pars.N.*100;
local_cum(:,end)=[];
ci(pars.scenario,:) = local_cum(:,59);


subplot(1,2,subfig);
hold on;
if pars.scenario == 1
local_cum_mean = mean(local_cum);
local_cum_bound = prctile(local_cum,[2.5 97.5]);
t0 = 9;
XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
h21 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
b21 = plot(tt2(t0:end),local_cum_mean(t0:end),'b:','linewidth',2);
set(h21,'facealpha',.5);
elseif pars.scenario == 2
    t1 = 25;
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    XX = [tt2(t1),  tt2(t1:end),  tt2(end),  fliplr(tt2(t1:end))];
    YY = [local_cum_bound(1,t1), local_cum_bound(2,t1:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t1:end))];
    h22 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b22 = plot(tt2(t1:end),local_cum_mean(t1:end),'b:','linewidth',2);
    set(h22,'facealpha',.5);
elseif pars.scenario == 3
    if subfig == 2
        %t0 = 26;
        t0 = 1;
    end
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h23 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b23 = plot(tt2(t0:end),local_cum_mean(t0:end),'b','linewidth',1.8);
    %b23 = plot(tt2(t0:end),local_cum_mean(t0:end),':','Color',[90/255 237/255 90/255],'linewidth',1.8);
    set(h23,'facealpha',.5);
    %legend([b21 b22 b23],'Initial tightening (up to 17 Feb)','Strengthened tightening','Current NPIs: strengthened tightening + Rapid antigen test');
    t0 = 1;
elseif pars.scenario == 4
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h24 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b24 = plot(tt2(t0:end),local_cum_mean(t0:end),'b','linewidth',1.8);
    set(h24,'facealpha',.5);
    %legend([b4 b3],'Daily weather maintained after incidence peak','Daily weather');
elseif pars.scenario == 5
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h25 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b25 = plot(tt2(t0:end),local_cum_mean(t0:end),'b','linewidth',1.8);
    set(h25,'facealpha',.5);
elseif pars.scenario == 6
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h26 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b26 = plot(tt2(t0:end),local_cum_mean(t0:end),'b--','linewidth',1.8);
    set(h26,'facealpha',.5);
    legend([b25 b26 b23],'Without booster since February','Booster rate (2x)','Current booster rate');
elseif pars.scenario == 7
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    t0 = 20;
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h27 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b27 = plot(tt2(t0:end),local_cum_mean(t0:end),'b--','linewidth',1.8);
    set(h27,'facealpha',.5);
    t0 = 1;
    %legend([b25 b26 b23],'Without booster since February','Booster rate (2x)','Current booster rate');
elseif pars.scenario == 9
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    t0 = 28;
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h29 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b29 = plot(tt2(t0:end),local_cum_mean(t0:end),'b:','linewidth',2);
    set(h29,'facealpha',.5);
    t0 = 1;
elseif pars.scenario == 10
    local_cum_mean = mean(local_cum);
    local_cum_bound = prctile(local_cum,[2.5 97.5]);
    t0 = 23;
    XX = [tt2(t0),  tt2(t0:end),  tt2(end),  fliplr(tt2(t0:end))];
    YY = [local_cum_bound(1,t0), local_cum_bound(2,t0:end), local_cum_bound(1,n), fliplr(local_cum_bound(1,t0:end))];
    h210 = fill(XX,YY,[218/255 218/255 237/255],'Linestyle','none');
    b210 = plot(tt2(t0:end),local_cum_mean(t0:end),'b:','linewidth',2);
    set(h210,'facealpha',.5);
end


date0 = ([0 28 59 89]);
date = ({'Feb','Mar','Apr','May'});
xlim([0 endaxis]);
set(gca,'xtick',[0 14 28 42 58],'XTickLabel',{'01/02','15/02','01/03','15/03','31/03'});
xtickangle(90);
ylim([0 65]);
ylabel({'Proportio of cumulative infections (%)'});
set(gca,'FontSize',14);
set(gca,'ytick',[0 20 40 60 65],'YTickLabel',{'0','20','40','60',''});
end

end
ci; % T1-3+RAT(3), with T1(1), T1,2(10), T1-3(2), T1-3+RAT-without large temperature drop or humidity increase(7)
T1_ci = ci(1,:);
T2_ci = ci(10,1);
T3_ci = ci(2,:);
T3_RAT_ci = ci(3,:);
T3_RAT_temp_ci = ci(7,:); 

T2_effect = (ci(10,:)-ci(1,:))./ci(1,:); %Without T2,3
T3_effect = (ci(3,:)-ci(1,:))./ci(1,:);  %Without T3
RAT_effect = (ci(3,:)-ci(2,:))./ci(3,:); %Without RAT
Vac_effect = (ci(3,:)-ci(5,:))./ci(3,:); 

end


