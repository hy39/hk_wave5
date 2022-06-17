% plotFigure3_population_immunity 
% produces antibody titres after vaccination
% Figure 3 in the paper
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [] =  plotFigure3_population_immunity()

load('mcmc/fullmod_mcmc_output');
ps = PosteriorSamples; %posterior distributions
sys_par = sys_par;%McMC system parameters 


%setup parameters
pars = InitParameters();
pars.alpha = 0/365; %immune decay
pars.proj = 'HK_lockdown';
%setup initial condition
%[yini age_arr] = make_ics_naive( pars, pars.arrSlu, pars.arrIlu, pars.arrC Ilu);
[yini age_arr pars] = make_ics( pars, pars.arrSlu, pars.arrIlu, pars.arrCIlu);
pars.age_arr = pars.age_arr;

%set up simulation time
times_vac = 0:1:90; %3 months full vaccination and 1 month booster
%options = odeset('RelTol',1e-20,'AbsTol',1e-4); %Setup the tolerence, otherwise get the wrong results.  
[t1 y_pre] = ode23(@(t1,y1)odef_prevac_mod(t1,y1,pars), times_vac, yini);
pars.y0 = y_pre(end,:);
% Plot pre-existing immunity
%figure('Renderer', 'painters', 'Position', [20 500 780 600]); % for Mac
figure('Renderer', 'painters', 'Position', [20 200 1050 800]); % for Windows
subplot(2,2,1);
s_time = 1; %60 days later at 31 March
plot_infecteds_distribution_prevac(y_pre, pars, times_vac, times_vac(end));


%% model parameter
y_pre = pars.y0;
burnIn = 100;
totaldays = 60;
times = 0:1:totaldays;
sample_time = 24; % orginal 40d
pars.sample_time = sample_time;
times_sim = 0:0.5:times(end);

R_mean = mean(ps.R(burnIn:end));
tqr_mean = mean(ps.Tqr(burnIn:end));
q_mean = mean(ps.q(burnIn:end));
IFD50_mean = mean(ps.IFD50(burnIn:end));
ka_mean = mean(ps.ka(burnIn:end));
kb_mean = mean(ps.kb(burnIn:end));
beta_temperature_mean = mean(ps.beta_temperature(burnIn:end));
beta_relhumid_mean = mean(ps.beta_relhumid(burnIn:end));
ratr_mean = mean(ps.ratr(burnIn:end));
alpha1_mean = mean(ps.alpha1(burnIn:end));
beta_mob_mean = mean(ps.beta_mob(burnIn:end));


%% 建立prediction intervals

pars.R = R_mean;
pars.tqr = tqr_mean;
pars.q = q_mean;
pars.IFD50 = IFD50_mean;
pars.det_rate_matrix = get_detection_rate(IFD50_mean);
pars.ka = ka_mean;
pars.kb = kb_mean;
pars.beta_temperature = beta_temperature_mean;
pars.beta_relhumid = beta_relhumid_mean;
pars.ratr = ratr_mean;
pars.alpha1 = alpha1_mean;
pars.beta_mob = beta_mob_mean;  

%% Preseting
% set susceptibility
pars.arrh = make_h(pars);

[t y] = ode23(@(t,y)odef_islmod_m1_12_151(t,y,pars), times_sim, y_pre);

%[t y] = ode23(@(t,y)odef_islmod_booster_mob(t,y,pars), times_sim, y_pre(end,:));

t = t(1:2:end);
y = y(1:2:end,:);
%transform y into proportion
y = y./pars.N;
sum(y(45, pars.arrVlu)+y(45, pars.arrVclu)) % How many people are vaccinated

% Plot the disease and herd immunity dynamics
%figure;
subplot(2,2,2);
%figure('Renderer', 'painters', 'Position', [20 500 480 560]);
%s_time = 47; %60 days later at 31 March
%Plot Population immunity at a particular day
%plot_infecteds_distribution_stacked(y, pars, times, s_time);
plot_infecteds_distribution_stacked(y, pars, t, sample_time);

%figure;
subplot(2,2,3);
%plot_protection(ps, sys_par, pars);
plot_seroprevalence( age_arr, y, pars, t );

subplot(2,2,4);
plot_grid( age_arr, y, pars, t ); % Plot the figure as objective1 in meeting
