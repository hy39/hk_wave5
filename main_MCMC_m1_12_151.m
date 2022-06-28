function [ output_args ] = main_MCMC_m1_12_151()
%main_MCMC The main script parameter estimation of the models using MCMC
% MCMC metropolis-hasting algorithm to estimate parameters
% The algorithm accepts or rejects the proposed state based on the density
% of the target distribution
% 1 March, 2022
% Hsiang-Yu Yuan

%% INITIALIZE THE MODEL AND ITS PARAMETERS
%ENVIRONMENTAL VARIABLES
%Assign global variables to local variables; Because global variables are not accessible when java is running.
clear all;close all

%Define model parameters
pars = InitParameters(); 

%Set running conditions
%nsteps = 1200000;
%nsteps = 200000;
%nsteps = 50000;
nsteps = 5000;
%nsteps = 1000;
%nsteps =  200;
%nsteps =  30;
burnIn = 100;
sampleno = nsteps - burnIn;
m = ['1_12_151']; % and contact tracing improvement


load('HK_virus');
%load('virus_hk07');
Cases = cases;
Cases_rat = cases_rat;
pars.Cases = Cases;
tic;

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
pars.t = [1:120]; % Total days of simulation
elapsed = toc;
disp(['Total number of iterations: ' num2str(nsteps)]);
disp(['Initial Parameters: ']);
%disp(pars);

%Initialize Prior density function for MCMC
if exist('nsteps')
    %[Par_stat] = initialModel_m1_12_14(nsteps);
    [Par_stat] = initialModel_m1_12_15(nsteps);
else
    %[Par_stat] = initialModel_m1_12_14();
    [Par_stat] = initialModel_m1_12_15();
end
pars.model = 1.12;
Par_stat.ode.mcmc.nSteps = nsteps;
Par_stat.ode.mcmc.burnIn = burnIn;
Par_stat.par = pars;
%Par_stat.mysir = @odeSIREmodel;
Par_stat.mysir = @odef_islmod_m1_12_151; % setup the model


[PosteriorSamples PriorPoint] = runMCMC2(pars.y0,Cases,Par_stat, pars.t,Cases_rat);
elapsed = toc;
disp(['Total number of iterations: ' num2str(nsteps)]);
disp(['Initial Parameters: ']);
%disp(pars);

%prompt = ['virus is selected. \n'];

disp(['virus is selected. \n']);
%% CREATE MCMC OUTPUT
mainoutdir = 'out/mcmc';

outfile = ['mcmc_output_m' num2str(m)];
modelparfile =  ['parameters_m' num2str(m) '.mat'];
out_dir_full = set_projectoutput( mainoutdir, ['/m' num2str(m)] );

%% Save the Parameters and Posterior into m file
sys_par.nSamples = sampleno;
sys_par.burnIn = burnIn;
sys_par.PriorPoint = PriorPoint;
sys_par.PriorMeta = Par_stat.ode.mcmc.Prior;
Par_stat.maxlikelihood = max(PosteriorSamples.LLH); %Add maximum likelihood

outfilename = [out_dir_full outfile '.mat'];
nofile = 0;
if exist(outfilename, 'file') == 2
   nofile = nofile + 1;
   newfileid = 2;
   for i=2:10
       outfilename = [out_dir_full outfile '(' num2str(i) ').mat'];
       if exist(outfilename) == 2
           newfileid = i+1;
       else
           break;
       end
   end
   outfilename = [out_dir_full outfile '(' num2str(newfileid) ').mat'];
end

par.out_dir = out_dir_full;

if nsteps < 10
    disp('not save output files');
else
    save([outfilename] ,'PosteriorSamples','par','sys_par','Par_stat','elapsed');
    disp('save output files');
end

%% Get figures output
    %PosteriorSamples = getOutput( out_dir_full, outfile, burnIn); % also remove MCMC results before burnIn
    %save([out_dir_full outfile '_final.mat'] ,'PosteriorSamples','Ab','par','sys_par','Par_stat','elapsed');
    save_posterior_figure(PosteriorSamples,pars,burnIn,sampleno, out_dir_full, m);
end





