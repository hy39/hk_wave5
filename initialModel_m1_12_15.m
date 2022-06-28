function [ Par_stat ] = initialModel_m1_12_148(nsteps)
%InitialModel3 Summary of this function goes here
%   Detailed explanation goes here
% Include a second tightening event
Par_stat = [];
PriorsPDF = {};


%% My ideal full model
    % if add a new varialble, then check and update variable names in 1) InitParameters, 2)odef_model, 
    % 3) getNegLLH, 4) plot_HK_dynamics 
    
% Initial conditions were taken from inc_qe1/mcmc/out/mcmc/m1/20200612
% [0.05-3]
    %9_10
    x0 =     [ 5         82.5     0.114      50000         0.5                  1.25        3.3561       -0.15               0.02              0.05       0.98             1];
    var =    [ 0.12      1        0.02       5000          0.018                0.45        0.255        0.0035              0.0025              0.003      0.02             0.02]; % This is standard deviation
    lb =     [ 3.0       1        0.05       2000          0.1                  0.05        2.5          -0.4                -0.05              0.005      0.8              0.9]; 
    ub =     [ 7.0       150      0.3        200000        0.5                  5.0         5.0          0.4                 0.05               0.06       1.2              1.1];  
  
    x_name = {'R'        'Tqr'    'q'        'IFD50'       'early_det_rate'     'ka'        'kb'        'beta_temperature'  'beta_relhumid'    'ratr'      'alpha1'         'beta_mob'};
    distr  = {'uniform'  'normal' 'normal'   'uniform'     'uniform'            'uniform'   'normal'    'uniform'           'uniform'          'uniform'   'normal'         'normal'};
    optb =   ub;
% control measures
% ct_b, contact tracing capacity
% social distancing relaxing
    
%% Update model metadata
    Par_stat.model = 'SEIR';
    Par_stat.ode.mcmc.x0 =     x0;
    Par_stat.ode.mcmc.var =    var;
    Par_stat.ode.mcmc.lb =     lb;
    Par_stat.ode.mcmc.ub =     ub;
    Par_stat.ode.mcmc.x_name = x_name;
    Par_stat.ode.mcmc.distr  = distr;
    
%% Define Prior and Log Likelihood 
    for i=1:length(var)
      Prior(i).init = x0(i);			% initial values
      Prior(i).name = char(x_name(i));
      Prior(i).distribution = char(distr(i));
      Prior(i).var = var(i);
      Prior(i).lb = lb(i);
      Prior(i).ub = ub(i);
      if(strcmp(Prior(i).distribution,'uniform') || strcmp(Prior(i).distribution,'Uniform'))
        Prior(i).pdf = @unifpdf;
        PriorPDF = @(x) Prior(i).pdf(x,Prior(i).lb,Prior(i).ub);
      end
      if(strcmp(Prior(i).distribution,'normal') || strcmp(Prior(i).distribution,'Normal'))
        Prior(i).norm_mu = x0(i);
        Prior(i).norm_sd = var(i);
        Prior(i).pdf = @truncate_normpdf;
        PriorPDF = @(x) Prior(i).pdf(x,Prior(i).norm_mu,Prior(i).norm_sd,Prior(i).lb,Prior(i).ub);
      end
      %PriorPDF = @(x) Prior(i).pdf(x,Prior(i).lb,Prior(i).ub);
      %%Call priorpdf in lib/model
      PriorsPDF(i) = {PriorPDF};
    end
    Par_stat.ode.mcmc.Prior = Prior;
    Par_stat.ode.mcmc.PriorsPDF = PriorsPDF;
    
    %par = setParameters(par,'tau',5.2);  %incubation period


end


 