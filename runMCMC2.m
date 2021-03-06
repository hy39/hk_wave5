%% The function of MCMC 
%run_MCMC(fitmodel, parameters, initstate, data, mcmc_parameters)
%fitmodel: a fitmodel in java component
%parameters: Values of fixed parameters for fitmodel
%theta: Parameters to be estimated for fitmodel
%initstate: initial antibody titres in the population
%data: observed antibody titres
%mcmc_parameters: prior and proposal distributions
function [ PosteriorSamples PriorPoint] = runMCMC(yini, observe, Par_stat, t, observe2)
rng('shuffle');
Priors = Par_stat.ode.mcmc.PriorsPDF; %Function handles
PriorMeta = Par_stat.ode.mcmc.Prior;

%Initialize proposal function
mu = 0;
var = [PriorMeta.var];
x_name = {PriorMeta.name};
novar = length(var);

%Intialize Prior point
%create table http://www.mathworks.co.uk/help/matlab/ref/tableproperties.html
PriorPoint = table(PriorMeta.init); 
PriorPoint.Properties.VariableNames = x_name;
PriorPoint.Properties.RowNames = {'Initial'};
PriorPoint(end+1,:)={PriorMeta.lb};
PriorPoint.Properties.RowNames(end) = {'lb'};
PriorPoint(end+1,:)={PriorMeta.ub};
PriorPoint.Properties.RowNames(end) = {'ub'};
%PriorPoint(end+1,:)={PriorMeta.discrete};
%PriorPoint.Properties.RowNames(end) = {'discrete'};
pp = table2array(PriorPoint('Initial',:));
par = Par_stat.par;
mysir = Par_stat.mysir;

%[y0 age_arr] = make_ics_fromtitres_byage( par, par.arrSlu, par.arrIlu, par.arrCIlu, Ab(1).Abl, Ab(1).age); %use round1 data as initial condition



% Get Posterior
LogPrior = 0;
%for i=1:length(Priors)
%    LogPrior = LogPrior + log(Priors{i}(pp(i)));
%end
% There are two places use getNegLLH()
LogLikelihood = -getNegLLH2(pp, x_name, mysir, yini, observe, par, observe2);
LogPosterior = LogLikelihood + LogPrior;
Posterior = exp(LogPosterior);
prevPosterior_exponent = LogPosterior;

%% INITIALIZE THE METROPOLIS-HASTINGS SAMPLER
%% DEFINE PROPOSAL DENSITY
%q = inline('exppdf(x,mu)','x','mu');
q = @(x,v) normpdf(x,v,1);

% SOME CONSTANTS
nSteps = Par_stat.ode.mcmc.nSteps;
burnIn = Par_stat.ode.mcmc.burnIn;
minn = 0.1; 
miaxx = 5;


% INTIIALZE SAMPLER
x = zeros(nSteps,length(pp));
LPr = zeros(nSteps,length(pp));% Log likelihood
%x(1) = PriorPoint.Beta;
x(1,:) = pp; 
LogLikelihood = -getNegLLH2(pp, x_name, mysir, yini, observe, par, observe2);
LPr(1,:) = LogLikelihood;
t = 1;

prevPosterior = Posterior; 

figsamples = figure;
%% RUN METROPOLIS-HASTINGS SAMPLER
alpha_list = []; % accetance ratio list

while t < nSteps
    t = t+1
    if t == 100
      t
    end
    % SAMPLE FROM PROPOSAL
    while 1
        xStar = x(t-1,:);
        delta_q = normrnd(mu,var,1,length(pp));
        xStar = xStar + delta_q;
%        if checkPriorCondition(xStar) == 0
%          continue;
%        end
        
        % Prior conditions
        flag = 1;
        for i=1:length(Priors)
            if Priors{i}(xStar(i)) ~= 0
                flag = flag*1;
            else
                flag = flag*0;
                %disp error
            end
        end
        if flag == 1
            break;
        end
        
    end
    % CORRECTION FACTOR
    c = q(x(t-1)-xStar,mu)/q(xStar-x(t-1),mu);
 
    % CALCULATE THE (CORRECTED) ACCEPTANCE RATIO
    pp = xStar;
    LogPrior = 0;
    for i=1:length(Priors)
        LogPrior = LogPrior + log(Priors{i}(pp(i)));
    end
    %x_name 
    %pp
    %LogLikelihood = -getNegLLH(pp, x_name, mysir, yini, observe, par);
    LogLikelihood = -getNegLLH2(pp, x_name, mysir, yini, observe, par, observe2);
    LLH = LogLikelihood;
    
    
    LogPosterior = LLH + LogPrior;
    Posterior_exponent = LogPosterior;
    Posterior = exp(LogPosterior);
    alpha = min([1, (exp(Posterior_exponent-prevPosterior_exponent))*c]);
    %alpha = min([1, (Posterior/prevPosterior)*c]);
    
    %system output
    alpha_list(end+1) = alpha;  
    if rem(t,100)==0
        meanalpha = mean(alpha_list);
        disp(['steps:' num2str(t) ', LLH:' num2str(LLH) ', acceptance ratio:' num2str(meanalpha) ]);
        alpha_list = [];
    end
    % ACCEPT OR REJECT?
    u = rand;
    if u < alpha  %accept
        x(t,:) = xStar;
        prevPosterior = Posterior;
        prevPosterior_exponent = Posterior_exponent;
        LPr(t,1) = LLH;
    else          %reject
        x(t,:) = x(t-1,:);
        LPr(t,1) = LPr(t-1,1);
    end
    %LPr(t,1) = LLH; % does it cause wrong result
    
%% DISPLAY MARKOV CHAIN
     figure(figsamples);
     %if rem(t,500) == 0
     if rem(t,10) == 0
     for v=1:novar
        subplot(novar+1,1,v);
        stairs(1:2:t, x(1:2:t,v), 'k');
        ylabel(x_name{v}); 
     end
     subplot(novar+1,1,v+1);
     stairs(1:2:t, LPr(1:2:t), 'k');
     ylabel('LLH'); 
     drawnow 
     title(num2str(par.model));
     end
end

%% CREATE SIMULATION OUTPUT
% OUTPUT SIMULATION
%PosteriorSamples = table();
noburnIn = 1;
for i=1:length(pp)
    t1(:,i) = x(noburnIn:end,i);
    %t2(i)=t1.Var1;
    %varname = PosteriorSamples.Properties.VariableNames(i);
    %PosteriorSamples.Properties.VariableNames{varname} = ColNames(i);    
end
t1(:,end+1) = LPr(noburnIn:end,1); %attach likelihood
x_name = [x_name 'LLH'];
PosteriorSamples = array2table(t1, 'VariableNames', x_name);

close(figsamples); %close the figure
end
