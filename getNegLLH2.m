function [nLLH]= getNegLLH2(theta, theta_name, mysir, x0, observe, par, observe2)

%theta: pp
%theta_name: 
%x: input parameter names
%par: parameter object
%meser: sera model in java
%yini: initial data
    
    %Setup parameters
    %OutbreakStartingDay = par.OutbreakStartingDay; %default 120 days 
    %SamplingLastDay = par.SamplingLastDay;
    
    %Update parameters
    
    %Read Obserbed Yt

    % run ODE
    %p = par.p;
    y0=par.y0; % SEIR
    %[R0, Tqr, q];
    for th = 1:length(theta)
      if (strcmp(theta_name{th},'R'))
        par.R = theta(th);
        %par.R = 8;
      end
      if (strcmp(theta_name{th},'Tqr'))
        par.tqr = theta(th);
      end
      if (strcmp(theta_name{th},'q'))
        par.q = theta(th);
      end
      if (strcmp(theta_name{th},'IFD50'))
        par.IFD50 = theta(th);
        par.det_rate_matrix = get_detection_rate(par.IFD50);
      end
      if (strcmp(theta_name{th},'IFD50Boost'))
        par.IFD50Boost = theta(th);
        par.late_det_rate_matrix = get_detection_rate(par.IFD50*par.IFD50Boost);
      end
      if (strcmp(theta_name{th},'early_det_rate'))
        par.early_det_rate = theta(th);
      end
      if (strcmp(theta_name{th},'det_rate'))
        par.det_rate = theta(th);
      end
      if (strcmp(theta_name{th},'ka'))
        par.ka = theta(th);
        par.arrh = make_h(par);
      end
      if (strcmp(theta_name{th},'kb'))
        par.kb = theta(th);
        par.arrh = make_h(par);
      end
      if (strcmp(theta_name{th},'beta_temperature'))
        par.beta_temperature = theta(th);
      end
      if (strcmp(theta_name{th},'beta_relhumid'))
        par.beta_relhumid = theta(th);
      end
      if (strcmp(theta_name{th},'ratr'))
        par.ratr = theta(th);
      end   
      if (strcmp(theta_name{th},'alpha1'))
        par.alpha1 = theta(th);
      end
      if (strcmp(theta_name{th},'beta_mob'))
        par.beta_mob = theta(th);
      end   
      if (strcmp(theta_name{th},'AbB3mu'))
        par.AbB3mu = theta(th);
      end
      if (strcmp(theta_name{th},'AbB3sigma'))
        par.AbB3sigma = theta(th);
      end
      if (strcmp(theta_name{th},'AbB3cmu'))
        par.AbB3cmu = theta(th);
      end
      if (strcmp(theta_name{th},'AbB3csigma'))
        par.AbB3csigma = theta(th);
      end
    end
	      
    %[R0];
    %par.p = p;
    times = par.t;
    %[T,Y]=ode45(mysir,times,y0,[],par);
%    Sol = ode45(mysir,[0 times(end)-1] ,y0,[],par);
    Sol = ode23(mysir,[0 times(end)-1] ,y0,[],par);
    starttime = 0;
    %starttime = 6; % max 7
    %totaltime = 43; % Social distancing 3 at day 43
    totaltime = length(times);
    tt = linspace(starttime, starttime+totaltime-1, totaltime-starttime);
    x = deval(Sol, tt)';
    
    %che22 = x(:,16); % With link: exposed
    %chi22 = x(:,17); % With link: infectious
    %chr22 = x(:,18); % With link: recovered
    %ch =  sum(x(:,par.arrHlu(1,:,1,1)),2); % NT
    ch =  sum(x(:,par.arrHNlu(1,:,1,1)),2); % NT
    Xt_nt = [diff(ch)];
    chr =  sum(x(:,par.arrHRlu(1,:,1,1)),2); % RAT
    Xt_rat = [diff(chr)];

    
    %%
    % Calculate likelihood for the local cases link(+)
    % Before t = 36
    llh = 0;
    t_point = 26;
    % Calculate likelihood for the imported
    Xt_nt; %Predicted NT
    Yt_nt = observe; %Observed NT
    Yt_rat = observe2; %Observed RAT
    llh;
    r = 40;
    Xt1 = round(Xt_nt(2:t_point-1));
    %Xt2 = round(Xt_rat(2:t_point-1));
    Xt = Xt1;
    mu = Xt;
    p = r./(mu+r);
    pr = nbinpdf(Yt_nt(2:t_point-1),r,p);
    llh = sum(log(pr));
    nLLH1 = -llh;
    if isinf(nLLH1) % change dispersion parameter if likelihood is 0
        r = 10;
        Xt1 = round(Xt_nt(2:t_point-1));
        %Xt2 = round(Xt_rat(2:t_point-1));
        Xt = Xt1;
        mu = Xt;
        p = r./(mu+r);
        pr = nbinpdf(Yt_nt(2:t_point-1),r,p);
        llh = sum(log(pr));
        nLLH1 = -llh;
    end
    % After day 26
    llh = 0;
    Xt_nt; %Predicted NT
    llh;
    r2 = 40;%40
    Xt2 = round(Xt_nt(t_point:length(Yt_nt)));
    mu2 = Xt2;
    p2 = r2./(mu2+r2);
    pr2 = nbinpdf(Yt_nt(t_point:end),r2,p2);
    llh = sum(log(pr2));
    nLLH2 = -llh;
    if isinf(nLLH2)
        r2 = 10;%40
        Xt2 = round(Xt_nt(t_point:length(Yt_nt)));
        mu2 = Xt2;
        p2 = r2./(mu2+r2);
        pr2 = nbinpdf(Yt_nt(t_point:end),r2,p2);
        llh = sum(log(pr2));
        nLLH2 = -llh;
    end
    
    llh = 0;
    Xt_rat; %Predicted RAT
    Yt_rat = observe2;
    llh;
    r3 = 40;%40
    Xt3 = round(Xt_rat(t_point:length(Yt_rat)));
    %Xt3(1) = sum(Xt_rat(1:t_point));
    mu3 = Xt3;
    p3 = r3./(mu3+r3);
    pr3 = nbinpdf(Yt_rat(t_point:end),r3,p3);
    llh = sum(log(pr3));
    nLLH3 = -llh;
    if isinf(nLLH3)
        r3 = 2;%40
        Xt3 = round(Xt_rat(t_point:length(Yt_rat)));
        %Xt3(1) = sum(Xt_rat(1:t_point));
        mu3 = Xt3;
        p3 = r3./(mu3+r3);
        pr3 = nbinpdf(Yt_rat(t_point:end),r3,p3);
        llh = sum(log(pr3));
        nLLH3 = -llh;
    end
    
    
    nLLH = nLLH1+nLLH2+nLLH3;
    
end
