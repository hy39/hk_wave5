%plot_grid
% plot seroprevalence over time
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ serodist_mat ] = plot_seroprevalence( agearr, y, pars, times, agegroup )
arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;

%Setup some standard auxilliary functions
tmin = min(times);
tmax = max(times);
x = times;

%define age groups
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
serodist = zeros(length(times),1);
serodist_mat = []; % [age x times x titre]
for a = 1:pars.maxa
    serodist = gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = reshape(serodist, [1, size(serodist)]);  % [age x times x titre]
end

% Total
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
local_daily_S = sum(serodist_tot(:,4:7),2)*100;
plot([0:58],local_daily_S(1:59),'LineWidth',2,'Color',[0,0,0]);

xlabel('Days');
ylabel({'Seroprevalence', 'in susceptible (%)'});
set(gca,'xtick',[0 14 28 42 58],'XTickLabel',{'01/02','15/02','01/03','15/03','31/03'},'XTickLabelRotation',90);
xlim([0 58]);
set(gca,'FontSize',16);
end

function [ titremat] = gen_strain_titres( y, times, agestates, pars )
  nots = times;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      for l=1:nostates
          a = agestates;
          for m = 1:pars.maxj
              for n = 1:pars.maxk
                %titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l,m,n)) + y(i,pars.arrRlu(a,l,m,n)) + y(i,pars.arrHlu(a,l,m,n));
                titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l,m,n))./sum(y(i,pars.arrSlu(a,1:end,m,n))); 
                %titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l,m,n));
                %titremat(i,l) = titremat(i,l) + y(i,pars.arrVlu(a,l,m,n)) + y(i,pars.arrVclu(a,l,m,n));
                %for X=1:pars.maxX
                %    titremat(i,l) = titremat(i,l) + y(i,pars.arrIlu(a,l,m,n));
                %end
              end
          end
      end
  end
end
