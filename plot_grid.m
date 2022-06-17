% plot_grid
% plot the antibody titres for susceptible individuals over time
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ serodist_mat ] = plot_grid( agearr, y, pars, times, agegroup )
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
%normialze 
nf = sum(serodist_tot,2);
serodist_tot = serodist_tot./nf(1);
serodist_tot = serodist;
% example data
titres = length(serodist_tot(1,:)); % number of titre levels 
x = linspace(min(x),max(x),length(x)); % return the index of x (times) on 2D grid
y = linspace(0,titres,titres+1); % return the index of y on 2D grid

[X,Y] = meshgrid(x,y); % 建立二維的Grid
%heat = X.*Y; % titre values
serodist_tot(:,end+1) = serodist_tot(:,end); %for using surf to plot, the last column should duplicate otherwise
                                             %will miss the last column                       
%if pars.maxi>2
    heat = serodist_tot(:,1:end)'*100; % ignore the naive group
%else
%    heat = serodist_tot(:,1:end-1)'; % include the naive group
%end
% the plot
%figure;
h = surf(X,Y,heat);
%h = mesh(X,Y,heat);
%h = surf(a)
set(h, 'edgecolor','none')
view(0,90);
colormap(hot);
colorbar;

    set(gca,'ylim',[0 titres]); % display the naive and other groups
    set(gca, 'YTick', 0:titres);

xlabel('Days');
ylabel('Antibody titre');
set(gca,'yticklabel',{'<1:10','1:10','1:20','1:40','1:80','1:160','\geq1:320',''},'XTickLabelRotation',90);
set(gca,'xtick',[0 14 28 42 58],'XTickLabel',{'01/02','15/02','01/03','15/03','31/03'});
zlabel('Proportion (%)');
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
