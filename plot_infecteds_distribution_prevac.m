% plot_infecteds__distribution_prevac
% plot pre-existing immunity after vaccination
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ serodist_mat ] = plot_infecteds_distribution_prevac( y, pars, times, sampletime, agegroup )
arrSlu = pars.arrSlu;
arrIlu = pars.arrIlu;

%Setup some standard auxilliary functions
tmin = min(times);
tmax = max(times);

%define age groups
arrIlu = pars.arrIlu; %I lookup
arrSlu = pars.arrSlu; %S lookup
serodist = zeros(length(times),1);
serodist_mat = []; % [age x times x titre]
for a = 1:pars.maxa
    serodist = gen_strain_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = serodist;
end

if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end

% Total vaccinated
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
if ndims(serodist_tot3d)>2
    serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
else
    serodist_tot = serodist_tot3d;
end
%serodist_pro_tmp = serodist_tot(sampletime,:)./sum(serodist_tot(sampletime,:));
serodist_pro_tmp = serodist_tot(sampletime,:);
Yv = serodist_pro_tmp';


%%
%%subplot(1,2,2);
%%--vaccinated titres for coronavac
for a = 1:pars.maxa
    serodist = get_vacc_titres(y, times, a, pars); % [times x titre]
    serodist_mat(a,:,:) = reshape(serodist, [1, size(serodist)]);  % [age x times x titre]
end

% Total vaccinated
if ~exist('agegroup')
    agegroup = 1:pars.maxa;
end
serodist_tot3d = sum(serodist_mat(agegroup,:,:),1);  % [1 x times x titre]
if ndims(serodist_tot3d)>2
    serodist_tot = reshape(serodist_tot3d, size(serodist_tot3d,2), []); % [times x titre]
else
    serodist_tot = serodist_tot3d;
end
%serodist_pro_tmp = serodist_tot(sampletime,:)./sum(serodist_tot(sampletime,:));
serodist_pro_tmp = serodist_tot(sampletime,:);
Yvc = serodist_pro_tmp';
%Yall = Yv;
Yv = Yv./pars.N;
Yv(1) = Yv(1)*0.1;
Yv(2:end) = Yv(2:end)*4;
%Yv2 = [Yv;0]
Yv2 = [Yv]
Yall = ([Yv2 zeros(7,1)]*100);
Ytot = [zeros(7,1);sum(Yv2(4:7))];
H = bar(1:pars.maxi,Yall,'stacked','FaceColor',[.2 .2 .2],'EdgeColor','none'); %ignore titres <1:5
maxBar = max(cellfun(@max, get(H, 'YData'))); 

xlabel('Antibody titre');
set(gca,'xticklabel',{'<1:10','1:10','1:20','1:40','1:80','1:160','\geq1:320'},'XTickLabelRotation',45);
set(gca,'ytick',[0:2:10]);
set(gca,'yticklabel',{'0','0.5','1','1.5','80','100'});
ylabel('Proportion (%)');
set(gca,'FontSize',16);
end

function [ titremat] = gen_strain_titres( y, times, agestates, pars )
  nots = times;
  a = agestates;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      for l=1:nostates
          titremat(i,l) = titremat(i,l) + y(i,pars.arrSlu(a,l));
      end
  end
end

function [ titremat] = get_vac_titres( y, times, agestates, pars )
  nots = times;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      for l=1:nostates
          a = agestates;
          for m = 1:pars.maxj
              for n = 1:pars.maxk
                titremat(i,l) = titremat(i,l) + y(i,pars.arrVlu(a,l,m,n));
                %for X=1:pars.maxX
                %    titremat(i,l) = titremat(i,l) + y(i,pars.arrIlu(X,a,l,m,n));
                %end
              end
          end
      end
  end
end

function [ titremat] = get_vacc_titres( y, times, agestates, pars )
  nots = times;
  nostates = pars.maxi;
  titremat = zeros(length(nots), nostates);
  %X = 1; %only plot for the first strain
  for i=1:length(nots)
      for l=1:nostates
          a = agestates;
          for m = 1:pars.maxj
              for n = 1:pars.maxk
                titremat(i,l) = titremat(i,l) + y(i,pars.arrVclu(a,l,m,n));
                %for X=1:pars.maxX
                %    titremat(i,l) = titremat(i,l) + y(i,pars.arrIlu(X,a,l,m,n));
                %end
              end
          end
      end
  end
end


