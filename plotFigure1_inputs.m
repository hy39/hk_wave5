% plotFigure1_inputs
% Show the daily changes in vaccination, mobility, temperature and relative
% humidity
% Figure 1C in the paper
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [] = plotFigure1_inputs()
load('HK_virus');

%figure;
figure('Renderer', 'painters', 'Position', [20 500 1080 400]); %for Windows
subplot(2,2,1);
hold on;
br3 = bar(1:length(dailybnt)-7, [dailybnt(8:end)./7.48E06.*100 dailysinovac(8:end)./7.48E06.*100],'stacked');
ylabel({'Daily booster', 'rate (%)'});
legend([br3],'BNT','CoronaVac');
xlim([0 60]);
set(gca,'xtick',[1 15 29 43 59],'XTickLabel',{});
set(gca,'FontSize',16);
box on;

subplot(2,2,2);
hold on;
scatter(1:length(temperature), temperature, 'o', 'MarkerFaceColor', 'g');
plot(temperature,'g', 'LineWidth', 1.2);
yline(mean(temperature(1:28)),':',{''}, 'LineWidth', 1.2);
yline(mean(temperature(29:end)),'--',{''}, 'LineWidth', 1);
ylabel({'Daily mean', 'temperature (°C)'});
xlim([0 60]);
set(gca,'xtick',[1 15 29 43 59],'XTickLabel',{});
set(gca,'FontSize',16);
box on;

subplot(2,2,3);
hold on;
scatter(1:length(dailymobility), -mobility_7d*100+100, 'o', 'MarkerFaceColor', 'r');
b1 = plot(-dailymobility*100+100,'b', 'LineWidth', 1.2);
b2 = plot(-mobility_7d*100+100,'r', 'LineWidth', 1.5);
hxl = xline(1,'--',{'T1'}, 'LineWidth', 1.2, 'Alpha', 0);
hxl.FontSize = 16;
hxl2 = xline(10,'--',{'T2'}, 'LineWidth', 1.2);
hxl2.FontSize = 16;
hxl3 = xline(24,'--',{'T3'}, 'LineWidth', 1.2);
hxl3.FontSize = 16;
ylabel({'Daily', 'mobility index (%)'});
xlim([0 60]);
date0 = ([1 15 29 43 59]);
date = ({'01/02','15/02','01/03','15/03','31/03'});
set(gca,'xtick',date0,'XTickLabel',date);
xtickangle(90);
legend([b1 b2],'Raw data','7 day moving averaged')
set(gca,'FontSize',16);
box on;

subplot(2,2,4);
hold on;
scatter(1:length(relhumid), relhumid, 'o', 'MarkerFaceColor', 'r');
plot(relhumid,'r', 'LineWidth', 1.2);
yline(mean(relhumid(1:28)),':',{}, 'LineWidth', 1.2);
yline(mean(relhumid(29:end)),'--',{}, 'LineWidth', 1);
ylabel({'Daily relative', 'humidity (%)'});
xlim([0 60]);

date0 = ([1 15 29 43 59]);
date = ({'01/02','15/02','01/03','15/03','31/03'});
set(gca,'xtick',date0,'XTickLabel',date);
xtickangle(90);
set(gca,'FontSize',16);
box on;
end


