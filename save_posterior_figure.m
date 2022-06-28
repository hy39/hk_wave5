function [  ] = save__posterior_figure(PosteriorSamples,par,BurnIn,sampleno, out_dir, m )
%% Save the figure and close the figure windows
disp('save posterior');
[FigH] = mcmc_posterior_hist(PosteriorSamples);
print('-dpng','-r0',[out_dir 'mcmc_posterior_histogram_m' num2str(m) '.png'])
savefig(FigH,[out_dir 'mcmc_posterior_histogram_m' num2str(m) '.fig']);
close(FigH);
if (exist('sampleno'))
    SampleNo = sampleno;
else
    SampleNo = 20;    
end
end