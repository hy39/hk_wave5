function [ hFig output_args ] = mcmc_posterior_hist(PosteriorSamples)
% DISPLAY MCMC posterior samples
hFig = figure;

vars = PosteriorSamples.Properties.VariableNames;
totalvars = length(vars);
posterior = table2array(PosteriorSamples);

for i=1:length(vars)-1 
    
    nBins = 61;
    sampleBins = min(posterior(:,i)):(max(posterior(:,i))-min(posterior(:,i)))./nBins:max(posterior(:,i));
    if totalvars<10
        subplot(totalvars-1,1,i);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
        xlim([min(posterior(:,i)) max(posterior(:,i))]);
    end

    if totalvars>=10 & totalvars<12
        subplot(6,2,i);
        minn = min(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        maxx = max(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        sampleBins = linspace(minn,maxx,nBins);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
    end
    if totalvars>=12 & totalvars<=14
        subplot(7,2,i);
        minn = min(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        maxx = max(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        if maxx == 0
          maxx = 1;
        end
        sampleBins = linspace(minn,maxx,nBins);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
    end
    if totalvars>=15 & totalvars<=16
        subplot(8,2,i);
        minn = min(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        maxx = max(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        sampleBins = linspace(minn,maxx,nBins);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
    end
    if totalvars>=17 & totalvars<=18
        subplot(9,2,i);
        minn = min(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        maxx = max(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        sampleBins = linspace(minn,maxx,nBins);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
    end
    if totalvars>=19 
        subplot(ceil(totalvars/2),2,i);
        minn = min(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        maxx = max(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        sampleBins = linspace(minn,maxx,nBins);
        counts = hist(table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        bar(sampleBins, counts/sum(counts), 'k');
        xlabel(char(vars(i)));
    end

        %minn = min(0.5+table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        %maxx = max(0.5+table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))));
        %sampleBins = linspace(minn,maxx,nBins);
        %counts = hist(0.5+table2array(PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i))), sampleBins);
        %bar(sampleBins, counts/sum(counts), 'k');
        %sxlabel(char(vars(i)));
end

set(hFig, 'Position', [10 10 800 120*totalvars])
end