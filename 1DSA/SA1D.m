clear;clc;

%%%% Pick these %%%%
tissue = "blood"; % "penile" or "blood" or "rectal"
disease = "HIV";
fcr1 = '2aH';
fcr2 = '3aV';
multRange = [.05 .1 .4 1 2.5 10 20];
needToCalc = false;
%%%%%%%%%%%%%%%%%%%%

% Change to user path
file_path = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/2FcRODEModel/1DSA/mat Files';

file_name = strcat('Sensitivity Analysis', '_', fcr1, '_', fcr2, '_',...
    tissue, '_', disease);
withPath = fullfile(file_path, file_name);

[p, y0, tspan, options] = MultipleFcRParameters(fcr1,fcr2,disease,tissue);

if needToCalc

    tic;
    parfor m = 1:length(p)
        for n = 1:7
            newP = p;
            newP(m) = p(m)*multRange(n);
            [~,y] = ode113(@MultipleFcRODEs, tspan, y0, options, newP);
            output2(n,m) = sum(y(end,15:24));
            output3(n,m) = sum(y(end,25:34));
            disp(m)
            disp(n)
        end
    end
    toc;

    names = callParamNames(fcr1, fcr2);

    save(withPath)

else
    load(withPath)
end

figure(1)
hm2 = heatmap(flipud(output2));
hm2.XData = names;
hm2.YData = fliplr(multRange);
hm2.YLabel = 'Parameter Multiplier';
hm2.Title = strcat(['FcR ',fcr1,' Complex Formation']);

figure(2)
hm3 = heatmap(flipud(output3));
hm3.XData = names;
hm3.YData = fliplr(multRange);
hm3.YLabel = 'Parameter Multiplier';
hm3.Title = strcat(['FcR ',fcr2,' Complex Formation']);

exporttoprism2 = log10(output2');
exporttoprism3 = log10(output3');
