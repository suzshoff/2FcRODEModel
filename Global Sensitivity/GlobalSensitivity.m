clear;clc;

rng(1)

%% Generate Patients
fcr1 = '2aH';
fcr2 = '3aV';
disease = "HIV";
bodycell = "blood";

[p, y0, tspan, options] = MultipleFcRParameters(fcr1,fcr2,disease,bodycell);

n = 100;

fakePatients = patientGenerator(n, .1, 5, fcr1, fcr2, disease, bodycell);

baselineThresholdMultiple1 = 1;
baselineThresholdMultiple2 = 1;

%% Solve for Outputs
fcr1Sum = zeros(n,1);
fcr2Sum = zeros(n,1);

tic;
parfor i = 1:n
% for i = 1:n
    [~,y] = ode113(@MultipleFcRODEs, tspan, y0, options, fakePatients(i,:));
    fcr1Sum(i) = sum(y(end,15:24));
    fcr2Sum(i) = sum(y(end,25:34));
    disp(i)
end
toc;

%% Calculate PRCC
clc
rhos = zeros(31,1);
pvals = zeros(31,1);
for i = 1:31
x = fakePatients(:,i);
y = fcr1Sum;
z = fakePatients;
z(:,i) = [];

[rhos(i),pvals(i)] = partialcorr(x,y,z,'Type','Pearson');
end

params = ["k1af" "k1ar" "k2af" "k2ar" "k3af" "k3ar" "k4af" "k4ar"...
    "k1r1f" "k1r1r" "k2r1f" "k2r1r" "k3r1f" "k3r1r" "k4r1f" "k4r1r"...
    "k1r2f" "k1r2r" "k2r2f" "k2r2r" "k3r2f" "k3r2r" "k4r2f" "k4r2r"...
    "igg1tot" "igg2tot" "igg3tot" "igg4tot" "agtot" "fcr1atot" "fcr2atot"];

sigIndex = pvals <= .05;
sigRhos = rhos(sigIndex);
sigNames = params(sigIndex);

bar(sigRhos)
xticks(1:length(sigIndex))
xticklabels(sigNames)

%% Find baselines
[~,y] = ode113(@MultipleFcRODEs, tspan, y0, options, p);

baseline1 = sum(y(end,15:24));
baseline2 = sum(y(end,25:34));

%% Split into above and below baseline

aboveIndex = fcr1Sum > baseline1*baselineThresholdMultiple1 & fcr2Sum > baseline2*baselineThresholdMultiple2;
belowIndex = ~aboveIndex;

abovePatients = fakePatients(aboveIndex,:);
belowPatients = fakePatients(belowIndex,:);

%% p and q values

pVals = zeros(31,1);

for i = 1:31
    pVals(i) = ranksum(abovePatients(:,i), belowPatients(:,i));
end

qVals = mafdr(pVals, 'BHFDR', true);

qPlot = -log(qVals);

%% rank difference

rankDiffs = zeros(31,1);

for i = 1:31
    firstHolder = [abovePatients(:,i) zeros(length(abovePatients(:,1)),1)];
    secondHolder = [belowPatients(:,i) ones(length(belowPatients(:,1)),1)];

    paramHolder = [firstHolder; secondHolder];

    [sorted,sortedIndex] = sortrows(paramHolder);

    aboveMean = mean(find(sorted(:,2) == 0));
    belowMean = mean(find(sorted(:,2) == 1));

    rankDiffs(i) = aboveMean - belowMean;
end

%% Plot

sgtitle(['Important Parameters for having', ' ', ...
    fcr1, ' ', 'output at', ' ', num2str(baselineThresholdMultiple1), ' ', 'times above baseline and', ' ',...
    fcr2, ' ', 'output at', ' ', num2str(baselineThresholdMultiple2), ' ', 'times above baseline'], fontsize=20)
subplot(1,3,1)

histogram(fcr1Sum)
xline(baseline1, LineWidth=3)
xline(baseline1*baselineThresholdMultiple1, LineWidth=3, Color='r')
legend('Model Output Distribution', 'Baseline', 'Threshold for Population Split')
title([fcr1, ' ', 'Complex Formation Distribution'], FontSize=15)
xlabel('Output', FontSize=15)
ylabel('Count', FontSize=15)

subplot(1,3,2)

histogram(fcr2Sum)
xline(baseline2, LineWidth=3)
xline(baseline2*baselineThresholdMultiple2, LineWidth=3, Color='r')
legend('Model Output Distribution', 'Baseline', 'Threshold for Population Split')
title([fcr2, ' ', 'Complex Formation Distribution'], FontSize=15)
xlabel('Output', FontSize=15)
ylabel('Count', FontSize=15)

params = ["k1af" "k1ar" "k2af" "k2ar" "k3af" "k3ar" "k4af" "k4ar"...
    "k1r1f" "k1r1r" "k2r1f" "k2r1r" "k3r1f" "k3r1r" "k4r1f" "k4r1r"...
    "k1r2f" "k1r2r" "k2r2f" "k2r2r" "k3r2f" "k3r2r" "k4r2f" "k4r2r"...
    "igg1tot" "igg2tot" "igg3tot" "igg4tot" "agtot" "fcr1atot" "fcr2atot"];


qthreshold = -log10(.0001);
subplot(1,3,3)
hold on
for i = 1:31
    if abs(rankDiffs(i)) > 50 && qPlot(i) > qthreshold
        scatter(rankDiffs(i), qPlot(i), 150, "red", "filled")
        text(rankDiffs(i)+.1, qPlot(i)+.1, params(i), FontSize=15)
    else
        scatter(rankDiffs(i), qPlot(i), 150, "black")
    end
end
xlabel('Average Rank Difference', FontSize=15)
ylabel('-log10(q)', FontSize=20)
title([fcr1, '/', fcr2, ' ', 'above x', num2str(baselineThresholdMultiple1), '/', ...
    num2str(baselineThresholdMultiple2), ' ', ...
    'baseline with n =', ' ', num2str(n), ' ',...
    'q-threshold of', ' ', num2str(10^(-qthreshold))], FontSize=15)
xline(0);
yline(qthreshold);
hold off
