clear;clc;

tissue = "blood"; % "blood" or "rectal" or "penile"
fcrCombos = ["2aH", "3aV";
    "2aH", "3aF";
    "2aR", "3aV";
    "2aR", "3aF"];

pointValues = zeros(4);

for k = 1:4
    fcr1 = fcrCombos(k,1);
    fcr2 = fcrCombos(k,2);

    % Change to user path
    file_path = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/2FcRODEModel/Surfaces/Surfaces mat files';
    name = strcat(tissue,'_',fcr1,'_',fcr2);
    withPath = fullfile(file_path, name);

    load(withPath)
    
    pointValues(k,1) = surf2(19,7);
    pointValues(k,2) = surf3(19,7);

    pointValues(k,3) = surf2(7,19);
    pointValues(k,4) = surf3(7,19);
end

pointValues = array2table(pointValues);
pointValues.Properties.VariableNames = ["PointA_FcR2", "PointA_FcR3", ...
    "PointB_FcR2", "PointB_FcR3"];
pointValues.Properties.RowNames = ["2aH_3aV", "2aH_3aF", "2aR_3aV", "2aR_3aF"];

ratios = zeros(4,2);

ratios(:,1) = pointValues.PointA_FcR3 ./ pointValues.PointA_FcR2;
ratios(:,2) = pointValues.PointB_FcR3 ./ pointValues.PointB_FcR2;

diff = zeros(4,2);

diff(:,1) = (pointValues.PointA_FcR3 - pointValues.PointA_FcR2) ./ ...
    (pointValues.PointA_FcR3 + pointValues.PointA_FcR2);
diff(:,2) = (pointValues.PointB_FcR3 - pointValues.PointB_FcR2) ./ ...
    (pointValues.PointB_FcR3 + pointValues.PointB_FcR2);

logRatio = log10(ratios);

excess = zeros(4,2);
excess(:,1) = (pointValues.PointA_FcR2./(pointValues.PointA_FcR3 + pointValues.PointA_FcR2))-0.5;
excess(:,2) = (pointValues.PointB_FcR2./(pointValues.PointB_FcR3 + pointValues.PointB_FcR2))-0.5;

figure(1)
b1 = bar(excess);
b1(1).FaceColor = 'k';
b1(2).FaceColor = 'r';
title(strcat(tissue, ': Excess of complex (0 indicates FcR2 = FcR3'), FontSize=20)
legend('Point A', 'Point B', fontSize=15)
set(gca, 'Xtick', 1:4, ...
    'XTickLabel', ["2aH 3aV", "2aH 3aF", "2aR 3aV", "2aR 3aF"], ...
    fontSize=15)
ylim([-0.5 0.5])
