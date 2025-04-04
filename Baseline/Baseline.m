
%%%% Pick these %%%%
tissue = "blood";
fcr1 = '2aH';
fcr2 = '3aV';

[p, y0, tspan, options] = MultipleFcRParameters(fcr1,fcr2,"HIV",tissue);


%%
[t,y] = ode113(@MultipleFcRODEs, tspan, y0, options, p);

names = callParamNames(fcr1, fcr2);

%%
ybase = zeros(1,34);
ybase(1:34) = y(end,:);

% Add on the sum of complexes
ybase(35:36) = [sum(ybase(15:24)) sum(ybase(25:34))];

rowfcr2Sums = sum(y(:, 15:24), 2);
rowfcr3Sums = sum(y(:, 25:34), 2);

complexname=["ag-IgG1", "ag-IgG2", "ag-IgG3", "ag-IgG4",...
    "ag-IgG1-IgG1", "ag-IgG1-IgG2","ag-IgG1-IgG3", "ag-IgG1-IgG4",...
    "ag-IgG2-IgG2", "ag-IgG2-IgG3","ag-IgG2-IgG4","ag-IgG3-IgG3",...
    "ag-IgG3-IgG4","ag-IgG4-IgG4","FcR2-ag-IgG1-IgG1",...
    "FcR2-ag-IgG1-IgG2","FcR2-ag-IgG1-IgG3", "FcR2-ag-IgG1-IgG4",...
    "FcR2-ag-IgG2-IgG2", "FcR2-ag-IgG2-IgG3","FcR2-ag-IgG2-IgG4",...
    "FcR2-ag-IgG3-IgG3","FcR2-ag-IgG3-IgG4","FcR2-ag-IgG4-IgG4",...
    "FcR3-ag-IgG1-IgG1",...
    "FcR3-ag-IgG1-IgG2","FcR3-ag-IgG1-IgG3", "FcR3-ag-IgG1-IgG4",...
    "FcR3-ag-IgG2-IgG2", "FcR3-ag-IgG2-IgG3","FcR3-ag-IgG2-IgG4",...
    "FcR3-ag-IgG3-IgG3","FcR3-ag-IgG3-IgG4","FcR3-ag-IgG4-IgG4",...
    "All FcR2 complexes","All FcR3 complexes"];

   % "IgG1", "IgG2", "IgG3", "IgG4", "ag", "FcR2",...
   % "FcR3", ...

X = categorical(complexname);
X = reordercats(X, complexname);

figure()
bar(X(1:36), ybase(1:36))
set(gca, 'YScale', 'log')
title("Baseline")

% Plot timeline for each complex
plot(t,y)
xlim([0 1000])
legend(complexname(1:34))

% Plot summed FcR2 and FcR3 complexes
plot(t,[rowfcr2Sums rowfcr3Sums])
xlim([0 1000])
legend(complexname(35:36))




