clear;clc;

% Change to user path
outputdir = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/2FcRODEModel/Boosting/Shared Boosting Excel Files/';

nM = readmatrix('Converted HIV Patient Data.xlsx', Sheet='nM', Range='B2:E106');

% Undiluted blood
nM = nM.*200;

boostlist = ["Baseline" "IgG1" "IgG3" "IgG1 and IgG3"];

% boost = "IgG1";
% boost = "IgG3";
% boost = "IgG1 and IgG3";
% boost = "Baseline";

tissue = "blood";
fcr1 = '2aR';
fcr2 = '3aF';

for k = 1:length(boostlist)

boost = boostlist(k);
% multOrAdd = "a";
multOrAdd = "Multiplication";
% multOrAdd = "Addition";
boostMult = 10; %200*100;
boostAdd = 1000;

[p, y0, tspan, options] = MultipleFcRParameters(fcr1,fcr2,"HIV",tissue);

% p([26 28]) = [median(nM(:,2)) median(nM(:,4))];

fcrFormation = zeros(length(nM),2);

parfor i = 1:length(nM)
    newp = p;
    newp([25 27]) = nM(i,[1 3]);

    if multOrAdd == "Multiplication"
        if boost == "IgG1"
            newp(25) = newp(25) * boostMult;
        elseif boost == "IgG3"
            newp(27) = newp(27) * boostMult;
        elseif boost == "IgG1 and IgG3"
            newp([25 27]) = newp([25 27]) * boostMult;
        elseif boost == "Baseline"
            newp([25 27]) = newp([25 27]);
        end
    elseif multOrAdd == "Addition"
        if boost == "IgG1"
            newp(25) = newp(25) + boostAdd;
        elseif boost == "IgG3"
            newp(27) = newp(27) + boostAdd;
        elseif boost == "IgG1 and IgG3"
            newp([25 27]) = newp([25 27]) + boostAdd;
        elseif boost =="Baseline"
            newp([25 27]) = newp([25 27]);
        end
    end

    [t,y] = ode113(@MultipleFcRODEs, tspan, y0, options, newp);

    yend = zeros(1,34);
    yend(1:34) = y(end,:);

    fcrFormation(i,:) = [sum(yend(15:24)) sum(yend(25:34))];
    disp(i);

    % fcrFormation(i,1) = totalFcr2Formation;
    % fcrFormation(i,2) = totalFcr3Formation;
end

%%
name = strcat(['Calculated HIV Patient Outputs Using Blood IgG2_4 ', fcr1, fcr2, '.xlsx']);

IgG1 = nM(:,1);
IgG2 = nM(:,2);
IgG3 = nM(:,3);
IgG4 = nM(:,4);

if multOrAdd == "Multiplication"
    if boost == "IgG1"
        IgG1 = nM(:,1) * boostMult;
    elseif boost == "IgG3"
        IgG3 = nM(:,3) * boostMult;
    elseif boost == "IgG1 and IgG3"
        IgG1 = nM(:,1) * boostMult;
        IgG3 = nM(:,3) * boostMult;
    end
elseif multOrAdd == "Addition"
    if boost == "IgG1"
        IgG1 = nM(:,1) + boostAdd;
    elseif boost == "IgG3"
        IgG3 = nM(:,3) + boostAdd;
    elseif boost == "IgG1 and IgG3"
        IgG1 = nM(:,1) + boostAdd;
        IgG3 = nM(:,3) + boostAdd;
    end
end

FcR2 = fcrFormation(:,1);
FcR3 = fcrFormation(:,2);

diff = FcR3 - FcR2;

compOutput = table(IgG1,IgG2,IgG3,IgG4,FcR2,FcR3,diff);

if multOrAdd == "Multiplication"
    if boost~= "Baseline"
        writetable(compOutput,strcat(outputdir,name),Sheet=strcat(boost, " Boost x", string(boostMult)))
    else
        writetable(compOutput,strcat(outputdir,name),Sheet=boost)
    end
elseif multOrAdd == "Addition"
    if boost~="Baseline"
        writetable(compOutput,strcat(outputdir,name),Sheet=strcat(boost, " Boost + ", string(boostAdd)))
    else
        writetable(compOutput,strcat(outputdir,name),Sheet=boost)
    end
elseif boost == "Baseline"
    writetable(compOutput,strcat(outputdir,name),Sheet=strcat(boost))
end

end