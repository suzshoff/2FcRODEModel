clear;clc;close all;

%%%% Pick these %%%%
% calculateAndSave = true;
calculateAndSave = false;

tissue = "blood";
fcr1 = '2aH'; % '2aH' or '2aR'
fcr2 = '3aV'; % '3aV' or '3aF'

wantBaseline =false;
%wantPoint = false;

%%
%%%%%%%%%%%%%%%%%%%%

% Change to user path
file_path = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/2FcRODEModel/Surfaces/Surfaces mat files';
name = strcat(tissue,'_',fcr1,'_',fcr2);
withPath = fullfile(file_path, name);

if calculateAndSave

    n = 25;
    multiplier = logspace(log10(.001), log10(1000), n);

    [p, y0, tspan, options] = MultipleFcRParameters(fcr1,fcr2,"HIV",tissue);

    surf2 = zeros(n);
    surf3 = zeros(n);

    igg1Parameter = 25; % IgG1
    igg3Parameter = 27; % IgG3

    % Create the waitbar for the parallel pool
    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    num_files = length(multiplier)*25;
    % Dummy call to nUpdateWaitbar to initialise
    nUpdateWaitbar(num_files, h);
    % Go back to simply calling nUpdateWaitbar with the data
    afterEach(D, @nUpdateWaitbar);


    parfor i = 1:length(multiplier)
        newP = p;

        newP(igg1Parameter) = p(igg1Parameter) * multiplier(i);

        for j = 1:25

            newP(igg3Parameter) = p(igg3Parameter) * multiplier(j);

            [~,y] = ode113(@MultipleFcRODEs, tspan, y0, options, newP);

            surf2(i,j) = sum(y(end,15:24));
            surf3(i,j) = sum(y(end,25:34));

            disp(i)
            disp(j)

            send(D, 1);

        end
    end

    save(withPath)
else
    load(withPath)
end

% Plot
%%
wantPoint = true;
calculateAndSave = false;

surf(p(igg3Parameter).*multiplier, p(igg1Parameter).*multiplier, surf2, FaceColor='#87BE7D', FaceAlpha=.5)
hold on
surf(p(igg3Parameter).*multiplier, p(igg1Parameter).*multiplier, surf3, FaceColor='#C19BF2', FaceAlpha=.5)

xlabel('[IgG3]','FontSize',18)
ylabel('[IgG1]','FontSize',18)

if tissue == 'blood'
    zlim([0 1.2e-9])
elseif tissue == 'penile'
    zlim([0 .6e-9])
elseif tissue == 'rectal'
    zlim([0 .4e-9])
end


if wantPoint
    xPointA = p(igg3Parameter)*multiplier(7);
    yPointA = p(igg1Parameter)*multiplier(19);
    if surf2(19,7) >= surf3(19,7)
        zPointA = surf2(19,7);
    else
        zPointA = surf3(19,7);
    end

    xPointB = p(igg3Parameter)*multiplier(19);
    yPointB = p(igg1Parameter)*multiplier(7);
    if surf2(7,19) >= surf3(7,19)
        zPointB = surf2(7,19);
    else
        zPointB = surf3(7,19);
    end

    plot3(xPointA,yPointA,zPointA,'o',MarkerSize=15,MarkerEdgeColor='k', MarkerFaceColor='k')
    plot3(xPointB,yPointB,zPointB,'o',MarkerSize=15,MarkerEdgeColor='r', MarkerFaceColor='r')
    legend(strcat(['FcR', fcr1, ' Complex Formation']), strcat(['FcR', fcr2, ' Complex Formation']), ...
        'Point A', 'Point B', Location='northwest',fontsize=18)
    
end

if wantBaseline
    xPointA = p(igg3Parameter);
    yPointA = p(igg1Parameter);
    z3PointA = surf3(13,13);
    z2PointA = surf3(13,13);
    plot3(xPointA,yPointA,z3PointA,'o',MarkerSize=5,MarkerEdgeColor='k', MarkerFaceColor='k')
    legend(strcat(['FcR', fcr1, ' Complex Formation']), strcat(['FcR', fcr2, ' Complex Formation']), ...
        'Baseline', 'Point B', Location='northwest',fontsize=18)
end

set(gca, 'XScale', 'log', 'YScale', 'log')

xlim([min(p(igg3Parameter).*multiplier), max(p(igg3Parameter).*multiplier)])
ylim([min(p(igg1Parameter).*multiplier), max(p(igg1Parameter).*multiplier)])

zlabel('FcR Complex Formation','FontSize',18)
%title([strcat('Tissue: ', {' '}, tissue), strcat('FcRs:', {' '}, fcr1, {', '}, fcr2)], FontSize=20)
title(mlreportgen.utils.capitalizeFirstChar(tissue), FontSize=20)

legend(strcat(['FcR', fcr1, ' Complex Formation']), strcat(['FcR', fcr2, ' Complex Formation']), Location='northwest',fontsize=18)
view(-45, 35);

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = p(igg3Parameter).*multiplier;
y = p(igg1Parameter).*multiplier;

%% Plot heatmaps (birds eye view)
figure()
surf(p(igg3Parameter).*multiplier, p(igg1Parameter).*multiplier, surf2)
if wantBaseline
    hold on
    plot3(xPointA,yPointA,z2PointA,'o',MarkerSize=5,MarkerEdgeColor='k', MarkerFaceColor='k')
end
xlabel('[IgG3]')
ylabel('[IgG1]')
set(gca, 'XScale', 'log', 'YScale', 'log')
title("FcRIIa")
colormap(parula)
view(2)
colorbar
xlim([min(p(igg3Parameter).*multiplier), max(p(igg3Parameter).*multiplier)])
ylim([min(p(igg1Parameter).*multiplier), max(p(igg1Parameter).*multiplier)])
set(gca,"FontSize",20)

figure()
surf(p(igg3Parameter).*multiplier, p(igg1Parameter).*multiplier, surf3)

if wantBaseline
    hold on
    plot3(xPointA,yPointA,z3PointA,'o',MarkerSize=5,MarkerEdgeColor='k', MarkerFaceColor='k')
end
xlabel('[IgG3]')
ylabel('[IgG1]')
set(gca, 'XScale', 'log', 'YScale', 'log')
title("FcRIIIa")
colormap(parula)
view(2)
colorbar
xlim([min(p(igg3Parameter).*multiplier), max(p(igg3Parameter).*multiplier)])
ylim([min(p(igg1Parameter).*multiplier), max(p(igg1Parameter).*multiplier)])
set(gca,"FontSize",20)


%% Plot cross section

% Calculate interpolation of surface
% At constant IgG3 (variable IgG1)

interpvalue = 100;

surf2interp_contstantg3 = interp2(x, y, surf2, ones(1,25)*interpvalue, y);
surf3interp_contstantg3 = interp2(x, y, surf3, ones(1,25)*interpvalue, y);

figure()
plot(y, surf2interp_contstantg3,'color','#87BE7D','linewidth',5)
hold on
plot(y, surf3interp_contstantg3,'color','#C19BF2','linewidth',5)

set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('[IgG1]')
ylabel('FcR complex formation')
set(gca,"FontSize",20)

xlim([min(y), max(y)])
title("At constant, intermediate IgG3")

% add baseline line
if wantBaseline
    hold on
    plot(yPointA*ones(1,2),ylim,'--k')
end
legend(["FcRIIa (ADCP)", "FcRIIIa (ADCC)", "Baseline"],'Location','northeast')

% At constant IgG1 (variable IgG3)
surf2interp_contstantg1 = interp2(x, y, surf2, x, ones(1,25)*interpvalue);
surf3interp_contstantg1 = interp2(x, y, surf3, x, ones(1,25)*interpvalue);

figure()
plot(x, surf2interp_contstantg1,'color','#87BE7D','linewidth',5)
hold on
plot(x, surf3interp_contstantg1,'color','#C19BF2','linewidth',5)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('[IgG3]')
ylabel('FcR complex formation')
set(gca,"FontSize",20)
xlim([min(x), max(x)])
title("At constant, intermediate IgG1")

% add baseline line
if wantBaseline
    hold on
    plot(xPointA*ones(1,2),ylim,'--k')
end

legend(["FcRIIa (ADCP)", "FcRIIIa (ADCC)", "Baseline"],'Location','northwest')



function p = nUpdateWaitbar(data, h)
persistent TOTAL COUNT H
if nargin == 2
    % initialisation mode
    H = h;
    TOTAL = data;
    COUNT = 0;
else
    % afterEach call, increment COUNT
    COUNT = 1 + COUNT;
    p = COUNT / TOTAL;
    waitbar(p, H);
end
end

