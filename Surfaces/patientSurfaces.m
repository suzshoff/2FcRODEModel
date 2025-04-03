clear;clc;clf;

%%%% Pick these %%%%
tissue = "blood";
fcr1 = '2aH';
fcr2 = '3aV';
%boost = "IgG1 and IgG3";
multOrAdd = "Multiplication";
boostMult = 10;
boostAdd = 1000;
boostlist = ["Baseline" "IgG1" "IgG3" "IgG1 and IgG3"];

%%%%%%%%%%%%%%%%%%%%

for k=1:length(boostlist)

    boost = boostlist(k);

    if boost == "Baseline"
        boostType = "Baseline";
    else
        if multOrAdd == "Multiplication"
            boostType = strcat(boost, " Boost x", string(boostMult));
        elseif multOrAdd == "Addition"
            boostType = strcat(boost, " Boost + ", string(boostAdd));
        end
    end

    vaccineedatapath = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/KadesCode/Boosting/Shared Boosting Excel Files/';

    nM = readtable(strcat([vaccineedatapath 'Calculated HIV Patient Outputs Using Blood IgG2_4 ', fcr1, fcr2, '.xlsx']), ...
        sheet = boostType);

    igg1Points = nM.IgG1;
    igg3Points = nM.IgG3;
    fcr2Points = nM.FcR2;
    fcr3Points = nM.FcR3;

    file_path = '/Users/suzieshoffner-beck/Library/CloudStorage/Dropbox/Arnold_Lab/KadesCode/Surfaces/Surfaces mat files';
    %name = strcat(tissue,'_',fcr1,'_',fcr2,'_ExtendedRange');
    name = strcat(tissue,'_',fcr1,'_',fcr2);
    withPath = fullfile(file_path, name);

    load(withPath)
    %
    % new names
    xParameter = igg1Parameter;
    yParameter = igg3Parameter;

    figure (k)

    hold on
    surf(p(yParameter).*multiplier, p(xParameter).*multiplier, surf2, FaceColor='#87BE7D', FaceAlpha=.5)
    surf(p(yParameter).*multiplier, p(xParameter).*multiplier, surf3, FaceColor='#C19BF2', FaceAlpha=.5)
    xlabel('[IgG3]','FontSize',18)
    ylabel('[IgG1]','FontSize',18)
    zlim([0 1.2e-9])
    % xlim([0 100])
    % ylim([5 1000])
    set(gca, 'XScale', 'log', 'YScale', 'log')

    zlabel('FcR Complex Formation','FontSize',18)
%     title(['Calculated Outputs, Using Blood Parameters for IgG2/4', ...
%         strcat('Tissue: ', {' '}, tissue), ...
%         strcat('FcRs:', {' '}, fcr1, {', '}, fcr2), ...
%         boostType], FontSize=20)
    title(boostType, FontSize=20)
    view(-45, 35);

    plot3(igg3Points,igg1Points,fcr2Points,'o',MarkerSize=7,MarkerEdgeColor='k', MarkerFaceColor='#87BE7D')
    plot3(igg3Points,igg1Points,fcr3Points,'o',MarkerSize=7,MarkerEdgeColor='k', MarkerFaceColor='#C19BF2')
    for i = 1:length(igg3Points)
        line([igg3Points(i) igg3Points(i)], [igg1Points(i) igg1Points(i)], [fcr2Points(i) fcr3Points(i)],...
            LineWidth=2, Color='k')
    end
    legend(strcat(['FcR', fcr1, ' Complex Formation']), strcat(['FcR', fcr2, ' Complex Formation']), Location='northwest',FontSize=18)
    hold off

end