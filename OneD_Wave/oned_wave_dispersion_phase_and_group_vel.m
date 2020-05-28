clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');

saveFigs = 0;

L = 1;
c = 330;
gamma = c/L;

fs = 44100;
k = 1/fs;

lambdaVec = [1,0.95,0.9,0.85,0.8];

legendStr1 = {};
legendStr2 = {};
for iLambda = 1:length(lambdaVec)

    lambda = lambdaVec(iLambda);

    h = gamma*k/lambda;


    beta = [0:0.001:418];

    omega = 2./k .* asin(lambda.*sin(beta.*h/2));

    figure(1);
    plot(omega./(2*pi),omega./beta./gamma,'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{\phi}/\gamma$ [-]','interpreter','latex')
    legendStr1 = cat(1, legendStr1, ['$\lambda = ',num2str(lambda),'$']);
    xlim([0,fs/2])
    hold all

    figure(2);
    hold all
    plot(omega(1:end-1)./(2*pi),diff(omega)./diff(beta)./gamma,'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{g}/\gamma$ [-]','interpreter','latex')
    legendStr2 = cat(1, legendStr2, ['$\lambda = ',num2str(lambda),'$']);
    xlim([0,fs/2])
    
end

figure(1);
legend(legendStr1,'interpreter','latex','location','best');


figure(2);
legend(legendStr2,'interpreter','latex','location','best');
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);

% % saveas(figure(1),fullfile(pathSave,'phase_vel_variation_with_lambda.png'))
% % saveas(figure(2),fullfile(pathSave,'group_vel_variation_with_lambda.png'))


