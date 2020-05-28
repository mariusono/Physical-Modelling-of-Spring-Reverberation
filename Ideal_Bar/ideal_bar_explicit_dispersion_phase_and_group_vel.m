clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');

L = 1;

kappa = 1;

fs = 44100;
k = 1/fs;

muVec = [0.5,0.45,0.4,0.35,0.3];

legendStr1 = {};
legendStr2 = {};
for iMu = 1:length(muVec)

    mu = muVec(iMu);
    % lambda = 0.85;

    h = sqrt(2*kappa*k);


    beta = [0:0.001:465];

    omega = 2*asin(2*mu*sin(beta.*h./2).^2)./k;

%     figure;
%     plot(beta,omega)

    figure(1);
    plot(omega./(2*pi),omega./beta,'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{\phi}$ [m/s]','interpreter','latex')
    legendStr1 = cat(1, legendStr1, ['$\mu = ',num2str(mu),'$']);
    xlim([0,fs/2])
    hold all

    figure(2);
    hold all
    plot(omega(1:end-1)./(2*pi),diff(omega)./diff(beta),'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{g}$ [m/s]','interpreter','latex')
    legendStr2 = cat(1, legendStr2, ['$\mu = ',num2str(mu),'$']);
    xlim([0,fs/2])
%     set(gca,'Ylim',0)
    
end

omega_real = kappa*beta.^2;

figure(1);
plot(omega_real./(2*pi),omega_real./beta,'k-','linewidth',2)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$v_{\phi}$ [m/s]','interpreter','latex')
legendStr1 = cat(1, legendStr1, ['Model Dispersion']);
xlim([0,fs/2])
hold all

figure(2);
hold all
plot(omega_real(1:end-1)./(2*pi),diff(omega_real)./diff(beta),'k-','linewidth',2)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$v_{g}$ [-]','interpreter','latex')
legendStr2 = cat(1, legendStr2, ['Model Dispersion']);
xlim([0,fs/2])

figure(1);
legend(legendStr1,'interpreter','latex','location','best');


figure(2);
legend(legendStr2,'interpreter','latex','location','best');
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);

% % saveas(figure(1),fullfile(pathSave,'ideal_bar_phase_vel_variation_with_mu.png'))
% % saveas(figure(2),fullfile(pathSave,'ideal_bar_group_vel_variation_with_mu.png'))


