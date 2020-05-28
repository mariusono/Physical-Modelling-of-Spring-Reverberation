clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');

L = 1;

kappa = 1;
% kappa = 0.06;
% kappa = 10;

fs = 44100;
k = 1/fs;

thetaVec = [0.52,0.57,0.7,1,2];
% thetaVec = [1];
% thetaVec = [0.7];

beta = [0:0.001:10000];

omega_real = kappa*beta.^2;

figure(111);
plot(beta,omega_real./(2*pi),'linewidth',3)
ylim([0,44100/2])

legendStr1 = {};
legendStr2 = {};

figure(1);
plot(omega_real./(2*pi),omega_real./beta,'k-','linewidth',2)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$v_{\phi}$ [m/s]','interpreter','latex')
legendStr1 = cat(1, legendStr1, ['Model Dispersion']);
xlim([0,fs/2])
% ylim([0,400])
hold all

figure(2);
hold all
plot(omega_real(1:end-1)./(2*pi),diff(omega_real)./diff(beta),'k-','linewidth',2)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$v_{g}$ [-]','interpreter','latex')
legendStr2 = cat(1, legendStr2, ['Model Dispersion']);
xlim([0,fs/2])
ylim([0,800])

for iTheta = 1:length(thetaVec)

    theta = thetaVec(iTheta);

    h = sqrt((2*kappa*k)/sqrt(2*theta-1));
    N = floor(1/h);
    h = 1/N;

% % Other solutions:    
%     omega = acos((2.0*kappa*k.^2*sin(beta.*h).^2/h.^4 + 4.0*kappa*k.^2*cos(beta.*h)./h.^4 - 4.0*kappa*k.^2./h.^4 - theta*cos(beta.*h) + theta + cos(beta.*h))./(-theta*cos(beta.*h) + theta + cos(beta.*h)))./k;
%     omega = acos((2.0*kappa*k.^2*sin(beta.*h).^2/h.^4 + 4.0*kappa*k.^2*cos(beta.*h)./h.^4 - 4.0*kappa*k.^2./h.^4 - theta*cos(beta.*h) + theta + cos(beta.*h))./(-theta*cos(beta.*h) + theta + cos(beta.*h)))./k;
%     omega = (6.28318530717959 - acos((2.0*kappa*k.^2*sin(beta.*h).^2./h.^4 + 4.0*kappa*k.^2*cos(beta.*h)./h.^4 - 4.0*kappa*k.^2./h.^4 - theta*cos(beta.*h) + theta + cos(beta.*h))./(-theta*cos(beta.*h) + theta + cos(beta.*h))))./k;
    
    omega = acos((2.0*kappa.^2*k.^2*sin(beta.*h).^2./h.^4 + 4.0*kappa.^2*k.^2*cos(beta.*h)./h.^4 - 4.0*kappa.^2*k.^2./h.^4 - theta*cos(beta.*h) + theta + cos(beta.*h))./(-theta*cos(beta.*h) + theta + cos(beta.*h)))./k;
    
    
    figure(111);
    hold all
    grid on
    plot(beta,omega./(2*pi))
%     ylim([0,44100/2])    
    
    phase_vel = omega./beta;
    
%     figure;
%     plot(omega./(2*pi),phase_vel)
%     
%     figure;
%     plot(phase_vel)
%     
%     figure;
%     plot(diff(phase_vel))

    [~,locEnd] = max(phase_vel); 
    

    figure(1);
    plot(omega(1:locEnd)./(2*pi),phase_vel(1:locEnd),'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{\phi}$ [m/s]','interpreter','latex')
    legendStr1 = cat(1, legendStr1, ['$\theta = ',num2str(theta),'$']);
    xlim([0,fs/2])
    hold all

    group_vel = diff(omega)./diff(beta);
%     figure;
%     plot(group_vel)
    
    [~,locEnd] = max(group_vel); 

    
    figure(2);
    hold all
%     plot(omega(1:end-1)./(2*pi),diff(omega)./diff(beta),'linewidth',2)
    plot(omega(1:locEnd)./(2*pi),group_vel(1:locEnd),'linewidth',2)
    grid on
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$v_{g}$ [-]','interpreter','latex')
    legendStr2 = cat(1, legendStr2, ['$\theta = ',num2str(theta),'$']);
    xlim([0,fs/2-10])
%     set(gca,'Ylim',0)
    
end



figure(1);
legend(legendStr1,'interpreter','latex','location','best');


figure(2);
legend(legendStr2,'interpreter','latex','location','best');
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);

% saveas(figure(1),fullfile(pathSave,'ideal_bar_implicit_phase_vel_variation_with_theta.png'))
% saveas(figure(2),fullfile(pathSave,'ideal_bar_implicit_group_vel_variation_with_theta.png'))


% % 
% % 
% % figure;
% % plot(1./group_vel(1:locEnd),omega(1:locEnd)./(2*pi),'k-','linewidth',3)
% % % plot(2./group_vel(1:locEnd),omega(1:locEnd)./(2*pi),'linewidth',2)
% % xlim([0,0.3])
% % hold all
% % plot(1./(2.*sqrt(kappa.*omega_real)),omega_real./(2*pi),'r-','linewidth',2)
% % for i = 1:26
% %     plot((2.*i-1)./(2.*sqrt(kappa.*omega_real)),omega_real./(2*pi),'r-','linewidth',2)
% % end
% % 
% % xlim([0,0.3])
% % ylim([0,44100/2])
% % 
% % 
% % 
% % 
% % 
% % figure;
% % plot(1./group_vel,omega(1:end-1)./(2*pi),'linewidth',2)
% % xlim([0,0.3])
% % hold all
% % % plot(1./(2.*sqrt(kappa.*omega_real)),omega_real./(2*pi),'linewidth',2)
% % plot(1./(2.*sqrt(kappa.*omega(1:end-1))),omega(1:end-1)./(2*pi),'linewidth',2)
% % xlim([0,0.3]);
% % 

data_slinky.group_vel = group_vel;
data_slinky.phase_vel = phase_vel;
data_slinky.omega = omega;
data_slinky.omega_real = omega_real;
data_slinky.beta = beta;

save('data_slinky_theta_0p7.mat','data_slinky')














