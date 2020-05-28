clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');
saveFigs = 0;


fs = 41000;
k = 1/fs;

L = 1; % It is the scaled case so always keep at 1 ! 

dur = 3;
NF = floor(dur*fs);
t = [0:1:NF-1].*k;

K = 1;

% theta = 0.52;
theta = 0.7;
% theta = 1;



h = sqrt((2*K*k)/sqrt(2*theta-1));
N = floor(L/h);
h = L/N;


x = linspace(0,1,N+1);

% create raised cosine
ctr = 0.5; wid = 0.2; % center location/width of excitation
xax = x.';
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));


figure(1);
plot(x,rc,'linewidth',2)
xlabel('Scaled Length (x'') [-]')
ylabel('Amplitude [m]')
grid on
if saveFigs
    saveas(figure(1),fullfile(pathSave,'ideal_bar_raised_cos_excitation.png'))
end


Dxxxx = (1/h^4).*toeplitz([6 -4 1 zeros(1,N-3)]);
Mx = (1/2).*toeplitz([0 1 zeros(1,N-2)],[0 1 zeros(1,N-2)]);

A = sparse((1/k^2).*(Mx-Mx.*theta+theta.*eye(N))); 
B = sparse((1/k^2).*(Dxxxx.*K^2.*k^2 + 2.*Mx.*theta - 2.*Mx - 2.*theta.*eye(N)));
C = sparse((1/k^2).*((Mx-Mx.*theta+theta.*eye(N)))); 

rp = [0.4]; % positions of readout(0-1)
rp_int = floor(N*rp); % rounded grid index for readout


% % initialize grid functions and output
u0 = 1; v0 = 0;
uPrev = rc.*u0;
u = (u0+k*v0).*rc; % same as uPrev but with some velocity ! 
uNext = zeros(N,1);
out = zeros(NF,1);
out2 = zeros(NF,1);

legendStr = {};
for n=1:NF

    uNext = A\(-B*u-C*uPrev);

    uNext(1) = 0;
    uNext(2) = 0;
    uNext(end) = 0;
    uNext(end-1) = 0;

    if ismember(n,[1,324,688,1000])
    
        figure(1234);
        plot(x,uNext,'linewidth',2)
        hold all
        legendStr = cat(1,legendStr,['u at time step ',num2str(n-1)]);
        grid on
        ylim([-1,1]) 
        xlabel('Scaled Length (x'') [-]')
        ylabel('Amplitude [m]')    
    end
    
%     out(n,:) = (1-rp_frac).*uNext(rp_int)'+rp_frac.*uNext(rp_int+1)'; % readout
    out(n,:) = uNext(rp_int); % readout
    
    
    uPrev = u;
    u = uNext;
end

limsy=get(gca,'YLim');

figure(1234);
plot([x(rp_int),x(rp_int)],limsy,'k-','linewidth',0.5)
legendStr = cat(1,legendStr,['Output location']);

legend(legendStr);
if saveFigs
    saveas(figure(1234),fullfile(pathSave,'ideal_bar_implicit_theta_0p7_displacement.png'))
end

figure(2);
subplot(2,1,1)
plot(t,out,'linewidth',2)
grid on
hold all
xlabel('Time [s]')
ylabel('Amplitude [m]')
ylim([-0.6,0.6])
subplot(2,1,2)
plot(t,out,'linewidth',2)
grid on
hold all
xlabel('Time [s]')
ylabel('Amplitude [m]')
xlim([0.2,0.3]);
ylim([-0.6,0.6]);
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Zoomed time interval')
if saveFigs
    saveas(figure(2),fullfile(pathSave,'ideal_bar_implicit_theta_0p7_vs_theta_1_output.png'))
end

f = [0:fs/NF:fs-fs/NF];

psdest = 1/(length(out)*fs)*abs(fft(out)).^2;

figure(123);
plot(f,20.*log10(abs(fft(out))),'linewidth',2)
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([0,fs/2])
set(gca,'xscale','log')
grid on
hold all
if saveFigs
    saveas(figure(123),fullfile(pathSave,'ideal_bar_implicit_theta_0p7_vs_theta_1_output_spectrum.png'))
end

% % soundsc(out,fs)

audiowrite(['ideal_bar_implicit_theta_1.wav'],out,fs)

