clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');

saveFigs = 0;

freq = 110;

fs = 44100;
k = 1/fs;
dur = 1;

NF = floor(dur*fs);
t = [0:1:NF-1].*k;


L = 1;
waveLength = 2*L;

c = freq*waveLength; %% Ask Silvin about this ?

h = c*k;
N = floor(L/h);
h = L/N;

C = c*k/h
lambda = c*k/h

x = linspace(0,1,N+1);

% create raised cosine
ctr = 0.5; wid = 0.1; % center location/width of excitation
xax = x.';
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));


figure(1);
plot(x,rc,'linewidth',2)
xlabel('Scaled Length (x'') [-]')
ylabel('Amplitude [m]')
grid on
if saveFigs
    saveas(figure(1),fullfile(pathSave,'oned_wave_raised_cos_excitation.png'))
end

% % initialize grid functions and output
u0 = 1; v0 = 0;
uPrev = rc.*u0;
u = (u0+k*v0).*rc; % same as uPrev but with some velocity ! 
uNext = zeros(N,1);
out = zeros(NF,1);



legendStr = {};
for n = 1:NF
        
    uNext(2:N-1) = lambda.^2.*(u(3:N)+u(1:N-2)) + ...
            2.*(1-lambda^2).*u(2:N-1) - ...
            uPrev(2:N-1);
    
    if ismember(n,[1,324,688,1000])
        figure(1234);
        plot(x,uNext,'linewidth',2)
        hold all
        legendStr = cat(1,legendStr,['u at time step ',num2str(n-1)]);
        grid on
        ylim([-1,1]) 
        xlim([0,1])
        xlabel('Scaled Length (x'') [-]')
        ylabel('Amplitude [m]')    
    end

    out(n) = uNext(floor(3*N/8));
    
    uPrev = u;
    u = uNext;
    
end


limsy=get(gca,'YLim');

figure(1234);
plot([x(floor(3*N/8)),x(floor(3*N/8))],limsy,'k-','linewidth',0.5)
legendStr = cat(1,legendStr,['Output location']);
legend(legendStr);
if saveFigs
    saveas(figure(1234),fullfile(pathSave,'oned_wave_displacement.png'))
end

figure(2);
subplot(2,1,1)
plot(t,out,'linewidth',2)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
ylim([-0.6,0.6])
subplot(2,1,2)
plot(t,out,'linewidth',2)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
xlim([0.2,0.22]);
ylim([-0.6,0.6]);
xlabel('Time [s]')
ylabel('Amplitude [m]')
title('Zoomed time interval')
if saveFigs
    saveas(figure(2),fullfile(pathSave,'oned_wave_output.png'))
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
    saveas(figure(123),fullfile(pathSave,'oned_wave_output_spectrum.png'))
end


% % soundsc(out,fs)

