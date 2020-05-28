clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');

makePlots = 0;

m = 30;
K = 5e6;
R = 100;

% initial state:
u0 = 0.5;
v0 = 0;

% Analytical sol
w0 = sqrt(K/m);
A = u0;
B = v0/w0;

k_max = 2/w0;

fs = 44100;
k = 1/fs;
dur = 5;

N = floor(dur*fs);

t = [0:1:N-1].*k;

t_new = linspace(0,1,10000);

out_real = A.*cos(w0.*t_new) + B.*sin(w0.*t_new);



uNext = 0;
u = u0;
uPrev = u0 - v0*k;


out = zeros(N,1);
for n = 1:N

    uNext = (-2*K*k.^2*u + R*k*uPrev + 4*m*u - 2*m*uPrev)/(R*k + 2*m);    
    out(n) = uNext;
    
    uPrev = u;    
    u = uNext;    
end


if makePlots

    figure(12443);
    hold all
    plot(t,out,'r-','linewidth',2)
    grid on
    xlabel('Time [s]');
    ylabel('Amplitude [m]');
    ylim([-1,1])
    saveas(figure(12443),fullfile(pathSave,'mass_spring_real_vs_fds2.png'))

    % soundsc(out)

    freq = [0:fs/N:fs-fs/N].*2.*pi;

    figure(1234);
    plot(freq,20*log10(abs(fft(out))),'r-','linewidth',2);
    grid on
    hold all
    xlim([0,fs/2])
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]')




    figure(12443);
    subplot(2,1,1)
    hold all
    plot(t,out,'r-','linewidth',2)
    grid on
    xlabel('Time [s]');
    ylabel('Amplitude [m]');
    saveas(figure(12443),fullfile(pathSave,'mass_spring_with_damping.png'))

    freq = [0:fs/N:fs-fs/N].*2.*pi;

    subplot(2,1,2)
    plot(freq,20*log10(abs(fft(out))),'r-','linewidth',2);
    grid on
    hold all
    xlim([0,fs/2])
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]')
    saveas(figure(12443),fullfile(pathSave,'mass_spring_with_damping_2.png'))

end

