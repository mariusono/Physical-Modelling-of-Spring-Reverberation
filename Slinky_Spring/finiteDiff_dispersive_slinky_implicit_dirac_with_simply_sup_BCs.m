clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Figures_Report');
saveFigs = 0;

fs = 44100;
k = 1/fs; % time step
L = 1;

% Loss terms
% sigma_0 = 0.01; % general loss, i.e. damping
sigma_0 = 0.1; % general loss, i.e. damping
% sigma_0 = 0.0; % general loss, i.e. damping
% sigma0 = 10.0; % general loss, i.e. damping
sigma_1 = 0.00001; % frequency dependent loss
% sigma1 = 0.0; % frequency dependent loss
% sigma1 = 0.0001; % frequency dependent loss


K = 0.06; % for slinky in Julian Parker paper
% kappa = 0.4; % for slinky in Julian Parker paper
% kappa = 1; % kappa = sqrt((E*I)/(rho*A*L^4))


% theta = 0.52;
theta = 0.7;
% theta = 1;

h = sqrt((2*K*k)/sqrt(2*theta-1));
N = floor(L/h);
h = L/N;

x = linspace(0,1,N+1);


mu = K*k/h^2



% You have to calculate from u_{1} to u_{N-1}
% Clamped ! 
Dxxxx = (1/h^4).*toeplitz([6 -4 1 zeros(1,N-4)]);
Dxxxx(1,:) = (1/h^4).*[5 -4 1 zeros(1,N-4)];
Dxxxx(end,:) = (1/h^4).*fliplr([5 -4 1 zeros(1,N-4)]);
Dxx = (1/h^2).*toeplitz([-2 1 zeros(1,N-3)]);
Mx = (1/2).*toeplitz([0 1 zeros(1,N-3)],[0 1 zeros(1,N-3)]); % Mx+ 

A = sparse((1/k^2).*(-Mx*theta + Mx + 2*k*sigma_0*eye(N-1) + theta*eye(N-1))); 
B = sparse((1/k^2).*(-2*Dxx*k*sigma_1 + Dxxxx*K^2*k^2 + 2*Mx*theta - 2*Mx - 2*theta*eye(N-1)));
C = sparse((1/k^2).*(-Mx*theta + Mx + 2*k*Dxx*sigma_1 - 2*k*sigma_0*eye(N-1) + theta*eye(N-1))); 


% % excitation = sweeptone(11,4,fs,'SweepFrequencyRange',[20 fs/2]);
% % excitation = sweeptone(2,1,fs,'SweepFrequencyRange',[20 fs/2]);
% % excitation = sweeptone(3,10,fs,'SweepFrequencyRange',[20 fs/2]);
% excitation = sweeptone(1,4,fs,'SweepFrequencyRange',[20 fs/2]);
% excitation = excitation./6;


dur = 5; % s
NF = floor(dur*fs);
excitation = zeros(NF,1);
excitation(1) = 1;

t = [0:1/fs:length(excitation).*1/fs - 1/fs];


figure(444);
plot(t,excitation,'linewidth',2)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
if saveFigs
    saveas(figure(444),fullfile(pathSave,'slinky_implicit_theta_0p7_excitation.png'))
end

% [y,fs_y] = audioread('sound.wav');
% excitation = y;
% soundsc(excitation,fs_y)


NF = length(excitation);


% % initialize grid functions and output
% Domain is from u_{1} to u_{N-1}
u0 = 1; v0 = 0;
uPrev = zeros(N-1,1); 
uPrev(1) = excitation(1);
u = uPrev; 
uNext = zeros(N-1,1);
out = zeros(NF,1);
out2 = zeros(NF,1);

%%%%%% start main loop
legendStr = {};
for n=1:NF

    out(n) = u(N-1);

    uNext = A\(-B*u-C*uPrev);
    
    
    % Add boundary conditions for simply supported case somehow ?? 
    uNext(1) = excitation(n) + uNext(1);    

    uPlot = [0;uNext;0];

    if ismember(n,[1,324,1000,44100])   
        figure(1234);
        plot(x,uPlot,'linewidth',2)
        hold all
        legendStr = cat(1,legendStr,['u at time step ',num2str(n-1)]);
        grid on
        xlabel('Scaled Length (x'') [-]')
        ylabel('Amplitude [m]')    
    end
    
    uPrev = u;
    u = uNext;
    
end

limsy=get(gca,'YLim');

figure(1234);
plot([x(N-1),x(N-1)],limsy,'k-','linewidth',0.5)
legendStr = cat(1,legendStr,['Output location']);
legend(legendStr);
if saveFigs
    saveas(figure(1234),fullfile(pathSave,'slinky_implicit_theta_0p7_displacement.png'))
end


t = 0:1/fs:(floor(length(out)/fs)-1/fs);

figure(2);
plot(t,out,'linewidth',2)
grid on
hold all
xlabel('Time [s]')
ylabel('Amplitude [m]')
if saveFigs
    saveas(figure(2),fullfile(pathSave,'slinky_implicit_theta_0p7_output.png'))
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
    saveas(figure(123),fullfile(pathSave,'slinky_implicit_theta_0p7_output_spectrum.png'))
end


% % soundsc(out,fs)



% irEstimate = impzest(excitation,out);
irEstimate = out;

% % % soundsc(irEstimate,fs)


t = [0:1/fs:length(irEstimate).*1/fs - 1/fs];

figure(333);
plot(t,irEstimate,'linewidth',2)
grid on
hold all
xlabel('Time [s]')
ylabel('Amplitude [m]')
if saveFigs
    saveas(figure(333),fullfile(pathSave,'slinky_implicit_theta_0p7_imp_resp.png'))
end


% % data_slinky = load('data_slinky_theta_0p7.mat');
% % data_slinky = data_slinky.data_slinky;

M = 256;
overlap = M/2; % overlap
% overlap = round(M*3/4); % overlap
K = M*2;

% win = ones(M,1);
win = window(@hann,M);

[stft,freq,time] = spectrogram(irEstimate,win,overlap,K,fs);
time = time - time(1);

mag = 20*log10(abs(stft));

figure(222);
p1 = surf(time,freq,mag);
view([0 90])
shading interp
xlabel('Time')
ylabel('Freq')
ylim([0,20000])
% set(gca,'yscale','log')
colorbar
caxis([-80,0])
figure(222);
xlim([0,5])

if saveFigs
    saveas(figure(222),fullfile(pathSave,'slinky_implicit_theta_0p7_imp_resp_spectrogram_full.png'))
end


xlim([0,0.3])

hold all


K = 0.06;

Td = (h.*N)./(2*sqrt(2.*pi.*freq.*K));

p3 = plot3(Td,freq,ones(size(Td)).*100,'k-','linewidth',3);
plot3(Td.*3,freq,ones(size(Td)).*100,'k-','linewidth',3)
plot3(Td.*5,freq,ones(size(Td)).*100,'k-','linewidth',3)

legend(p3,{'Echoes - Model'})

if saveFigs
    saveas(figure(222),fullfile(pathSave,'slinky_implicit_theta_0p7_imp_resp_spectrogram.png'))
end

