clc;
clear all;
close all;

[sound,fs_sound] = audioread('repeated_E.wav');
sound = sound(:,1);

inpSeq = sound;
inpSeq = cat(1,inpSeq,inpSeq);

x = inpSeq;

ana = load('C:\Users\msoi\UNIVERSITY\SEMESTER_2\Main_Project\Spring_implicit\irEstimate_K_0p16_q_600_gamma_1800.mat');
ir = ana.irEstimate;
% ir = ir(2000:2400);
% ir = ir.*linspace(1,0,length(ir)).';
ir = ir./max(abs(ir));

figure;
plot(ir)

% fftsize = 2048;
K = 512;
L = 2^ceil(log(2*K)/log(2));



N = length(ir);

noBins = ceil(N/K);
P = noBins;

ir = cat(1,ir,zeros(noBins*K-N,1));
N = length(ir);

figure;
plot(ir)

% irfft = fft(ir_dirac);

position = 0;
count = 1;
irBins = zeros(noBins,L);
while position + K <= length(ir)
    irBins(count,:) = [ir(position + (1:K));zeros(L-K,1)];
    position = position + K;
    count = count + 1;
end

S = zeros(size(irBins));
for iBin = 1:size(irBins,1)
    S(iBin,:) = fft(irBins(iBin,:));
end

S = S.'; % transpose it for consistency with the input blocks

xold = x;
x = [zeros(K,1);x]; % if you don't do this you get some lag between solutions.

position = 0;
count = 1;
step_size = L-K;
filterBlocks = [];
outBlock_real = {};
while position + L <= length(x)
    
    filterBlocks = cat(2,fft(x(position+(1:L))),filterBlocks); % stack the new bank on top of the old ones. for consistency with S.
%     startIdx = max(size(filterBlocks,2)-P,1);    
    endIdx = min(size(filterBlocks,2),P);    
%     multiplication = S(:,1:endIdx).*filterBlocks(:,startIdx:size(filterBlocks,2));
%     multiplication = S(:,1:endIdx).*filterBlocks(:,size(filterBlocks,2)-endIdx+1:size(filterBlocks,2));
    multiplication = S(:,1:endIdx).*filterBlocks(:,1:endIdx);
    outBlock = ifft(sum(multiplication,2));
    
    outBlock_real{count} = outBlock(L-K+1:end);
    
    position = position + step_size;
    count = count + 1;
end


out = cat(1,outBlock_real{:});

size(out)

soundsc(out,fs_sound)



out_real = conv(xold,ir);

figure;
hold all
plot(out_real)
plot(out)
grid on
