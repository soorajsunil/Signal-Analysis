close all; clear; clc;

% experiment parameters 
Nsnap      = 50; % batch size 

%  signal parameters ......................................................
theta    = [-20; 45]*pi/180; % true direction of arrival [degrees] 
Fc       = 1.5e3;   % centre frequency 
Fs       = DOA1.5e4;   % sampling frequency 
c        = 1500;    % speed of sound [meter/seconds] 
lambda   = c/Fc;    % wavelength 

% noise parameters ........................................................
amp_var  = 1;    % variance of amplitude noise 
SNR      = 1;   % signal to noise ratio [dB]

% sensor parameters .......................................................
Nsensors  = 8;                % number of hydrophone sensors 
d         = (1/2)*(lambda);   % hydrophone spacing 

% simulate zero-mean signals ..............................................
Ntargets = length(theta);             % Number of narrow band targets
amp_var  = amp_var*rand(Ntargets, 1); % set random noise values for each target or manually set 
st       = sqrt(amp_var).*randn(Ntargets, Nsnap) ... % complex signal 
                     .*exp(1j*(2*pi*Fc*repmat((1:Nsnap)/Fs, Ntargets, 1)));
% plot signals 
figure('name', 'Simulated signals')
for n = 1:Ntargets
    subplot(Ntargets, 2,  2*n-1)
    plot(real(st(n,:)), 'r-'); hold on; 
    plot(ones(1,Nsnap)*sqrt(amp_var(n,:)), 'k--')
    plot(-ones(1,Nsnap)*sqrt(amp_var(n,:)), 'k--')
    ylabel('Re s(t)')
    xlabel('t (sec)')
    title(['Signal: ', num2str(n)])

    subplot(Ntargets, 2, 2*n)
    plot(imag(st(n,:)), 'r-'); hold on; 
    plot(ones(1,Nsnap)*sqrt(amp_var(n,:)), 'k--')
    plot(-ones(1,Nsnap)*sqrt(amp_var(n,:)), 'k--')
    ylabel('Im s(t)')
    xlabel('t (sec)')
    title(['Signal: ', num2str(n)])
end 

% simulate measurements....................................................
steerMat = exp(-1j*(0:Nsensors-1)'*2*pi*(d/lambda)*sin(theta')); % steering matrix 
amp_cov  = diag(amp_var);   % covariance of amplitude noise  
z_var    = trace(amp_cov)/10^(SNR/10); % measurement noise deviation 
zt       = zeros(Nsensors,Nsnap);             % measurement vector 

for t = 1:Nsnap
    zt(:,t)    = steerMat*st(:,t) + sqrt(z_var/2)*(randn(Nsensors,1) + 1j*randn(Nsensors,1)) ;      
end 

% plot array measurements 
figure('Name','Simulated Measuremens')
subplot(2, 1, 1)
lgd = cell(1,Nsensors); 
for n = 1:Nsensors
    plot(real(zt(n,:))); hold on; 
    lgd{n} = ['sensor : ' num2str(n)]; 
end 
legend(lgd); 
ylabel('Re z(t)')
xlabel('t (sec)')
title('Array Measurements')
subplot(2, 1, 2)
for n = 1:Nsensors
    plot(imag(zt(n,:))); hold on; 
end  
legend(lgd);
ylabel('Im z(t)')
xlabel('t (sec)')
title('Array Measurements')

% MUSIC ...................................................................

% Form sensor covariance matrix 
Zhat              = (zt*zt')/Nsnap ;   

% Eigen decomposition 
[eigVector,eigValues] = eig(Zhat); 
 Nsub                 = eigVector(:, 1:(Nsensors-Ntargets)) ; % Noise sub-space 
sweep                 = (-90:1:90)*pi/180; % brute force sweep 

P                     = zeros(1,length(sweep)) ; 
 
for t = 1: length(sweep)
    SS = (steermat(Nsensors, d, lambda, (sweep(t))))'; 
    P(t)=abs(1/ (SS*(Nsub*Nsub')*SS'));
 
end

P=10*log10(P/max(P));

figure('name','MUSIC')
grid on; box on; hold on
set(gca, 'fontsize', 14)
plot(sweep*180/pi, P,'b-', 'LineWidth',2) 
xlabel('Arrival Angle (deg)');
ylabel('Magnitude (dB)');

function [A] = steermat(M, d, lambda, theta)
% calculate the steering matrix 
%
% M      ~ number of hydrophones 
% d      ~ hydrophone spacing 
% lambda ~ wavelength 
% theta  ~ direction of arrival 

A = exp(-1j*(0:M-1)'*2*pi*(d/lambda)*sin(theta')); % steering matrix 
end



