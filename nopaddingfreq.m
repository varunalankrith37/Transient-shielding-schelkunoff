clear all; close all; clc;

%% Time setup
dt = 1e-8;              % Time step (1 nanosecond)
t_max = 100e-6;           % Total time (5 microseconds)
t = 0:dt:t_max;         % Time vector

%% Gaussian pulse parameters
V0 = 1;                 % Peak amplitude (1 Volt)
t0 = 2e-6;              % Center time (2 microseconds)
sigma = 0.3e-6;         % Standard deviation (0.3 microseconds)

%% Generate Gaussian pulse
v_in = V0 * exp(-((t - t0).^2) / (2 * sigma^2));

%% Plot input signal
figure;
plot(t*1e6, v_in, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input Gaussian Pulse');

%% Convert to frequency domain using FFT
%{
N = length(t);              % Number of time points
V_freq = fft(v_in);         % FFT of input signal

% Create frequency vector
df = 1/t_max;               % Frequency resolution
f = (0:N-1) * df;           % Frequency vector 
%}
[f,V_freq]=fourier(t,v_in);

%% Plot frequency spectrum
figure;
plot(f/1e6, abs(V_freq), 'r-', 'LineWidth', 2);
%plot(f/1e6, real(V_freq), 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Frequency Spectrum of Gaussian Pulse');
xlim([0 50]);  

%% Physical constants
mu_0 = 4*pi*1e-7;           % Permeability of free space 
eps_0 = 8.854e-12;          % Permittivity of free space 
Z_0 = sqrt(mu_0/eps_0);     % Impedance of free space 

%% Copper shield properties
sigma_copper = 5.96e7;             % Conductivity 
mu_r = 1;                   % Relative permeability
eps_r = 1;                  % Relative permittivity
thickness = 0.5*1e-3;           % Shield thickness 
%%  Shield transfer function H(f)
% Create frequency vector with negative frequencies properly identified
f_sym = f;                                                 % Start with standard frequency vector [0, df, 2df, ..., (N-1)df]
%f_sym(f_sym > 1/(2*dt)) = f_sym(f_sym > 1/(2*dt)) - 1/dt;  % Now f_sym = [0, df, ..., fs/2, -fs/2, ..., -2df, -df] (properly unwrapped)
%f_sym(f_sym > 1/(2*dt)) = f_sym(f_sym > 1/(2*dt)) - 1/dt -df;  % Now f_sym = [0, df, ..., fs/2, -fs/2, ..., -2df, -df] (properly unwrapped)

% Angular frequency using abs(f_sym) for symmetric transfer function
omega = 2*pi*f_sym;
%omega = 2*pi*abs(f_sym);               % Use absolute value so +100MHz and -100MHz get same |H|
                                       % Without abs(): we'd only apply correct physics to positive freqs, output amplitude cut in half!
%omega(omega < 2*pi*1e-6) = 2*pi*1e-6;  % Avoid division by zero at DC and near-zero freqs
omega(abs(omega) < 2*pi*1e-6) = 2*pi*1e-6;  % Avoid division by zero at DC and near-zero freqs

% Material properties
mu = mu_r * mu_0;              % Absolute permeability (H/m)
eps = eps_r * eps_0;           % Absolute permittivity (F/m)

% Complex propagation constant (accounts for both conduction and displacement currents)
gamma = sqrt(1j*omega*mu.*(sigma_copper + 1j*omega*eps));

% Intrinsic impedance of shield material
Z_s = sqrt(1j*omega*mu./(sigma_copper + 1j*omega*eps));

% Schelkunoff's shielding theory: accounts for multiple reflections
% Reflection coefficient at air-shield interface
rho_01 = (Z_s - Z_0)./(Z_s + Z_0);

% Transmission coefficients
tau_01 = 2*Z_s./(Z_s + Z_0);    % Air to shield
tau_10 = 2*Z_0./(Z_s + Z_0);    % Shield to air

% Attenuation through shield thickness
A = exp(-gamma*thickness);

% Complete transfer function including reflections and transmission
% Denominator accounts for infinite series of internal reflections
H = (tau_01 .* tau_10 .* A) ./ (1 - rho_01.^2 .* A.^2);

% Diagnostic output
fprintf('=== Physical Parameters ===\n');
fprintf('Z_0 = %.2f Ohms\n', Z_0);
fprintf('Max |H| = %.6e (at low frequencies)\n', max(abs(H)));
fprintf('Min |H| = %.6e (at high frequencies)\n', min(abs(H)));
fprintf('Z_s at 1 MHz = %.6e Ohms\n\n', abs(Z_s(100)));

%%  Plot transfer function 
figure;
subplot(2,1,1);
%semilogx(f/1e6, abs(H), 'b-', 'LineWidth', 2);
%loglog(f/1e6, abs(H), 'b-', 'LineWidth', 2);
%semilogy(f/1e6, abs(H), 'b-', 'LineWidth', 2);
semilogy(f_sym/1e6,abs(H), 'b-', 'LineWidth', 2);
%plot(f_sym/1e6,real(H), 'b-', 'LineWidth', 2);
%plot(f_sym/1e6,imag(H), 'b-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('|H(f)|');
title('Transfer Function Magnitude');
%xlim([0.1 500]);

subplot(2,1,2);
SE_dB = -20*log10(abs(H));
%semilogx(f/1e6, SE_dB, 'r-', 'LineWidth', 2);
plot(f_sym/1e6, SE_dB, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('Shielding Effectiveness (dB)');
title('Shielding Effectiveness');
%xlim([0.1 500]);
%ylim([0 500]);

%% Multiplication in frequency domain
V_out_freq = V_freq .* H;   % Multiply input spectrum by transfer function

%% Convert back to time domain
%v_out = ifft(V_out_freq);     % IFFT to get time-domain output
%v_out = v_out - mean(v_out);        % Remove DC offset
%v_out = ifft(V_out_freq);          % Without real()

[t2,v_out]=invfourier(f,V_out_freq);

%fprintf('Max imaginary part: %.6e\n', max(abs(imag(v_out))));

% Check the output statistics
fprintf('Output mean (DC offset): %.6e V\n', mean(v_out));
fprintf('Output max: %.6e V\n', max(v_out));
fprintf('Output min: %.6e V\n', min(v_out));
fprintf('Output peak-to-peak: %.6e V\n', max(v_out) - min(v_out));

%% Plot input vs output
figure;
plot(t*1e6, v_in, 'b-', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
plot(t2*1e6, v_out, 'r-', 'LineWidth', 2, 'DisplayName', 'Output');
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Gaussian Pulse: Input vs Output (after 1mm Copper Shield)');
legend('Location', 'best');
%xlim([0 5]);

figure;
plot(t2*1e6, v_out, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output Signal (DC removed, Zoomed In)');
%xlim([0 5]);

%%  Final comparison
figure;
subplot(2,1,1);
plot(t*1e6, v_in, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input: Gaussian Pulse');
%xlim([0 5]);

subplot(2,1,2);
plot(t2*1e6, v_out, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output: After 1mm Copper Shield');
%xlim([0 5]);

sgtitle('Schelkunoff Shielding: Output vs Intput', 'FontSize', 14, 'FontWeight', 'bold');
