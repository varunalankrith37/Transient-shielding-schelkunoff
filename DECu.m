%% Double-Exponential Pulse - Frequency vs Time Domain Methods
% Using custom fourier/invfourier functions from Dr. Kochetov
clear all; close all; clc;

%% Time setup
dt = 1e-9;              % Time step (1 nanosecond)
t_max = 5e-6;           % Total time (5 microseconds)
t = (0:dt:t_max)';      % Time vector (COLUMN vector for fourier.m)
N = length(t);

%% Physical constants
mu_0 = 4*pi*1e-7;           % Permeability of free space (H/m)
eps_0 = 8.854e-12;          % Permittivity of free space (F/m)
Z_0 = sqrt(mu_0/eps_0);     % Impedance of free space (~377 Ohms)

%% Copper shield properties
sigma_copper = 5.96e7;      % Conductivity (S/m)
mu_r = 1;                   % Relative permeability (non-magnetic)
eps_r = 1;                  % Relative permittivity
thickness = 1e-3;           % Shield thickness (1 mm)

%% Define Double-Exponential Pulse (Lightning-Type)
alpha1 = 1.5e7;             % Fast decay rate (1/s)
alpha2 = 3e6;               % Slow decay rate (1/s)
V0 = 1;                     % Peak amplitude (1 V)

% Generate pulse
v_in = V0 * (exp(-alpha2*t) - exp(-alpha1*t)) .* (t>=0);

% Plot input
figure;
plot(t*1e6, v_in, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input: Double-Exponential Pulse');
xlim([0 5]);
fprintf('Input peak: %.6f V\n\n', max(v_in));

%% METHOD 1: Frequency-Domain Approach

fprintf('=== METHOD 1: Frequency-Domain ===\n');

% Forward Fourier transform of input
[f, V_freq] = fourier(t, v_in, 'pulse');
fprintf('Number of frequency points: %d\n', length(f));

% Calculate Schelkunoff transfer function H(f)
omega = 2*pi*f;
omega(1) = 2*pi*1e-6;  % Avoid division by zero at DC

% Material properties
mu = mu_r * mu_0;
eps = eps_r * eps_0;

% Schelkunoff equations
gamma = sqrt(1j*omega*mu.*(sigma_copper + 1j*omega*eps));
Z_s = sqrt(1j*omega*mu./(sigma_copper + 1j*omega*eps));
rho_01 = (Z_s - Z_0)./(Z_s + Z_0);
tau_01 = 2*Z_s./(Z_s + Z_0);
tau_10 = 2*Z_0./(Z_s + Z_0);
A = exp(-gamma*thickness);
H = (tau_01 .* tau_10 .* A) ./ (1 - rho_01.^2 .* A.^2);

fprintf('Max |H| = %.6e\n', max(abs(H)));
fprintf('Min |H| = %.6e\n', min(abs(H)));

% Apply shield in frequency domain
V_out_freq = V_freq .* H;

% Inverse Fourier transform to get output
[t_out_freq, v_out_freq] = invfourier(f, V_out_freq, 'pulse');

% Remove DC offset
v_out_freq = v_out_freq - mean(v_out_freq);

fprintf('Freq-domain output peak: %.6e V\n', max(abs(v_out_freq)));
fprintf('Output length: %d points\n\n', length(v_out_freq));

%% METHOD 2: Time-Domain Convolution

fprintf('=== METHOD 2: Time-Domain Convolution ===\n');

% Get impulse response h(t) from H(f)
[t_impulse, h_impulse] = invfourier(f, H, 'pulse');
h_impulse = h_impulse - mean(h_impulse);  % Remove DC offset

fprintf('Impulse response length: %d points\n', length(h_impulse));

% Check where impulse peaks
[h_peak, h_peak_idx] = max(abs(h_impulse));
fprintf('Impulse peak at t = %.3e s (index %d)\n', t_impulse(h_peak_idx), h_peak_idx);
fprintf('Impulse peak value: %.6e\n', h_peak);

% The key insight: invfourier returns 2N points
% We need to handle this carefully for convolution

% Approach: Resample impulse response to original time grid
h_resampled = interp1(t_impulse, h_impulse, t, 'linear', 0);

fprintf('Resampled impulse length: %d points\n', length(h_resampled));

% Perform linear convolution with zero-padding
N_fft = length(v_in) + length(h_resampled) - 1;
V_in_fft = fft(v_in, N_fft);
H_fft = fft(h_resampled, N_fft);
v_out_conv_full = real(ifft(V_in_fft .* H_fft));

% Extract 'same' portion (centered)
center_start = floor(length(h_resampled)/2);
v_out_conv = v_out_conv_full(center_start + (1:length(v_in)));

% CRITICAL: Apply dt scaling for proper convolution integral
% Convolution: y(t) = ∫ x(τ)h(t-τ)dτ
% Discrete: y[n] ≈ Σ x[k]h[n-k]Δt
v_out_conv = v_out_conv * dt;

% Remove DC offset
v_out_conv = v_out_conv - mean(v_out_conv);

fprintf('Conv-domain output peak: %.6e V\n', max(abs(v_out_conv)));
fprintf('Output length: %d points\n\n', length(v_out_conv));

%% Compare Methods

fprintf('=== COMPARISON ===\n');

% Both outputs should be on original time grid
% Interpolate freq-domain output to match conv-domain output
v_out_freq_interp = interp1(t_out_freq, v_out_freq, t, 'linear', 0);

% Calculate error
error_abs = abs(v_out_freq_interp - v_out_conv);
max_error = max(error_abs);
rms_error = sqrt(mean(error_abs.^2));

fprintf('Max absolute error: %.6e V\n', max_error);
fprintf('RMS error: %.6e V\n', rms_error);
fprintf('Relative error: %.2f%%\n\n', 100*max_error/max(abs(v_out_freq_interp)));

% Check if methods match
if max_error < 1e-12
    fprintf('✓ METHODS MATCH! Both give the same result.\n\n');
else
    fprintf('✗ Methods do NOT match - still debugging needed.\n\n');
end

%% Plot Results

% Input
figure('Position', [100 100 1400 900]);

subplot(2,3,1);
plot(t*1e6, v_in, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input: Double-Exponential Pulse');
xlim([0 5]);

% Impulse Response
subplot(2,3,2);
plot(t_impulse*1e9, h_impulse, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (ns)');
ylabel('h(t)');
title('Impulse Response (short scale)');
xlim([0 100]);

subplot(2,3,3);
plot(t_impulse*1e6, h_impulse, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('h(t)');
title('Impulse Response (full scale)');
xlim([0 10]);

% Frequency-Domain Output
subplot(2,3,4);
plot(t_out_freq*1e6, v_out_freq, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output: Frequency-Domain Method');
xlim([0 10]);

% Time-Domain Output
subplot(2,3,5);
plot(t*1e6, v_out_conv, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output: Time-Domain Convolution');
xlim([0 5]);

% Error
subplot(2,3,6);
plot(t*1e6, error_abs, 'm-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Absolute Error (V)');
title(sprintf('Error Between Methods\nMax: %.2e V', max_error));
xlim([0 5]);

sgtitle('Double-Exponential Pulse: Method Comparison', 'FontSize', 14, 'FontWeight', 'bold');

% Overlay comparison
figure;
plot(t*1e6, v_out_conv, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Time-Domain Conv');
hold on;
plot(t_out_freq*1e6, v_out_freq, 'b--', 'LineWidth', 2, 'DisplayName', 'Freq-Domain');
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output Comparison: Both Methods Overlaid');
legend('Location', 'best');
xlim([0 10]);

%% Transfer Function Plot
figure;
subplot(2,1,1);
semilogx(f/1e6, abs(H), 'b-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('|H(f)|');
title('Transfer Function Magnitude');
xlim([0.1 1000]);

subplot(2,1,2);
SE_dB = -20*log10(abs(H));
semilogx(f/1e6, SE_dB, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('Shielding Effectiveness (dB)');
title('Shielding Effectiveness');
xlim([0.1 1000]);
ylim([0 500]);

sgtitle('1mm Copper Shield - Transfer Function', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('=== ANALYSIS COMPLETE ===\n');