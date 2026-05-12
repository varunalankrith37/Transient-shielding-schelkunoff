clear all; close all; clc;

%% Shield thickness
d = 0.10;   % [m] 10 cm

%% Time setup
dt = 1e-8;
t_max = 20e-6;
t = 0:dt:t_max;

%% Gaussian pulse
V0 = 1;
t0 = 2e-6;
sigma_t = 0.3e-6;

v_in = V0 * exp(-((t - t0).^2) / (2*sigma_t^2));

%% Fourier transform
[f, V_freq] = fourier(t, v_in);

%% Constants
mu_0 = 4*pi*1e-7;
eps_0 = 8.854187817e-12;
Z_0 = sqrt(mu_0/eps_0);

%% Poor conductor / concrete properties
sigma_mat = 1e-4;     % try 1e-5, 1e-4, 1e-3
eps_r = 6;
mu_r = 1;

mu = mu_r * mu_0;
eps = eps_r * eps_0;

%% Transfer function
omega = 2*pi*f;

gamma = zeros(size(omega));
Z_s = zeros(size(omega));

idx = 2:length(omega);

gamma(idx) = sqrt(1j*omega(idx)*mu .* ...
                  (sigma_mat + 1j*omega(idx)*eps));

% Force positive attenuation
gamma(idx) = abs(real(gamma(idx))) + 1j*imag(gamma(idx));

Z_s(idx) = sqrt(1j*omega(idx)*mu ./ ...
                (sigma_mat + 1j*omega(idx)*eps));

gamma(1) = 0;
Z_s(1) = Z_s(2);

rho_01 = (Z_s - Z_0) ./ (Z_s + Z_0);
rho_10 = (Z_0 - Z_s) ./ (Z_0 + Z_s);

tau_01 = 2*Z_s ./ (Z_s + Z_0);
tau_10 = 2*Z_0 ./ (Z_s + Z_0);

A = exp(-gamma*d);

H = (tau_01 .* tau_10 .* A) ./ ...
    (1 - rho_01 .* rho_10 .* A.^2);

% Poor conductor DC handling
H(1) = H(2);

%% Impulse response
[t_h, h_full] = invfourier(f, H);

h = real(h_full(1:length(t)));
t_h = t_h(1:length(t));

% Normalized impulse response for visualization
h_norm = abs(h) / max(abs(h));

%% Time-domain convolution
v_out_conv_full = conv(v_in, h) * dt;
t_conv = (0:length(v_out_conv_full)-1) * dt;

v_out_conv = v_out_conv_full(1:length(t));
t_out = t_conv(1:length(t));

%% Frequency-domain validation
V_out_freq = V_freq .* H;
[t_fd, v_out_fd] = invfourier(f, V_out_freq);

v_out_fd = real(v_out_fd(1:length(t)));
t_fd = t_fd(1:length(t));

%% Shielding effectiveness
SE_peak_td = 20*log10(max(abs(v_in)) / max(abs(v_out_conv)));
SE_energy_td = 10*log10(sum(v_in.^2) / sum(v_out_conv.^2));

f_dominant = 1 / (2*pi*sigma_t);
[~, f_idx] = min(abs(f - f_dominant));

SE_fft_dB = -20*log10(abs(H(f_idx)));

%% Diagnostics
fprintf('=== Poor Conductor / Concrete Shield Simulation ===\n');
fprintf('Conductivity = %.3e S/m\n', sigma_mat);
fprintf('Relative permittivity = %.2f\n', eps_r);
fprintf('Thickness = %.3e m\n', d);
fprintf('Z0 = %.2f Ohm\n', Z_0);
fprintf('Input peak = %.3e V\n', max(abs(v_in)));
fprintf('Impulse response peak = %.3e\n', max(abs(h)));
fprintf('Output max (conv) = %.3e V\n', max(v_out_conv));
fprintf('Output min (conv) = %.3e V\n', min(v_out_conv));
fprintf('Dominant frequency = %.3e Hz\n', f_dominant);
fprintf('SE from H(f) at f_dom = %.1f dB\n', SE_fft_dB);
fprintf('Time-domain peak SE = %.1f dB\n', SE_peak_td);
fprintf('Time-domain energy SE = %.1f dB\n', SE_energy_td);

%% Figure 1: Input pulse
figure;
plot(t*1e6, v_in, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Input Gaussian Pulse');
xlim([0 10]);

%% Figure 2: Improved impulse response visualization
figure;

subplot(3,1,1);
plot(t_h*1e6, h, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time [\mus]');
ylabel('h(t)');
title('Impulse Response: Full Linear Scale');
xlim([0 10]);

subplot(3,1,2);
plot(t_h*1e6, h, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time [\mus]');
ylabel('h(t)');
title('Impulse Response: Zoomed View');
xlim([0.05 5]);

subplot(3,1,3);
semilogy(t_h*1e6, h_norm, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time [\mus]');
ylabel('Normalized |h(t)|');
title('Normalized Impulse Response Magnitude');
xlim([0 5]);
ylim([1e-6 1]);

%% Figure 3: Input and output
figure;
plot(t*1e6, v_in, 'b', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2, 'DisplayName', 'Output');
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title(sprintf('Gaussian Pulse Through Concrete Shield, d = %.0f cm', d*100));
legend('Location', 'best');
xlim([0 10]);

%% Figure 4: Frequency-domain validation
figure;
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2, ...
     'DisplayName', 'Time-domain convolution');
hold on;
plot(t_fd*1e6, v_out_fd, 'b--', 'LineWidth', 2, ...
     'DisplayName', 'Frequency-domain multiplication');
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Validation: Time Domain vs Frequency Domain');
legend('Location', 'best');
xlim([0 10]);

%% Figure 5: Frequency-domain SE
figure;

SE_dB_plot = -20*log10(abs(H));
f_mask = f >= 1e4 & f <= 50e6;

semilogx(f(f_mask)/1e6, SE_dB_plot(f_mask), 'b', 'LineWidth', 1.5);
hold on;
semilogx(f_dominant/1e6, SE_fft_dB, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Frequency [MHz]');
ylabel('SE [dB]');
title(sprintf('SE(f): f_{dom}=%.2f MHz, SE=%.1f dB', ...
    f_dominant/1e6, SE_fft_dB));
legend('SE(f)', 'f_{dominant}', 'Location', 'best');

%% Figure 6: Thesis-style 4-subplot summary
figure('Position', [100, 100, 1400, 450]);

subplot(4,1,1);
plot(t*1e6, v_in, 'b', 'LineWidth', 1.5);
grid on;
ylabel('v_{in}(t) [V]');
title(sprintf('Gaussian Pulse Through Concrete Shield, d = %.0f cm', d*100));
xlim([0 10]);

subplot(4,1,2);
semilogy(t_h*1e6, h_norm, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Norm. |h(t)|');
title('Normalized Impulse Response Magnitude');
xlim([0 5]);
ylim([1e-6 1]);

subplot(4,1,3);
semilogx(f(f_mask)/1e6, SE_dB_plot(f_mask), 'b', 'LineWidth', 1.5);
hold on;
semilogx(f_dominant/1e6, SE_fft_dB, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
ylabel('SE [dB]');
xlabel('Frequency [MHz]');
title(sprintf('SE(f): at f_{dom}=%.2f MHz → %.1f dB', ...
    f_dominant/1e6, SE_fft_dB));

subplot(4,1,4);
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 1.5);
grid on;
ylabel('v_{out}(t) [V]');
xlabel('Time [\mus]');
xlim([0 10]);

figure('Position', [100, 100, 1400, 450]);

subplot(2,1,1);
plot(t*1e6, v_in, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Input: Gaussian Pulse');
xlim([0 10]);

subplot(2,1,2);
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Output: After 10cm Concrete Shield');
xlim([0 10]);

exportgraphics(gcf, 'C:\Users\Varun\Desktop\gaussian_concrete.png', 'Resolution', 150);