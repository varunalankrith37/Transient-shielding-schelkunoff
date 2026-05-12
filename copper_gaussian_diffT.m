clear all; close all; clc;

%% Shield thickness 
d = 2e-3;  % [m]  0.5e-3, 1e-3, 2e-3

%% Time setup 
dt = 1e-8;
t_max = max(100e-6, (d/1e-3)^2 * 100e-6);
t = 0:dt:t_max;

%% Gaussian pulse
V0 = 1;
t0 = 2e-6;
sigma_t = 0.3e-6;
v_in = V0 * exp(-((t - t0).^2) / (2*sigma_t^2));

%% Fourier transform of input 
[f, V_freq] = fourier(t, v_in);

%% Constants
mu_0 = 4*pi*1e-7;
eps_0 = 8.854187817e-12;
Z_0 = sqrt(mu_0/eps_0);

%% Copper shield properties
sigma_cu = 5.8e7;
mu_r = 1;
eps_r = 1;
mu = mu_r * mu_0;
eps = eps_r * eps_0;

%% Schelkunoff transfer function
omega = 2*pi*f;
omega(1) = 2*pi*1e-6;

gamma = sqrt(1j*omega*mu .* (sigma_cu + 1j*omega*eps));
Z_s = sqrt(1j*omega*mu ./ (sigma_cu + 1j*omega*eps));
rho_01 = (Z_s - Z_0) ./ (Z_s + Z_0);
tau_01 = 2*Z_s ./ (Z_s + Z_0);
tau_10 = 2*Z_0 ./ (Z_s + Z_0);
A = exp(-gamma*d);

H = (tau_01 .* tau_10 .* A) ./ (1 - rho_01.^2 .* A.^2);
H(1) = 0;

%% Impulse response via invfourier + tail correction
[t_h, h_full] = invfourier(f, H);
h = real(h_full(1:length(t)));
t_h = t_h(1:length(t));

tail_start = round(0.8 * length(h));
h = h - mean(h(tail_start:end));

%% conv() output
v_out_conv_full = conv(v_in, h) * dt;
t_conv = (0:length(v_out_conv_full)-1) * dt;
v_out_conv = v_out_conv_full(1:length(t));
t_out = t_conv(1:length(t));

%% SE values
f_dominant = 1 / (2*pi*sigma_t);
[~, f_idx] = min(abs(f - f_dominant));
SE_fft_dB = -20*log10(abs(H(f_idx)));

%% Dynamic plot window
t_show = min((d/1e-3)^2 * 40e-6 * 1e6, t_max*1e6);  % in µs

%% Diagnostics
fprintf('=== Copper Shield Simulation ===\n');
fprintf('Conductivity = %.3e S/m\n', sigma_cu);
fprintf('Thickness = %.3e m\n', d);
fprintf('t_max = %.0f us\n', t_max*1e6);
fprintf('Z0 = %.2f Ohm\n', Z_0);
fprintf('Input peak = %.3e V\n', max(abs(v_in)));
fprintf('Impulse response peak = %.3e\n', max(abs(h)));
fprintf('h(t) tail mean (last 20%%) = %.3e\n', mean(h(tail_start:end)));
fprintf('Output max (conv) = %.3e V\n', max(v_out_conv));
fprintf('Output min (conv) = %.3e V\n', min(v_out_conv));
fprintf('Dominant frequency = %.3e Hz\n', f_dominant);
fprintf('SE from H(f) at f_dom = %.1f dB\n', SE_fft_dB);

%% Figure 1: Input Gaussian Pulse
figure;
plot(t*1e6, v_in, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Input Gaussian Pulse');
xlim([0 10]);  % input pulse is always in first 10 µs

%% Figure 2: Impulse Response
figure;
plot(t_h*1e6, h, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('h(t)');
title(sprintf('Impulse Response of Copper Shield (d=%.1f mm)', d*1e3));
xlim([0 t_show]);

%% Figure 3: Input vs Output
figure;
plot(t*1e6, v_in, 'b', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2, 'DisplayName', 'Output (conv)');
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title(sprintf('Gaussian Pulse Through Copper Shield (d=%.1f mm)', d*1e3));
legend('Location', 'best');
xlim([0 t_show]);

%% Figure 4: Zoomed output
figure;
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title(sprintf('Output Pulse After Copper Shield (d=%.1f mm)', d*1e3));
xlim([0 t_show]);

%% Figure 5: 4-subplot figure
figure;

subplot(4,1,1);
plot(t*1e6, v_in, 'b', 'LineWidth', 1.5);
grid on;
ylabel('v_{in}(t) [V]');
title(sprintf('Validation: Gaussian Pulse Through Copper Shield (d=%.1f mm)', d*1e3));
xlim([0 t_show]);

subplot(4,1,2);
plot(t_h*1e6, h, 'r', 'LineWidth', 1.5);
grid on;
ylabel('h(t)');
xlim([0 t_show]);

subplot(4,1,3);
SE_dB_plot = -20*log10(abs(H));
f_mask = f >= 1e4 & f <= 50e6;
semilogx(f(f_mask)/1e6, SE_dB_plot(f_mask), 'b', 'LineWidth', 1.5);
hold on;
semilogx(f_dominant/1e6, SE_fft_dB, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
ylabel('SE [dB]');
ylim([0 400]);
xlim([0.01 50]);
xlabel('Frequency [MHz]');
title(sprintf('SE(f): at f_{dom}=%.2f MHz → %.1f dB', ...
    f_dominant/1e6, SE_fft_dB));
legend('SE(f)', 'f_{dominant}', 'Location', 'northwest');

subplot(4,1,4);
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 1.5);
grid on;
ylabel('v_{out}(t) [V]');
xlabel('Time [\mus]');
xlim([0 t_show]);

%% Input vs Output 
figure('Position', [100, 100, 1400, 450]);

subplot(2,1,1);
plot(t*1e6, v_in, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude (V)');
title('Input: Gaussian Pulse');
xlim([0 60]);  % <-- changed from [0 10]

subplot(2,1,2);
plot(t_out*1e6, v_out_conv, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude (V)');
title('Output: After 5mm Copper Shield');
xlim([0 t_show]);

%sgtitle('Schelkunoff Shielding: Output vs Input (Time Domain)', 'FontSize', 14, 'FontWeight', 'bold');

figure('Position', [100, 100, 1400, 450]);

subplot(2,1,1);
plot(t*1e6, v_in, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude [V]');
title('Input: Gaussian Pulse');
xlim([0 60]);

subplot(2,1,2);
plot(t_h*1e6, h, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [\mus]');
ylabel('h(t)');
title('Impulse Response: Copper Shield (d = 2mm)');
xlim([0 60]);

exportgraphics(gcf, 'C:\Users\Varun\Desktop\copper_impulse2mm.png', 'Resolution', 150);

SE_peak_td = 20*log10(max(abs(v_in)) / max(abs(v_out_conv)));
SE_energy_td = 10*log10(sum(v_in.^2) / sum(v_out_conv.^2));

fprintf('Time-domain peak SE = %.1f dB\n', SE_peak_td);
fprintf('Time-domain energy SE = %.1f dB\n', SE_energy_td);