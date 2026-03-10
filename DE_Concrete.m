%% Damped Sinusoidal Pulse - Concrete Shield Analysis
clear all; close all; clc;

%% Time setup
dt = 1e-9;                  % Time step (1 nanosecond)
t_max = 5e-6;               % Total time (5 microseconds)
t = 0:dt:t_max;             % Time vector
N = length(t);              % Number of points

%% Define Double-Exponential Pulse (Lightning-Type)
alpha1 = 1.5e7;             % Fast decay rate (1/s)
alpha2 = 3e6;               % Slow decay rate (1/s)
V0 = 1;                     % Peak amplitude (1 V)

% Generate pulse
v_in = V0 * (exp(-alpha2*t) - exp(-alpha1*t)) .* (t>=0);

%% Plot input signal
figure;
plot(t*1e6, v_in, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input: Damped Sinusoidal Pulse');
xlim([0 5]);
fprintf('Input peak: %.6f V\n\n', max(abs(v_in)));

%% Physical constants
mu_0 = 4*pi*1e-7;           % Permeability of free space (H/m)
eps_0 = 8.854e-12;          % Permittivity of free space (F/m)
Z_0 = sqrt(mu_0/eps_0);     % Impedance of free space (~377 Ohms)


%% Concrete shield properties
sigma_concrete = 0.01;      % Conductivity (S/m)
mu_r = 1;                   % Relative permeability
eps_r = 6;                  % Relative permittivity (typical for concrete)
thickness = 0.2;            % Shield thickness (200 mm)

fprintf('=== Shield Material: CONCRETE ===\n');
fprintf('Conductivity: %.3f S/m (vs Concrete: 5.96e7 S/m)\n', sigma_concrete);
fprintf('Relative permittivity: %.1f (vs Copper: 1)\n', eps_r);
fprintf('Thickness: %.0f mm (vs typical metal shield: 1 mm)\n\n', thickness*1000);

%% Create frequency vector
df = 1/t_max;               % Frequency resolution
f = (0:N-1)*df;             % Frequency vector

%% Shield transfer function H(f) 

% Create frequency vector with negative frequencies properly identified
f_sym = f;
f_sym(f_sym > 1/(2*dt)) = f_sym(f_sym > 1/(2*dt)) - 1/dt;

% Angular frequency using abs(f_sym) for symmetric transfer function
omega = 2*pi*abs(f_sym);
omega(omega < 2*pi*1e-6) = 2*pi*1e-6;  % Avoid division by zero at DC

% Material properties (CONCRETE)
mu = mu_r * mu_0;              % Absolute permeability (H/m)
eps = eps_r * eps_0;           % Absolute permittivity (F/m) - HIGHER than copper

% Complex propagation constant
% For concrete: displacement current (jωε) matters more than conduction (σ)!
gamma = sqrt(1j*omega*mu.*(sigma_concrete + 1j*omega*eps));

% Intrinsic impedance of concrete
Z_s = sqrt(1j*omega*mu./(sigma_concrete + 1j*omega*eps));

% Schelkunoff's shielding theory components
rho_01 = (Z_s - Z_0)./(Z_s + Z_0);  % Reflection coefficient
tau_01 = 2*Z_s./(Z_s + Z_0);        % Air to concrete
tau_10 = 2*Z_0./(Z_s + Z_0);        % Concrete to air
A = exp(-gamma*thickness);          % Attenuation through 20 cm concrete

% Complete transfer function
H = (tau_01 .* tau_10 .* A) ./ (1 - rho_01.^2 .* A.^2);

% Diagnostic output
fprintf('=== Physical Parameters ===\n');
fprintf('Z_0 (free space) = %.2f Ohms\n', Z_0);
fprintf('Z_s (concrete, 1 MHz) = %.2f Ohms\n', abs(Z_s(100)));
fprintf('Max |H| = %.6e (at low frequencies)\n', max(abs(H)));
fprintf('Min |H| = %.6e (at high frequencies)\n', min(abs(H)));
fprintf('\nNote: Concrete Z_s >> Copper Z_s (concrete is poor conductor)\n\n');

%% Calculate impulse response (time domain)
h_impulse = real(ifft(H));
h_impulse = h_impulse - mean(h_impulse);    % Remove DC offset

% Check impulse response location
[peak_value, peak_idx] = max(abs(h_impulse));
fprintf('Impulse response peak at index: %d (time = %.3e s)\n', peak_idx, t(peak_idx));
fprintf('Peak value: %.6e\n\n', peak_value);

% Plot impulse response
figure;
subplot(2,1,1);
plot(t*1e9, h_impulse, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (ns)');
ylabel('h(t)');
title('Concrete Shield Impulse Response (short time scale)');
xlim([0 100]);

subplot(2,1,2);
plot(t*1e6, h_impulse, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('h(t)');
title('Concrete Shield Impulse Response (full time scale)');
xlim([0 5]);

%% Time domain calculation - LINEAR CONVOLUTION

% Method 1: Direct convolution
v_out_conv = conv(v_in, h_impulse, 'same');
v_out_conv = v_out_conv - mean(v_out_conv);

% Method 2: FFT-based LINEAR convolution (with zero-padding)
N_fft = length(v_in) + length(h_impulse) - 1;
V_in_freq_padded = fft(v_in, N_fft);
H_padded = fft(h_impulse, N_fft);
V_out_freq = V_in_freq_padded .* H_padded;
v_out_fft_full = real(ifft(V_out_freq));

% Extract 'same' portion
center_start = floor(length(h_impulse)/2);
v_out_fft = v_out_fft_full(center_start + (1:length(v_in)));
v_out_fft = v_out_fft - mean(v_out_fft);

% Verify both methods match
fprintf('=== Convolution Comparison ===\n');
fprintf('Max error between methods: %.6e\n\n', max(abs(v_out_conv - v_out_fft)));

% Output statistics
fprintf('=== Shielding Performance ===\n');
fprintf('Input peak: %.6f V\n', max(abs(v_in)));
fprintf('Output peak: %.6e V\n', max(abs(v_out_conv)));
fprintf('Attenuation factor: %.2e\n', max(abs(v_in))/max(abs(v_out_conv)));
fprintf('Attenuation (dB): %.1f dB\n\n', 20*log10(max(abs(v_in))/max(abs(v_out_conv))));

%% Plot comparison
figure('Position', [100 100 1200 800]);

subplot(3,1,1);
plot(t*1e6, v_in, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Input: Damped Sinusoidal Pulse');
xlim([0 5]);

subplot(3,1,2);
plot(t*1e6, v_out_conv, 'b-', 'LineWidth', 2.5, 'DisplayName', 'conv() method');
hold on;
plot(t*1e6, v_out_fft, 'g--', 'LineWidth', 2, 'DisplayName', 'FFT method (padded)');
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude (V)');
title('Output: Both Methods (20 cm Concrete Shield)');
legend('Location', 'best');
xlim([0 5]);

subplot(3,1,3);
plot(t*1e6, abs(v_out_conv - v_out_fft), 'm-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Error (V)');
title('Absolute Error Between Methods');
xlim([0 5]);

sgtitle('Concrete Shield: Time-Domain Convolution Validation', 'FontSize', 14, 'FontWeight', 'bold');

%% Plot Transfer Function
figure('Position', [100 100 1000 600]);
subplot(2,1,1);

% Plot only positive frequencies
idx_plot = 2:floor(N/2);
f_plot = f(idx_plot);
H_plot = abs(H(idx_plot));

semilogx(f_plot/1e6, H_plot, 'b-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('|H(f)|');
title('Concrete Shield: Transfer Function Magnitude');
xlim([0.01 1000]);

subplot(2,1,2);
SE_dB = -20*log10(H_plot);
semilogx(f_plot/1e6, SE_dB, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (MHz)');
ylabel('Shielding Effectiveness (dB)');
title('Concrete Shield: Shielding Effectiveness (20 cm thickness)');
xlim([0.01 1000]);
ylim([0 100]);

sgtitle('Concrete Shield Performance', 'FontSize', 14, 'FontWeight', 'bold');

%% Boundary check
fprintf('=== Boundary Check ===\n');
fprintf('Input at t=0: %.6e V\n', v_in(1));
fprintf('Input at t=end: %.6e V\n', v_in(end));
fprintf('Impulse at t=0: %.6e\n', h_impulse(1));
fprintf('Impulse at t=end: %.6e\n', h_impulse(end));

fprintf('\n=== VALIDATION COMPLETE ===\n');
fprintf('Concrete shield analysis complete!\n');
fprintf('Note: Concrete provides MUCH LESS shielding than copper.\n');

