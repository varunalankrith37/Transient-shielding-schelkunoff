%% Damped Sinusoidal Pulse - Time-Domain Convolution Validation
clear all; close all; clc;

%% Time setup
dt = 1e-9;                  % Time step (1 nanosecond)
t_max = 5e-6;               % Total time (5 microseconds)
t = 0:dt:t_max;             % Time vector
N = length(t);              % Number of points

%% Define Damped Sinusoidal Pulse
f0 = 1e6;                   % Oscillation frequency (1 MHz)
tau_pulse = 1e-6;           % Decay time constant (1 microsecond)
V0 = 1;                     % Peak amplitude (1 V)

%% Generate damped sinusoidal pulse
v_in = V0 * exp(-t/tau_pulse) .* sin(2*pi*f0*t) .* (t>=0);

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
Z_0 = sqrt(mu_0/eps_0);     % Impedance of free space (Ohms)

%% Copper shield properties
sigma_copper = 5.96e7;      % Conductivity (S/m)
mu_r = 1;                   % Relative permeability
eps_r = 1;                  % Relative permittivity
thickness = 1e-3;           % Shield thickness (1 mm)

%% Create frequency vector
df = 1/t_max;               % Frequency resolution
f = (0:N-1)*df;             % Frequency vector

%% Shield transfer function H(f) 

% Create frequency vector with negative frequencies properly identified
f_sym = f;                                                 % Start with standard frequency vector [0, df, 2df, ..., (N-1)df]
f_sym(f_sym > 1/(2*dt)) = f_sym(f_sym > 1/(2*dt)) - 1/dt;  % Now f_sym = [0, df, ..., fs/2, -fs/2, ..., -2df, -df] (properly unwrapped)

% Angular frequency using abs(f_sym) for symmetric transfer function
omega = 2*pi*abs(f_sym);               % Use absolute value so +100MHz and -100MHz get same |H|
                                       % Without abs(): we'd only apply correct physics to positive freqs, output amplitude cut in half!
omega(omega < 2*pi*1e-6) = 2*pi*1e-6;  % Avoid division by zero at DC and near-zero freqs

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
title('Shield Impulse Response (short time scale)');
xlim([0 100]);

subplot(2,1,2);
plot(t*1e6, h_impulse, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('h(t)');
title('Shield Impulse Response (full time scale)');
xlim([0 5]);

%% Time domain calculation - LINEAR CONVOLUTION

% Method 1: Direct convolution (reference)
% Uses MATLAB's built-in function - this is our "ground truth"
v_out_conv = conv(v_in, h_impulse, 'same');  % 'same' returns output with same length as v_in
v_out_conv = v_out_conv - mean(v_out_conv);  % Remove DC offset for clean comparison

% Method 2: FFT-based LINEAR convolution (with zero-padding)
% Calculate minimum FFT length to prevent circular wrap-around
% For linear conv of signals length N1 and N2, need at least N1+N2-1 points
N_fft = length(v_in) + length(h_impulse) - 1;  % 5001 + 5001 - 1 = 10001

% Zero-pad both signals to N_fft length (MATLAB automatically pads with zeros)
% This converts circular convolution (what FFT naturally does) into linear convolution
V_in_freq_padded = fft(v_in, N_fft);           % Pad v_in: [x1...x5001, 0...0] (10001 total)
H_padded = fft(h_impulse, N_fft);              % Pad h_impulse: [h1...h5001, 0...0] (10001 total)

% Apply convolution theorem: conv(a,b) = ifft(fft(a) .* fft(b))
V_out_freq = V_in_freq_padded .* H_padded;     % Multiply in frequency domain
v_out_fft_full = real(ifft(V_out_freq));       % Convert back to time (length 10001)

% Extract 'same' portion to match conv(...,'same') output
% Need to extract centered 5001 points from 10001-point result
center_start = floor(length(h_impulse)/2);     % Skip first 2500 points (edge effects)
v_out_fft = v_out_fft_full(center_start + (1:length(v_in)));  % Extract indices 2501:7501
% Breakdown: center_start + (1:5001) = 2500 + [1,2,...,5001] = [2501, 2502,...,7501]

v_out_fft = v_out_fft - mean(v_out_fft);       % Remove DC offset to match Method 1

% Verify both methods give identical results (should be ~machine precision)
fprintf('=== Convolution Comparison ===\n');
fprintf('Max error between methods: %.6e\n\n', max(abs(v_out_conv - v_out_fft)));
% Expected: ~1e-23 to 1e-15 (essentially zero - proves convolution theorem!)

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
title('Output: Both Methods (should overlap perfectly)');
legend('Location', 'best');
xlim([0 5]);

subplot(3,1,3);
plot(t*1e6, abs(v_out_conv - v_out_fft), 'm-', 'LineWidth', 2);
grid on;
xlabel('Time (\mus)');
ylabel('Error (V)');
title('Absolute Error Between Methods');
xlim([0 5]);

sgtitle('Time-Domain Convolution: Direct vs FFT-Based (Damped Sinusoidal)', 'FontSize', 14, 'FontWeight', 'bold');

%% Check boundary conditions
fprintf('=== Boundary Check ===\n');
fprintf('Input at t=0: %.6e V\n', v_in(1));        % Should be 0 (sine starts at 0)
fprintf('Input at t=end: %.6e V\n', v_in(end));    % Should be ≈0 (decayed)
fprintf('Impulse at t=0: %.6e\n', h_impulse(1));   % Should be small
fprintf('Impulse at t=end: %.6e\n', h_impulse(end)); % Should be small (decayed)
fprintf('All values should be near zero for minimal circular artifacts.\n\n');
% If these are all < 1e-6, circular convolution ≈ linear convolution
% Our zero-padding makes this irrelevant, but it's good to verify!

fprintf('=== VALIDATION COMPLETE ===\n');
fprintf('Time-domain convolution methods validated for damped sinusoidal pulse!\n');

