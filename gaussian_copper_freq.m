%% Multi-Material Shielding Comparison - COMPLETE CODE
% Compares Copper, Aluminum, Soil, and Concrete shields
clear all; close all; clc;

%% Time setup
dt = 1e-8;              % Time step (10 nanoseconds)
t_max = 100e-6;         % Total time (100 microseconds)
t = (0:dt:t_max)';      % Time vector (column)

%% Gaussian pulse parameters
V0 = 1;                 % Peak amplitude (1 Volt)
t0 = 2e-6;              % Center time (2 microseconds)
sigma = 0.3e-6;         % Standard deviation (0.3 microseconds)

%% Generate Gaussian pulse
v_in = V0 * exp(-((t - t0).^2) / (2 * sigma^2));

%% Physical constants
mu_0 = 4*pi*1e-7;           % Permeability of free space (H/m)
eps_0 = 8.854e-12;          % Permittivity of free space (F/m)
Z_0 = sqrt(mu_0/eps_0);     % Impedance of free space (~377 Ohms)

%% Define Material Properties
materials = struct();

% 1. COPPER
materials(1).name = 'Copper';
materials(1).sigma = 5.96e7;    % Conductivity (S/m)
materials(1).mu_r = 1;          % Relative permeability
materials(1).eps_r = 1;         % Relative permittivity
materials(1).thickness = 0.5e-3; % 0.5 mm
materials(1).color = 'r';

% 2. ALUMINUM
materials(2).name = 'Aluminum';
materials(2).sigma = 3.77e7;    % Conductivity (S/m) - ~63% of copper
materials(2).mu_r = 1;          % Relative permeability (non-magnetic)
materials(2).eps_r = 1;         % Relative permittivity
materials(2).thickness = 0.5e-3; % 0.5 mm
materials(2).color = 'b';

% 3. SOIL (Dry)
materials(3).name = 'Soil (Dry)';
materials(3).sigma = 0.001;     % Conductivity (S/m) - typical for dry soil
materials(3).mu_r = 1;          % Relative permeability
materials(3).eps_r = 4;         % Relative permittivity (typical range: 2-6)
materials(3).thickness = 0.2;   % 20 cm (thick layer)
materials(3).color = [0.6 0.4 0.2]; % Brown

% 4. CONCRETE
materials(4).name = 'Concrete';
materials(4).sigma = 0.01;      % Conductivity (S/m) - typical concrete
materials(4).mu_r = 1;          % Relative permeability
materials(4).eps_r = 6;         % Relative permittivity (typical: 4-10)
materials(4).thickness = 0.2;   % 20 cm
materials(4).color = [0.5 0.5 0.5]; % Gray

%% Convert input to frequency domain
[f, V_freq] = fourier(t, v_in, 'pulse');

fprintf('=== INPUT SIGNAL ===\n');
fprintf('Peak amplitude: %.3f V\n', max(v_in));
fprintf('Time points: %d\n', length(t));
fprintf('Frequency points: %d\n\n', length(f));

%% Process each material
n_materials = length(materials);
outputs = cell(n_materials, 1);
t_outputs = cell(n_materials, 1);
H_all = zeros(length(f), n_materials);

for m = 1:n_materials
    fprintf('=== PROCESSING: %s ===\n', materials(m).name);
    fprintf('Conductivity: %.3e S/m\n', materials(m).sigma);
    fprintf('Thickness: %.0f mm\n', materials(m).thickness * 1000);
    
    % Material properties
    sigma_mat = materials(m).sigma;
    mu_r = materials(m).mu_r;
    eps_r = materials(m).eps_r;
    thickness = materials(m).thickness;
    
    % Angular frequency
    omega = 2*pi*f;
    omega(abs(omega) < 2*pi*1e-6) = 2*pi*1e-6;  % Avoid division by zero
    
    % Material constants
    mu = mu_r * mu_0;
    eps = eps_r * eps_0;
    
    % Schelkunoff transfer function
    gamma = sqrt(1j*omega*mu.*(sigma_mat + 1j*omega*eps));
    Z_s = sqrt(1j*omega*mu./(sigma_mat + 1j*omega*eps));
    
    rho_01 = (Z_s - Z_0)./(Z_s + Z_0);
    tau_01 = 2*Z_s./(Z_s + Z_0);
    tau_10 = 2*Z_0./(Z_s + Z_0);
    A = exp(-gamma*thickness);
    
    H = (tau_01 .* tau_10 .* A) ./ (1 - rho_01.^2 .* A.^2);
    H_all(:, m) = H;
    
    % Apply shielding
    V_out_freq = V_freq .* H;
    
    % Convert back to time domain
    [t_out, v_out] = invfourier(f, V_out_freq, 'pulse');
    
    % Store results
    outputs{m} = v_out;
    t_outputs{m} = t_out;
    
    % Statistics
    fprintf('Output peak: %.6e V\n', max(abs(v_out)));
    fprintf('Attenuation factor: %.3e\n', max(abs(v_in))/max(abs(v_out)));
    fprintf('Attenuation (dB): %.1f dB\n', 20*log10(max(abs(v_in))/max(abs(v_out))));
    fprintf('Max |H| at 1 MHz: %.6e\n', abs(H(min(find(f >= 1e6)))));
    fprintf('\n');
end

%% Summary Table
fprintf('================================================================================\n');
fprintf('                              SUMMARY TABLE\n');
fprintf('================================================================================\n');
fprintf('%-15s | Thickness | Input (V) | Output (V) | Atten (dB)\n', 'Material');
fprintf('--------------------------------------------------------------------------------\n');

for m = 1:n_materials
    fprintf('%-15s | %5.0f mm  | %9.3f | %10.3e | %10.1f\n', ...
        materials(m).name, ...
        materials(m).thickness * 1000, ...
        max(abs(v_in)), ...
        max(abs(outputs{m})), ...
        20*log10(max(abs(v_in))/max(abs(outputs{m}))));
end
fprintf('================================================================================\n\n');

% FIGURES - FULLY UNIFORM (LaTeX compatible)

% Set default font sizes for all figures
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLegendFontSize', 10);

% UNIFORM FIGURE SIZE FOR ALL PLOTS
fig_width = 1000;
fig_height = 450;

%% Figure 1: All Outputs Overlaid
figure('Position', [100 100 fig_width fig_height]);

plot(t*1e6, v_in, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Input');
hold on;

for m = 1:n_materials
    plot(t_outputs{m}*1e6, outputs{m}, ...
        'Color', materials(m).color, 'LineWidth', 2, ...
        'DisplayName', materials(m).name);
end

grid on;
xlabel('Time (\mus)', 'FontSize', 12);
ylabel('Amplitude (V)', 'FontSize', 12);
title('Shielding Comparison: All Materials', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
xlim([0 10]);
set(gca, 'FontSize', 11);

text(0.98, 0.02, 'Note: Metal outputs (Cu, Al) not visible at this scale (~10^{-8} V)', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'FontSize', 9, 'BackgroundColor', 'white');

%% Figure 2: Metals vs Dielectrics 
figure('Position', [100 100 fig_width fig_height]);

% Left: Metals (LOG scale)
subplot(1,2,1);
semilogy(t*1e6, v_in, 'k--', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
for m = 1:2
    semilogy(t_outputs{m}*1e6, abs(outputs{m}), ...
        'Color', materials(m).color, 'LineWidth', 2.5, ...
        'DisplayName', materials(m).name);
end
grid on;
xlabel('Time (\mus)', 'FontSize', 12);
ylabel('Amplitude (V)', 'FontSize', 12);  % Shortened label
title('Good Conductors (Metals)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
xlim([0 10]);
ylim([1e-10 2]);
set(gca, 'FontSize', 11, 'YScale', 'log');  % Explicitly set log scale

% Annotation box with peak values
text(0.5, 0.45, sprintf('Cu: ~%.1e V\nAl: ~%.1e V', ...
    max(abs(outputs{1})), max(abs(outputs{2}))), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', 'white', 'EdgeColor', 'k', ...
    'HorizontalAlignment', 'center');

% Right: Dielectrics (LINEAR scale)
subplot(1,2,2);
plot(t*1e6, v_in, 'k--', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
for m = 3:4
    plot(t_outputs{m}*1e6, outputs{m}, ...
        'Color', materials(m).color, 'LineWidth', 2.5, ...
        'DisplayName', materials(m).name);
end
grid on;
xlabel('Time (\mus)', 'FontSize', 12);
ylabel('Amplitude (V)', 'FontSize', 12);
title('Poor Conductors (Dielectrics)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
xlim([0 10]);
ylim([0 1.05]);
set(gca, 'FontSize', 11);

% Annotation box with peak values
text(0.5, 0.55, sprintf('Soil: ~%.2f V\nConcrete: ~%.2f V', ...
    max(abs(outputs{3})), max(abs(outputs{4}))), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', 'white', 'EdgeColor', 'k', ...
    'HorizontalAlignment', 'center');

% Add note about log scale in the title or as text
subplot(1,2,1);
text(0.02, 0.98, 'Log Scale', 'Units', 'normalized', ...
    'FontSize', 9, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'yellow', 'EdgeColor', 'k');

%% Figure 3: Transfer Functions
figure('Position', [100 100 fig_width fig_height]);

% Left: Transfer Function Magnitude
subplot(1,2,1);
for m = 1:n_materials
    loglog(f/1e6, abs(H_all(:,m)), ...
        'Color', materials(m).color, 'LineWidth', 2.5, ...
        'DisplayName', materials(m).name);
    hold on;
end
grid on;
xlabel('Frequency (MHz)', 'FontSize', 12);
ylabel('|H(f)|', 'FontSize', 12);
title('Transfer Function Magnitude', 'FontSize', 13);
legend('Location', 'southwest', 'FontSize', 10);
xlim([0.01 100]);
ylim([1e-12 1]);
set(gca, 'FontSize', 11);

% Right: Shielding Effectiveness
subplot(1,2,2);
for m = 1:n_materials
    SE = -20*log10(abs(H_all(:,m)));
    semilogx(f/1e6, SE, ...
        'Color', materials(m).color, 'LineWidth', 2.5, ...
        'DisplayName', materials(m).name);
    hold on;
end
grid on;
xlabel('Frequency (MHz)', 'FontSize', 12);
ylabel('Shielding Effectiveness (dB)', 'FontSize', 12);
title('Shielding Effectiveness', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 10);
xlim([0.01 100]);
ylim([0 250]);
set(gca, 'FontSize', 11);

%% Figure 4: Individual Materials 
figure('Position', [100 100 fig_width 650]);

titles_list = {'Copper (0.5 mm)', 'Aluminum (0.5 mm)', 'Soil (200 mm)', 'Concrete (200 mm)'};

for m = 1:4
    subplot(2,2,m);
    
    if m <= 2  % Metals - use log scale
        semilogy(t*1e6, v_in, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Input');
        hold on;
        semilogy(t_outputs{m}*1e6, abs(outputs{m}), ...
            'Color', materials(m).color, 'LineWidth', 2.5, ...
            'DisplayName', materials(m).name);
        ylabel('Amplitude (V) - Log', 'FontSize', 11);
        ylim([1e-10 2]);  % Adjusted to match your data
    else  % Dielectrics - use linear scale
        plot(t*1e6, v_in, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Input');
        hold on;
        plot(t_outputs{m}*1e6, outputs{m}, ...
            'Color', materials(m).color, 'LineWidth', 2.5, ...
            'DisplayName', materials(m).name);
        ylabel('Amplitude (V)', 'FontSize', 11);
        ylim([0 1.05]);
    end
    
    grid on;
    xlabel('Time (\mus)', 'FontSize', 11);
    title(titles_list{m}, 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 9);
    xlim([0 10]);
    set(gca, 'FontSize', 10);
    
    % Peak value and attenuation annotation
    peak_val = max(abs(outputs{m}));
    atten_dB = 20*log10(max(abs(v_in))/peak_val);
    text(0.98, 0.05, sprintf('Peak: %.2e V\nAtten: %.0f dB', peak_val, atten_dB), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'FontSize', 9, ...
        'BackgroundColor', 'white', 'EdgeColor', 'k');
end

%% Print summary
fprintf('\n=== FIGURE GENERATION COMPLETE ===\n');
fprintf('All figures generated with uniform formatting:\n');
fprintf('  - Figures 1-3: %d x %d pixels (1×2 layout)\n', fig_width, fig_height);
fprintf('  - Figure 4:    %d x 650 pixels (2×2 layout)\n', fig_width);
fprintf('\nRecommended LaTeX usage:\n');
fprintf('  \\includegraphics[width=1.0\\textwidth]{figure_name.png}\n');