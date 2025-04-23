clc
close all
clear all

%% Globals

% Urban micro (LoS)
env_type = 'UMi_LoS';
% Frequency (in GHz)
freq_ghz = 3;
% Number of MPCs
num_mpcs = 100;

%% Generate RMS-DS
% RMS-DS statistics
[rmsds_mean, rmsds_std, r_prop_param] = get_delay_stats( freq_ghz, env_type );
% Generate random RMS-DS (Log scale)
rmsds_log = normrnd( rmsds_mean, rmsds_std, [1,1] );
% Convert to seconds 
rmsds = 10^( rmsds_log );


%% MPCs delays

% Average delay
mu_delay = r_prop_param * rmsds;
% Delay terms
tau_n_pp = exprnd( mu_delay, [num_mpcs, 1] );
tau_n_p = tau_n_pp - min( tau_n_pp );
tau_n = sort( tau_n_p );

% ------------------------------------------------
%% Plot PDF Delay
% ------------------------------------------------
% figure(1)
% num_delay_samples = 1e6;
% tau_n_pdf = exprnd( mu_delay, [num_delay_samples, 1] );
% tau_n_us = tau_n_pdf;

% % max(tau_n_pdf_pdf)
% delay = linspace(0, max(tau_n_pdf), 100) ;
% pdf_delay = ( 1 / mu_delay ) * exp( - delay ./ mu_delay );

% histogram(tau_n_us, 100, 'normalization', 'pdf')
% hold on
% plot(delay , pdf_delay, 'LineWidth', 2, 'Color', 'r')

% title('Histograma de Atrasos e PDF Teórica');
% xlabel('Atraso');
% ylabel('Dsitribuição');
% grid on

%% MPCs powers

% Power stats
[rf_mean, rf_std, sh_mean, sh_std] = get_power_stats( freq_ghz, env_type );

% Generate shadowing samples
xi_n_db = normrnd( sh_mean, sh_std, [num_mpcs, 1] );
xi_n = 10.^( -xi_n_db / 10.0 ); % Convert to linear

% Power terms
power_c = ( r_prop_param - 1 ) / mu_delay;
alpha2_n_p = exp( - tau_n * power_c ) .* xi_n;

% Generate a Rice factor sample
rf_db = normrnd( rf_mean, rf_std, [1,1] );
rf = 10.^( rf_db / 10 ); % Convert to linear

% Scattered power
omega_c_s = sum( alpha2_n_p(2:end) );

% Normalize power components based on Rice factor
alpha2_n_p(2:end) = ( 1 / ( rf + 1 ) ) * ( alpha2_n_p(2:end) / omega_c_s );
alpha2_n_p(1) = ( rf / ( rf + 1 ) );
% Final powers after normalization
alpha2_n = alpha2_n_p;

% Reverse computation of the Rice factor
rf_c = alpha2_n_p(1) / ( sum( alpha2_n_p( 2 : end ) ) );

% ------------------------------------------------
%% Plot PDP
% ------------------------------------------------
% figure(2)
% stem( tau_n / 1e-6, alpha2_n, '^', 'Color', 'k', 'Linewidth', 1.5 );
% hold on
% stem( tau_n(1) / 1e-6, alpha2_n(1), '^', 'Color', 'blue', 'Linewidth', 1.5 );

% ylim([1e-10, 1]);
% % xlim([0, 10]);

% xlabel( 'Atraso multipercurso -- $\tau$ ($\mu$s)', 'Interpreter', 'Latex', 'Fontsize', 13 );
% ylabel( 'Pot{\^{e}}ncia multipercurso', 'Interpreter', 'Latex', 'Fontsize', 13 );
% grid on

% ax = gca;
% ax.TickLabelInterpreter = 'Latex';
% ax.FontSize = 14;

% % Put y-scale on log
% set(ax,'yscal','log');

%% MPC Angles

% Angle stats
[azimuth_mean, azimuth_std, elevation_mean, elevation_std] = get_angular_stats( freq_ghz, env_type );

%% Azimuth angular scattering
sigma_theta = normrnd( azimuth_mean, azimuth_std, [1, 1] );
sigma_theta_rad = (pi / 180) * (10 ^ sigma_theta);

% Azimuth samples
alpha_max = max( alpha2_n );
theta_n_p = 1.42 * sigma_theta_rad * sqrt( -log( alpha2_n ./ alpha_max ) );

% Random azimuth varing samples
un = randsample( [-1, 1], num_mpcs, true)';
% Flutuacoes angulares
yn = normrnd(0, sigma_theta_rad/7, [num_mpcs, 1]);
% Final azimuth angles
theta_n = un .* theta_n_p + yn;
% Los
theta_n = theta_n - theta_n(1);

% ------------------------------------------------
%% Plot Azimuth Angles
% ------------------------------------------------
% figure(3)

% for n = 1 : length( theta_n )
% %for n = 1 : 5   

%     power_x = 10 * log10( alpha2_n(n) / min( alpha2_n ) );
%     if n == 1 
%         polarplot( [theta_n(n),theta_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'blue' );
%     else
%         polarplot( [theta_n(n),theta_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'red' );
%     end
    
%     hold on
% end
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = 2;
% ax.ThetaTickLabel = {'0'; '30'; '60'; '90'; '120'; '150'; '180'; '210'; '240'; '270'; '300'; '330';};
% title('Espalhamento angular em azimute', 'Interpreter', 'Latex');
% grid on

% figure(4)
% stem( rad2deg( theta_n ), alpha2_n, '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
% hold on
% stem( rad2deg( theta_n(1) ), alpha2_n(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue' );
% ylim([1e-8,1]);
% xlim([-6 * (10 .^ sigma_theta), 6 * (10 .^ sigma_theta)]);
% set(gca,'yscal','log');
% grid on
% xlabel('{\^{A}}ngulos de chegada em azimute ($^{\circ}$)', 'Interpreter', 'Latex');
% ylabel('Pot{\^{e}}ncia', 'Interpreter', 'Latex');
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% grid on

%% Elevation angular scattering
sigma_phi = normrnd( elevation_mean, elevation_std, [1, 1] );
sigma_phi_rad = (pi / 180) * (10 ^ sigma_phi);

% Elevation samples
phi_n_p = -sigma_phi_rad * log(alpha2_n ./ alpha_max);

% Random elevation variation
un = randsample( [-1, 1], num_mpcs, true)';
% Angular flutuation
yn = normrnd(0, sigma_theta_rad/7, [num_mpcs, 1]);
% Medium elevation
a = 0;
b = pi/2;
phi_bar = a + (b-a) .* rand(1, 1);

% Final elevation angles
phi_n = un .* phi_n_p + yn;
% Los
phi_n = phi_n - phi_n(1) + phi_bar;

% ------------------------------------------------
%% Plot Elevation Angles
% ------------------------------------------------
% figure(5)

% for n = 1 : length( phi_n )
% %for n = 1 : 5   

%     power_x = 10 * log10( alpha2_n(n) / min( alpha2_n ) );
%     if n == 1 
%         polarplot( [phi_n(n),phi_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'blue' );
%     else
%         polarplot( [phi_n(n),phi_n(n)], [0,power_x], '-o', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'red' );
%     end
    
%     hold on
% end

% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = 2;
% ax.ThetaTickLabel = {'0'; '30'; '60'; '90'; '120'; '150'; '180'; '210'; '240'; '270'; '300'; '330';};
% title('Espalhamento angular em elevacao', 'Interpreter', 'Latex');
% grid on

% figure(6)
% stem( rad2deg( phi_n ), alpha2_n, '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
% hold on
% stem( rad2deg( phi_n(1) ), alpha2_n(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue' );
% ylim([1e-8,1]);
% xlim([-6 * (10 .^ sigma_theta), 6 * (10 .^ sigma_theta)]);
% set(gca,'yscal','log');
% grid on
% % title()
% xlabel('{\^{A}}ngulos de chegada em Elevacao ($^{\circ}$)', 'Interpreter', 'Latex');
% ylabel('Pot{\^{e}}ncia', 'Interpreter', 'Latex');
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize = 14;
% grid on


% Angles vectors
rn = [cos(theta_n) .* sin(phi_n), sin(theta_n) .* sin(phi_n), cos(phi_n)];

%% Doppler efects
v_rx = 10; % m/s
c = 299792458;
lambda = c / (freq_ghz * (10^9));

a = 0;
b = 2*pi;
theta_v = a + (b-a) .* rand(1, 1); % number between 0 and 2*pi
a = 0;
b = pi;
phi_v = a + (b-a) .* rand(1, 1); % number between 0 and pi

vector_rx = v_rx * [cos(theta_v) * sin(phi_v), sin(theta_v) * sin(phi_v), cos(phi_v)];

% %% Plot angle vectors
% figure(7)
% for i = 1 : num_mpcs
%     if i == 1
%         style = "b-^";
%     else
%         style = "r-o";
%     end
%     plot3( [0 rn(i, 1)], [0 rn(i, 2)], [0 rn(i, 3)], style, 'LineWidth', 1.5 );
%     hold on;
% end
% xlabel('x'), ylabel('y'), zlabel('z')
% set(gca,'CameraPosition',[1 2 3]);
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize = 14;
% grid on
% style = "g->";
% plot3( [0 vector_rx(1, 1)]/v_rx, [0 vector_rx(1, 2)]/v_rx, [0 vector_rx(1, 3)]/v_rx, style, 'LineWidth', 1.5);
% blanck_labels = repmat({''}, 1, num_mpcs-2);
% legend(['LoS','NLoS', blanck_labels, 'Direção Velocidade'])

% % Doppler shift
vn = (1/lambda) .* (rn * vector_rx');
vn = sum(vn, 2);

%% Plot doppler shift
% figure(8)
% stem( vn(1), alpha2_n(1), '^', 'Linewidth', 2.0, 'Color', 'blue', 'MarkerFacecolor', 'blue');
% hold on
% stem( vn(2:num_mpcs) , alpha2_n(2:num_mpcs), '^', 'Linewidth', 1.0, 'Color', 'k', 'MarkerFacecolor', 'k' );
% ylim([1e-8,1]);
% % xlim([-6 * (10 .^ sigma_theta), 6 * (10 .^ sigma_theta)]);
% set(gca,'yscal','log');
% grid on
% title("$v_{rx}="+num2str(v_rx)+"$ m/s, $f_c = "+num2str(freq_ghz)+"$ GHz", 'Interpreter', 'Latex')
% xlabel('Desvio Doppler - $\nu$ (Hz)', 'Interpreter', 'Latex');
% ylabel('Potencia', 'Interpreter', 'Latex');
% legend('LoS', 'NLoS')
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize = 14;
% grid on

%% MPCs Phases
varphi_n = 2*pi* ((freq_ghz + vn) .* tau_n) - 2*pi*vn;






%% Signal transmitted
delta = 1e-7;
pulse_width = 1*delta;
num_samples = 1e5;
t = linspace(0, 5*delta, num_samples);
[signal_tx] = generate_pulse(0, t, pulse_width);

signal_rx = zeros(num_samples, num_mpcs);
for i = 1 : num_mpcs
    delayed_signal = generate_pulse(tau_n(i), t, pulse_width);
    signal_rx(:, i) = alpha2_n(i) * exp(-varphi_n(i) * 1j) * delayed_signal;
end

scattered_signal_rx = sum(signal_rx, 2);

figure(9)
plot(t, abs(signal_tx), 'Color', 'b', 'Linewidth', 1.5)
hold on
plot(t, abs(scattered_signal_rx), 'Color', 'r', 'Linewidth', 1.2)
xticks(-2*delta:delta:max(t))
xticklabels({'$-2\delta_t$', '$-\delta_t$', '$0$', '$\delta_t$', '$2\delta_t$', '$3\delta_t$', '$4\delta_t$', '$5\delta_t$', '$6\delta_t$'})
title('Sinal Transmitido', 'Interpreter', 'Latex')
ylabel('$s(t)$', 'Interpreter', 'Latex')
xlabel('$t$ [s]', 'Interpreter', 'Latex')
title("$\delta=10^{-7} s, \sigma_{\tau}="+num2str(rmsds * 1e9)+"ns$", 'Interpreter', 'Latex')
legend('TX', 'RX')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
grid on








% ---------------------------------------
%% Auxiliar Functions
% ---------------------------------------
function [rmsds_mean, rmsds_std, r_prop_param] = get_delay_stats( freq_ghz, env_type )
    switch env_type
        case 'UMi_LoS'
            rmsds_mean = -0.24 * log10( 1 + freq_ghz ) - 7.14;
            rmsds_std = 0.38;
            r_prop_param = 3;

        case 'UMi_NLoS'
            rmsds_mean = -0.24 * log10( 1 + freq_ghz ) - 6.83;
            rmsds_std = -0.16 * log10( 1 + freq_ghz ) + 0.28;
            r_prop_param = 2.1;
            
        case 'UMa_LoS'
            rmsds_mean = -0.0963 * log10( 1 + freq_ghz ) - 6.955;
            rmsds_std = 0.66;
            r_prop_param = 2.5;
            
        case 'UMa_NLoS'
            rmsds_mean = -0.204 * log10( 1 + freq_ghz ) - 6.28;
            rmsds_std = 0.39;
            r_prop_param = 2.3;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [rf_mean, rf_std, sh_mean, sh_std] = get_power_stats( freq_ghz, env_type )
    switch env_type
        case 'UMi_LoS'
            rf_mean = 9;
            rf_std = 5;
            sh_mean = 0;
            sh_std = 4;
            
        case 'UMi_NLoS'
            rf_mean = -inf;
            rf_std = -inf;
            sh_mean = 0;
            sh_std = 7.82;
            
        case 'UMa_LoS'
            rf_mean = 9;
            rf_std = 3.5;
            sh_mean = 0;
            sh_std = 4;
            
        case 'UMa_NLoS'
            rf_mean = -inf;
            rf_std = -inf;
            sh_mean = 0;
            sh_std = 6;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [azimuth_mean, azimuth_std, elevation_mean, elevation_std] = get_angular_stats( freq_ghz, env_type )
    switch env_type
        case 'UMi_LoS'
            azimuth_mean = -0.08 * log10( 1 + freq_ghz ) + 1.73;
            azimuth_std = 0.014 * log10( 1 + freq_ghz ) + 0.28;
            elevation_mean = -0.1 * log10( 1 + freq_ghz ) + 0.73;
            elevation_std = -0.04 * log10( 1 + freq_ghz ) + 0.34;

        case 'UMi_NLoS'
            azimuth_mean = -0.08 * log10( 1 + freq_ghz ) + 1.81;
            azimuth_std = 0.05 * log10( 1 + freq_ghz ) + 0.3;
            elevation_mean = -0.04 * log10( 1 + freq_ghz ) + 0.92;
            elevation_std = -0.07 * log10( 1 + freq_ghz ) + 0.41;
            
        case 'UMa_LoS'
            azimuth_mean = 1.81;
            azimuth_std = 0.2;
            elevation_mean = 0.95;
            elevation_std = 0.16;
            
        case 'UMa_NLoS'
            azimuth_mean = -0.27 * log10( freq_ghz ) + 2.08;
            azimuth_std = 0.11;
            elevation_mean = -0.3236 * log10( freq_ghz ) + 1.512;
            elevation_std = 0.16;
            
        otherwise
            error( 'Ambiente inválido.' )
    end
end

function [signal] = generate_pulse( delay, t, pulse_width )
    signal = zeros(length(t), 1);
    for i = 1 : length(t)
        if t(i) >= delay && t(i) <= pulse_width + delay
            signal(i) = 1;
        end
    end
end