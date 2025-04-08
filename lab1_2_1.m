clc;
close all;
% parámetros de configuración
A = 1; %Amplitud
fm = 100000; % Hz
tm = 1/fm; % segundos
ls = 200; % largo de la señal
f_c = 1000; % Hz
f_s = 5000; % Hz
t_s = 1/f_s; % segundos
tau = 0.5*t_s; % segundos
d = tau/t_s; % ciclo de trabajo

% vectores
t = (0:ls-1)*tm;
m_t = A*sin(2*pi*f_c*t);

% auxiliaries
r = floor(t_s/tm);
s = floor(tau/tm);
disp(r)

% muestreo natural
s_nat = zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
s_nat(i:i+s) = 1;
end
end
s_nat = s_nat(1:length(t));
m_t_nat = m_t.*s_nat;

% muestreo instantaneo
m_t_inst = zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
m_t_inst(i:i+s) = m_t(i);
end
end
m_t_inst = m_t_inst(1:length(t));

% Fourier Transform
M_t = fft(m_t);
M_t_nat = fft(m_t_nat);
M_t_inst = fft(m_t_inst);
% Frequency axis for plotting
f_axis = (0:(length(t) - 1)) * (1 / (ls * tm));
% Create the first figure with original, natural, and instantaneous signals
figure;
subplot(3, 1, 1);
plot(t, m_t);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3, 1, 2);
plot(t, m_t_nat);
title('Muestreo Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3, 1, 3);
plot(t, m_t_inst);
title('Muestreo Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
% Create the second figure with Fourier Transforms
figure;
subplot(3, 1, 1);
plot(f_axis, abs(M_t));
title('Transformada de Fourier Señal Original');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
subplot(3, 1, 2);
plot(f_axis, abs(M_t_nat));
title('Transformada de Fourier Muestreo Natural');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
subplot(3, 1, 3);
plot(f_axis, abs(M_t_inst));
title('Transformada de Fourier Muestreo Instantáneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
% PCM Parameters
bit_depth = 8; % Number of bits for PCM
pcm_levels = 2^bit_depth; % Total PCM levels
% Quantize the instantaneous signal using PCM
pcm_signal_inst = round((m_t_inst + 1) * (pcm_levels - 1) / 2); % Quantization
grid on;