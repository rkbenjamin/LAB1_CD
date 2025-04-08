clc;
close all;
clear all;
% parámetros de configuración
A = 1; % Amplitud
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
m_t = sin(2*pi*f_c*t);

% auxiliaries
r = floor(t_s/tm);
s = floor(tau/tm);
disp(r)

% muestreo instantaneo
m_t_inst = zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
m_t_inst(i:i+s) = m_t(i);
end
end
m_t_inst = m_t_inst(1:length(t));

plot(t,m_t);
hold on;
%plot(t,r)
%plot(t,m_t_nat,'-r');
plot(t,m_t_inst);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal PAM con Muestreo instantaneo');
grid on;