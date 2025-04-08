clc;
close all;
A=1;
fm=100000;
tm=1/fm;
ls=200;
f_c=1000;
f_s=5000;
t_s=1/f_s;
tau=0.5*t_s;
d=tau/t_s; %

t=(0:ls-1)*tm;
m_t=A*sin(2*pi*f_c*t);
r=floor(t_s/tm);
s=floor(tau/tm);
disp(r)

s_nat=zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
s_nat(i:i+s)=1;
end
end
s_nat=s_nat(1:length(t));
m_t_nat=m_t.*s_nat;
% muestreo instantaneo

m_t_inst=zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
m_t_inst(i:i+s) = m_t(i);
end
end
m_t_inst=m_t_inst(1:length(t));
N =64; %cantidad de bits para pcm
pcm_levels=2ˆN; % Total de niveles PCM
pcm_signal_inst=round((m_t_inst+1)*(pcm_levels-1)/2); %cuantizacion
m_t_norm=(m_t-min(m_t))/(max(m_t)-min(m_t));
m_t_inst_norm=(m_t_inst-min(m_t_inst))/(max(m_t_inst)-min(m_t_inst));
pcm_signal_inst_norm = (pcm_signal_inst-min(pcm_signal_inst))/(max(pcm_signal_inst)-min(pcm_signal_inst));

%calculo del error de cuantizacion
quantization_error_inst=m_t_inst-((2*pcm_signal_inst/(pcm_levels-1))-1);

figure;
plot(t, m_t_norm, 'b', 'LineWidth', 1.5);
hold on;
plot(t, m_t_inst_norm, 'r', 'LineWidth', 1.5);
stem(t, pcm_signal_inst_norm, 'g', 'Marker','o', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud Normalizada');
title('señal original', 'señal PAM
instantanea', 'señal PAM Cuantificada(PCM)');
legend('señal original', 'señal PAM
instantanea', 'señal PAM Cuantificada(PCM)');
grid on;

figure;
plot(t, quantization_error_inst, 'k--',
'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Error de Cuantizaci n');
title('Error de Cuantizaci n para la S e a l
PAM Cuantificada (PCM)');
grid on;
