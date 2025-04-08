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
d=tau/t_s;

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

m_t_inst=zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
m_t_inst(i:i+s)=m_t(i);
end
end
m_t_inst=m_t_inst(1:length(t));

M_t=fft(m_t);
M_t_nat=fft(m_t_nat);
M_t_inst=fft(m_t_inst);
f_axis=(0:(length(t)-1))*(1/(ls*tm));

figure;
subplot(3, 1, 1);
plot(t, m_t);
title('Señal original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3, 1, 2);
plot(t, m_t_nat);
title('Muestreo natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3, 1, 3);
plot(t, m_t_inst);
title('muestreo instantaneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
figure;
subplot(3, 1, 1);
plot(f_axis, abs(M_t));
title('Transformada de Fourier señal original');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
subplot(3, 1, 2);
plot(f_axis, abs(M_t_nat));
title('Transformada de Fourier muestreo natural');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
subplot(3, 1, 3);
plot(f_axis, abs(M_t_inst));
title('Transformada de Fourier muestreo instantaneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');

bit_depth=8;
pcm_levels=2^bit_depth;

pcm_signal_inst=round((m_t_inst+1)*(pcm_levels-1)/2);
grid on;
