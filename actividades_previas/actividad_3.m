clc;
close all;
clear all;
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
m_t=sin(2*pi*f_c*t);

r=floor(t_s/tm);
s=floor(tau/tm);
disp(r)
m_t_inst=zeros(1,length(t));
for i=1:length(m_t)
if mod(i,r)==0
m_t_inst(i:i+s)=m_t(i);
end
end
m_t_inst=m_t_inst(1:length(t));

plot(t,m_t);
hold on;
%plot(t,r)
%plot(t,m_t_nat,'-r');
plot(t,m_t_inst);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal PAM con Muestreo instantaneo');
grid on;
