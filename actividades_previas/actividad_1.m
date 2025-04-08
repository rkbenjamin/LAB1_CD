clc;
clear;
close all;

A=1;             
fc=1000;        
Ts=1/100000;    
t=0:Ts:5/fc;    

m_t=A*sin(2*pi*fc*t);

figure;
plot(t, m_t, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Sinusoidal m(t)');
grid on;
