N=3; %cantidad de bits configurables
Fs=1000;
t=0:1/Fs:1;
m_t=sin(2*pi*50*t);
L=2ˆN; %numero de niveles de cuantificacion

minst_t=m_t;
m_quant=round(minst_t*(L-1));
m_quant=m_quant/(L-1);

%calculo del error de cuantificacion
error_quant=minst_t-m_quant;

figure;
subplot(3,1,1);
plot(t, m_t);
title('señal original m(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3,1,2);
stem(t, minst_t, 'filled');
title('señal PAM Muestreada Instantaneamente minst(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(3,1,3);
stem(t, m_quant, 'filled');
title('señal PAM Cuantificada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
figure;
plot(t, error_quant);
title('error de cuantizacion');
xlabel('Tiempo (s)');
ylabel('Error');
