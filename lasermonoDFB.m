%%%% FBP MULTICAPA %%%%% - ANTONIO DELGADO BEJARANO

clear, close,clc;

% PROGRAMA DE SIMULACION PARA PERTURBACIONES PERIODICAS CON n COMPLEJO

% UNIDADES SIST. INTERNACIONAL

% CONSTANTES
c = 2.99793e8;

%% Ganancia optica:

%g = 0;
%g = 15000;
%g = 20000;
%g = 37500;
g = 34500

% Longitud de onda pedida y frecuencia correspondiente
lambda0 = 1300e-9; f0 = c/lambda0;

n0 = 3.5; % Parte real
n1 = 2.069e-3; 
n2 = c*inv(4*pi*f0)*g; % Parte imaginaria

%% Definicion y muestreo de la Perturbacion que origina la FG
     
Lper = 125e-6; % Longitud de la perturbacion

% Periodo LAMBDA MAYUSCULA
LAMBDA = lambda0/2/n0;

% Muestreo de la perturbacion

    % Periodo de muestreo - con 11 o 13 basta
    z_muestreo = LAMBDA/25;
    
    vector_z = linspace(0,Lper,Lper/z_muestreo);
  
    % Introducimos salto de pi en el centro
    n_z = [n0-1j*0, ...
        n0 + n1*sin(2*pi*inv(LAMBDA).*vector_z(1:floor(length(vector_z)/2))) + 1j*n2, ...
        n0 + n1*sin(2*pi*inv(LAMBDA).*vector_z(floor(length(vector_z)/2)+1:end)+pi) + 1j*n2, n0-1j*0];

% Representacion indice de refraccion muestreado
subplot(211), plot(z_muestreo.*(0:1:length(n_z)-1).',real(n_z)); ylabel('Parte real Indice Refraccion')
subplot(212), plot(z_muestreo.*(0:1:length(n_z)-1).',imag(n_z)); xlabel('Eje z'), ylabel('Parte Imaginaria Indice Refraccion')



%% Caracterizacion macroscopica dispositivo:
    
    % Componentes de frecuencia optica a las que se analiza el dispositivo
    N_frec = 2.^10;  % Num. de frecuencias.
    
    f_i = linspace(f0-2e12, f0+2e12, N_frec); % Vector fila
    landa_i = 2.99793e8./f_i;
    
    % Estructura de datos en el dominio del tiempo:
    f_max = f_i(length(f_i))-f_i(1);
    %delta_f = f_i-f_i(1,1);
    f_muestreo = f_i(2)-f_i(1);
    
    t_muestreo = inv(2.*f_max);
    t_i = t_muestreo.*(0:1:N_frec);
    
    % Obtencion de la matriz de transferencia:
    % Matriz orden: 2*2*N_frec:

    MT = formTFF_MT(n_z,z_muestreo,f_i);

    % Funcion de transferencia en reflexion, redimensionando para no tener
    % 1x1xfrecuencias
    r_0L = reshape(MT(2,1,:)./MT(1,1,:),1,N_frec); 
        
    % Funcion de transferencia en transmisi?n
    t_0L = reshape(1./MT(1,1,:),1,N_frec); 

%% Reflectividad
figure(2)
plot(f_i, abs(r_0L).^2)
title('Reflectividad [u. n.]: Respuesta espectral en amplitud')
xlabel('Frecuencia ?ptica'); ylabel('Reflectividad [u. n.]');

%% Transmitividad
figure(7)
plot(f_i, abs(t_0L).^2)
title('Transmitividad [u. n.]: Respuesta espectral en amplitud')
xlabel('Frecuencia ?ptica'); ylabel('Transmitividad [u. n.]');

