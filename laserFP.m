%%%% FBP MULTICAPA %%%%% - ANTONIO DELGADO BEJARANO

clear, close,clc;

% PROGRAMA DE SIMULACION PARA PERTURBACIONES PERIODICAS CON n COMPLEJO

% UNIDADES SIST. INTERNACIONAL

% CONSTANTES
c = 2.99793e8;

%% Ganancia optica:

g = 0;
%g = 2e-4;
%g = 3.5e-4;
%g = 4.8336e-4;
%g = 5000;
f = 193e12; % Frecuencia de emision deseada
nC = c*inv(4*pi*f)*g;

%% Definicion y muestreo de la Perturbacion que origina la FG

perturbacion = (3.5 + 1j*nC); % Cavidad fria - Resonador
     
Lper = 300e-6; % Longitud de la perturbacion

% Muestreo de la perturbacion

    % Periodo de muestreo
    z_muestreo = Lper/200;
    
    % Muestreo de la perturbacion
    [capas] = perturbacion.*ones(floor(Lper/z_muestreo),1);
    
% Creacion de las muestras que originan la cavidad FP
% (Espejos)
 n_z = [(1-1j*0); (1-1j*0); capas; (1-1j*0); (1-1j*0)]; % Vector columna

 % Representacion indice de refraccion
subplot(211), plot(z_muestreo.*(0:1:length(n_z)-1).',real(n_z)); ylabel('Parte real Indice Refracci?n')
subplot(212), plot(z_muestreo.*(0:1:length(n_z)-1).',imag(n_z)); xlabel('Eje z'), ylabel('Parte Imaginaria Indice Refracci?n')

%% Caracterizacion macroscopica dispositivo:
    % Funciones de transferencia y respuestas impulsivas en Rx y Tx
    
    % Componentes de frecuencia optica a las que se analiza el dispositivo
    N_frec = 2.^12;  % Num. de frecuencias.
    
    f_i = linspace(191e12, 195e12, N_frec); % Vector fila
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

    % Funcion de transferencia en reflexion
    r_0L = reshape(MT(2,1,:)./MT(1,1,:),1,N_frec); 
        
    % Funcion de transferencia en transmision
    t_0L = reshape(1./MT(1,1,:),1,N_frec); 

%% Reflexividad
figure(2)
plot(f_i, abs(r_0L).^2)
title('Reflectividad [u. n.]: Respuesta espectral en amplitud')
xlabel('Frecuencia Optica'); ylabel('Reflectividad [u. n.]');

%% Ref: Fase
phi_rH = unwrap(angle(r_0L));
figure(3)
plot(f_i, phi_rH)
title('Reflexi?n: Respuesta espectral de fase')
xlabel('Frecuencia Optica');
ylabel('Respuesta Espectral: Fase');

%% Ref: Tiempo grupo
t_rg = -diff(phi_rH)./(2*pi.*diff(f_i));

figure(4)
plot(f_i(1:end-1), t_rg)
title('Reflexi?n: Respuesta espectral de retardo de grupo')
xlabel('Frecuencia ?ptica');
ylabel('Retardo de Grupo');

%% Ref: Respuesta espectral
figure(5)
plot(f_i(1:end-2), diff(t_rg)./(2*pi.*diff(f_i(1:end-1))))
title('Reflexi?n: Respuesta espectral de dispersi?n')
xlabel('Frecuencia ?ptica');
ylabel('Dispersi?n');

%% Ref: Respuesta temporal
hr_0L = real(ifft([r_0L; conj(r_0L(1,(N_frec):-1:1))]));
hr_0L = hr_0L(1:N_frec);
figure(1), plot(t_i(1:end-1), abs(hr_0L)); xlabel('Tiempo (ns)'); 
ylabel('Optical magnitude [n. u.]');

%% Transmitividad
figure(7)
plot(f_i, abs(t_0L).^2)
title('Transmitividad [u. n.]: Respuesta espectral en amplitud')
xlabel('Frecuencia ?ptica'); ylabel('Transmitividad [u. n.]');

%% Trans: Fase
phi_tH = unwrap(angle(t_0L));
figure(8)
plot(f_i, phi_tH)
title('Transmisi?n: Respuesta espectral de fase')
xlabel('Frecuencia ?ptica');
ylabel('Respuesta Espectral: Fase');

%% Trans: Tiempo de grupo
t_tg = -diff(phi_tH)./(2*pi.*diff(f_i));

figure(9)
plot(f_i(1:end-1), t_tg)

title('Transmisi?n: Respuesta espectral de retardo de grupo')
xlabel('Frecuencia ?ptica');
ylabel('Retardo de Grupo');

%% Trans: Dispersion
figure(10)
plot(f_i(1:end-2), diff(t_tg)./(2*pi.*diff(f_i(1:end-1))))
title('Transmisi?n: Respuesta espectral de dispersi?n')
xlabel('Frecuencia ?ptica');
ylabel('Dispersi?n');

%% Tans: Respuesta temporal
hr_0L = real(ifft([t_0L; conj(t_0L(1,(N_frec):-1:1))]));
hr_0L = hr_0L(1:N_frec);
figure(1), plot(t_i(1:end-1), abs(hr_0L)); xlabel('Tiempo (ns)'); 
ylabel('Optical magnitude [n. u.]');
