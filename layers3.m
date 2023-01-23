function MTt = layers3(n_z, z_muestreo, f_i)

c = 2.99793e8;

% Conjunto de matrices de transferencia 2x2, una para cada f_i
MTt = zeros(2,2,length(f_i));

% Valor de beta en funcion de f sin contar con el valor variable de n_z. 
% Posteriormente se anadira el valor de n_z para cada beta. La longitud 
% corresponde con la z_muestreo
beta_n = 1j*z_muestreo.*(2*pi*f_i)./c;

    % Recorremos todas las frecuencias, para cada frecuencia recorremos
    % todo el vector n_z
    for f = 1:length(f_i)
        
        % Primera capa de interfaz entre el primer valor del vector n_z
        % (que corresponde a ind. refraccion en aire) y el siguiente valor
        % de n_z (que ha de corresponder a ind.refraccion del material).
        z=1;
        
        Mdiel = [exp(beta_n(f)*n_z(z)) 0; 0 exp(-beta_n(f)*n_z(z))];
       
        % Matriz de interfaz entre capas de dielectrico:           
        Mint = inv(2*n_z(z))*[n_z(z)+real(n_z(z+1)) n_z(z)-real(n_z(z+1)); ...
                              n_z(z)-real(n_z(z+1)) n_z(z)+real(n_z(z+1))];
        
        % Hacemos que esta capa de interfaz sea la primera y sera a esta a
        % la que iremos multiplicando la posterior estructura de capas:
        MTt(:,:,f) = Mdiel*Mint;                            
        
        % Recorremos el resto del vector n_z. El ultimo valor de n_z ha de 
        % corresponder al de aire fuera del dielectrico. Puesto que en ese
        % caso ha de calcularse una matriz de interfaz entre el penulitimo
        % y el ultimo elemento y este calculo va implicito en la forma en 
        % que definimos las matrices de interfaces,se recorre el vector n_z 
        % hasta el penultimo elemento:
        for z = 2:(length(n_z)-1)
            
            % Matriz de una capa de dielectrico:
            Mdiel = [exp(beta_n(f)*n_z(z)) 0; 0 exp(-beta_n(f)*n_z(z))];
            
            % Matriz de interfaz entre capas de dielectrico:           
            Mint = inv(2*n_z(z))*[n_z(z)+real(n_z(z+1)) n_z(z)-real(n_z(z+1)); ...
                                  n_z(z)-real(n_z(z+1)) n_z(z)+real(n_z(z+1))];
            
            % Agrupamos las capas de dielectrico y las de interfaces:                     
            MTt(:,:,f) = MTt(:,:,f)*Mdiel*Mint;

        end
        
        z = length(n_z);
        Mdiel = [exp(beta_n(f)*n_z(z)) 0; 0 exp(-beta_n(f)*n_z(z))];
        MTt(:,:,f) = MTt(:,:,f)*Mdiel;

    end
    
end

