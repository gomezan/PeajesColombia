% Peajes en Colombia
rng('default') % For reproducibility

%El codigo asume que las primeras posiciones de los vectores corresponden a
%los cajeros tipo C (solo tag), las del medio son tipo B y las primeras
%tipo A

%N, total de casetas | #B numero de casetas tipo B | #C, numero de casetas tipo C
peaje1=[15,3,2];
peaje2=[8,2,1];
peaje3=[4,1,0];


%[tEsperaProm, nCarros, tUserTag, tsA] = simularPico(peaje2, 0.1);
%tEsperaProm, 
%nCarros

%[tEsperaPromedio,nCarrosEnCola ]=simularCompleto(peaje2)

%Imprimir resultados finales 
%display(tEsperaPromedio)
%display(nCarrosEnCola)

%[finalt,finalc]=busqueda();

%graficar(finalt,"Tiempo de espera")
%graficar(finalc,"Número de carros")
%graficarSecciones(finalt)
%graficarSecciones(finalc)


%*******************************************
% Graficar secciones de los 52 fines de semana
function graficarSecciones(resultado)

% Crear tres arreglos entre 1 y 10
arreglo1 = 1:10;
arreglo2 = 1:10;
arreglo3 = 1:10;

% Obtener todas las combinaciones posibles
combinaciones = combvec(arreglo1, arreglo2, arreglo3);

figure;

% Colorbar
colormap('hot');
colorbarHandle = colorbar;

for i = 1:10
    % Seleccionar valores correspondientes a una sección
    indices = (combinaciones(3, :) == i);

    % Crear subgráfico
    subplot(2, 5, i);
    scatter(combinaciones(1, indices), combinaciones(2, indices), 50, resultado(indices), 'filled');

    xlabel('Casetas A');
    ylabel('Casetas B');
    title(['Sección casetas C ', num2str(i)]);
    
    % mismo colorbar
    caxis([min(resultado), max(resultado)]);
end

end


%*****************************************************************
%Función 1 fin de semana
function [tEsperaProm, nCarros, tUserTag, tsA] = simularPico(peaje, tagInicial)
%Tasas de servicio
t_s_tag = 0.5;
t_s_efe = 1.2;

N=peaje(1); % # de casetas
n_C=peaje(3);  % # de casetas tipo C
n_B=peaje(2);  % # de casetas tipo B

%T_c arreglo que representa el tiempo de espera acumulado por cada
%caseta.
%colas arreglo que representa el tamaño de colas de cada una de las
%casetas.
%t_last arreglo que almacena el tiempo de llegada del último auto de cada caseta. 
%T arreglo que almacenan los tiempos de espera puntuales por cada auto que llega a cada caseta.
%tsTag es un arreglo que almacena los tiempos de espera de usuarios tag.
[T_c, colas, t_last] = deal(zeros(1,N));
T = cell(1,N);
C = cell(1,N);
tsTag = [];


n_a = [random('Poisson',0.5*30*N) random('Poisson',1.2*60*N) random('Poisson',0.5*30*N)]; % number of arrivals
t_a = [30*rand(1,n_a(1)) 60*rand(1,n_a(2))+30 30*rand(1,n_a(3))+90]; % times of arrivals

noTag=1-tagInicial;
tag_flag = rand(1,size(t_a,2)) > noTag; % does the arrival have a tag?

%Se crea el vector de eventos t_ contiene los tiempos, tipo (llegada o salida) y si es de tag o no. 
t_ = sort(t_a)'; % sorted arrival times
t_(:,2) = zeros(size(t_,1),1); % arrival flag
t_(:,3)=tag_flag;
t = t_(1);

fin = 120; % tiempo de simulación

tam=size(t_,1)*2;

for i =  1:tam

    %Revisar si entra en la ventana de análisis, es decir:
    %Se encuentra dentro del vector de eventos
    long=size(t_,1);
    if i>long
        break
    end

    %t es el tiempo actual, type el tipo de evento y tags indica si el usuario tiene tag 
    t = t_(i,1);
    type = t_(i,2); 
    tags= t_(i,3);

    %Se encuentre dentro de los 120 minutos de tiempo?
    if t>fin
        break
    end

    if ~type % arrival event
        %Selecciona la caseta con la cola más corta de acuerdo a su tipo de pago
        if tags
            [~, m] = min(colas); % Shortest queue index
            t_s = random('Exponential',t_s_tag);
        else
            [~, m] = min(colas(n_C + 1:N)); % Shortest cash payment queue index
            m = m + n_C;
            t_s = random('Exponential',t_s_efe);
        end
        colas(m) = colas(m) + 1;
        if t_last(m) == 0
            t_last(m) = t;
        end
        
        %El tiempo acumulado por caseta tiene en cuenta el tiempo
        %transcurrido hasta el tiempo actual
        T_c(m) = T_c(m) - (t - t_last(m)) + t_s; % cumulative time for each queue
        t_(end+1,:) = [t + T_c(m), m, tags]; %#ok<SAGROW> % new departure time %si es una salida el tipo es el indice de la caseta que lo atiende. 
        t_ = sortrows(t_,1);

        %Busca almacenar los tiempos de espera promedio de los usuarios de tag  
        if tags
         tsTag = [tsTag, T_c(m)];
        end
        
        %concatena el contenido de T_c(m) al contenido previo de T{m}
        T{m} = [T{m} T_c(m)]; % waiting time
        t_last(m) = t;
        % Se repite el procedimiento pero esta vez con el tamaño de las colas 
        C{m} = [C{m} colas(m)];
        
    else % departure event
        colas(type) = colas(type) - 1;
        if colas(type) == 0
            t_last(type) = 0;
            T_c(type) = 0;
        end
    end    
end

tEsperaProm=cellfun(@mean,T);
nCarros=cellfun(@mean,C);
tUserTag=tsTag;
tsAtent=[];

%Busca almacenar los tiempos de espera de las casetas tipo A
soloA=T(n_C + n_B+ 1:N);
for i = 1:length(soloA)
    tsAtent = [tsAtent, soloA{i}];
end

tsA=tsAtent;

end


%**************************************************************
%Simluar 52 fines de semana
function [tEsperaProm, nCarros] =simularCompleto(peaje)

tamTag=0.1;

for i=1:54

[tEsperaPromedio, nCarrosEnCola, tEsperaTag, tEsperaA] = simularPico(peaje,tamTag);

%Calcular tiempos de espera promedio
promA=mean(tEsperaA);
promTag=mean(tEsperaTag);

%Si se cumple la primera regla
if promA> 10*promTag
 tamTag=tamTag+0.05;
else
    if promA> 3*promTag
     tamTag=tamTag+0.02;
   end
end

end 

tEsperaProm=tEsperaPromedio;
nCarros=nCarrosEnCola;

end

%*********************************************************************

function graficar(resultado, titulo)

 % Crea tres arreglos entre 1 y 10
    casetaA = 1:10;
    casetaB = 1:10;
    casetaC = 1:10;

% Obtiene todas las combinaciones posibles
    combinaciones = combvec(casetaA, casetaB, casetaC);

% Graficar en 3D
figure;
scatter3(combinaciones(1, :), combinaciones(2, :), combinaciones(3, :), 50, resultado, 'filled');
xlabel('Casetas A');
ylabel('Casetas B');
zlabel('Casetas C');
title(titulo);

% Establecer un mapa de colores
%colormap([1 0 0; 0 0 0]);
colormap('hot');

% Añadir barra de colores
h = colorbar;
ylabel(h, 'Resultado'); 

end

%Busqueda en malla
%*********************************************************************************
function [tproms,cproms]= busqueda()
    % Crea tres arreglos entre 1 y 10
    casetaA = 1:10;
    casetaB = 1:10;
    casetaC = 1:10;
   
    %Vectores para almacenar resultados de simulación
    promts=[];
    promCr=[];

    % Obtiene todas las combinaciones posibles
    combinaciones = combvec(casetaA, casetaB, casetaC);

    for i = 1:size(combinaciones, 2)
        i
        % Extrae los tres valores de cada combinación
        A = combinaciones(1, i);
        B = combinaciones(2, i);
        C = combinaciones(3, i);
        N=A+B+C;
        [tEsperaPromedio,nCarrosEnCola ]=simularCompleto([N,B,C])
        promts=[promts mean(tEsperaPromedio)];
        promCr=[promCr mean(nCarrosEnCola)];
    end

    tproms=promts;
    cproms=promCr;
end

