% Peajes en Colombia
rng('default') % For reproducibility
total = 8; % # de casetas
n_solo_c = 1; % # de casetas tipo C;
[tEsperaPromedio, nCarrosEnCola] = simularPeaje(total, n_solo_c)

display(tEsperaPromedio)
display(nCarrosEnCola)

function [tEsperaProm, nCarros] = simularPeaje(N, n_C)
%Tasas de servicio
t_s_tag = 0.5;
%t_s_tag = random('Exponential',1000,1,s_tag);
t_s_efe = 1.2;
%t_s_efe = random('Exponential',1000,1,s_efe);
%T_c arreglo que representa el tiempo de espera acumulado por cada
%caseta.
%colas arreglo que representa el tamaño de colas de cada una de las
%casetas.
%t_last arreglo que almacena el tiempo de llegada del último auto de cada caseta. 
%T arreglo que almacenan los tiempos de espera puntuales por cada auto que llega a cada caseta. 
[T_c, colas, t_last] = deal(zeros(1,N));
T = cell(1,N);

%tasas_pico = [0.5 1 0.3];
n_a = [random('Poisson',0.5*30*N) random('Poisson',1.2*60*N) random('Poisson',0.5*30*N)]; % number of arrivals
t_a = [30*rand(1,n_a(1)) 60*rand(1,n_a(2))+30 30*rand(1,n_a(3))+90]; % times of arrivals
tag_flag = rand(1,size(t_a,2)) > 0.9; % does the arrival have a tag?

%Se crea el vector t_ de eventos contiene tiempos y tipo. 
t_ = sort(t_a)'; % sorted arrival times
t_(:,2) = zeros(size(t_,1),1); % arrival flag
t = t_(1);
fin = 120; % tiempo de simulación
c = 0; % contador
for i =  1:size(t_,1)
    %t es el tiempo actual y type el tipo de evento
    t = t_(i,1);
    type = t_(i,2);    
    
    if ~type % arrival event
        %Selecciona la caseta con la cola más corta de acuerdo a su tipo de pago
        if tag_flag(i)
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
        t_(end+1,:) = [t + T_c(m), m]; %#ok<SAGROW> % new departure time %si es una salida el tipo es el indice de la caseta que lo atiende. 
        t_ = sortrows(t_,1);
        
        %concatena el contenido de T_c(m) al contenido previo de T{m}
        T{m} = [T{m} T_c(m)]; % waiting time
        t_last(m) = t;
        
    else % departure event
        colas(type) = colas(type) - 1;
        if colas(type) == 0
            t_last(type) = 0;
            T_c(type) = 0;
        end
    end    
end

tEsperaProm=cellfun(@mean,T)
nCarros=colas

end
