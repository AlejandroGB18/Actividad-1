% Limpieza de pantalla
clear all
close all
clc

% Declaración de variables simbólicas
syms th1(t) th2(t) th3(t) l1 l2 l3 t

% Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP = [0, 0, 0];

% Creamos el vector de coordenadas articulares
Q = [th1, th2, th3];
disp('Cordenadas articulares')
pretty(Q);

% Creamo el vector de velocidades articulares
Qp = diff(Q, t); % Utilizo diff para derivadas cuya variable de referencia no depende de otra: ejemplo el tiempo
disp('Velocidades articulares')
pretty(Qp);


%% ETAPA 1. Grados de libertad
% Número de grados de libertad del robt
GDL = size(RP,2) % Siempre se coloca 2, ya que indica la dimensión de las columnas
GDL_str = num2str(GDL); % Convertimos el valor numérico a una cadena de carácteres tipo string


%% Articulación 1
% Posición de la junta 1 respecto a 0
P(:,:,1) = [l1*cos(th1);
            l1*sin(th1);
                      0]; % Vector de posición indexado por página

% Matriz de rotación de la articulación 1 respecto a 0
R(:,:,1) = [cos(th1) -sin(th1) 0; % Análisis de robot péndulo
            sin(th1)  cos(th1) 0;
            0         0        1];


%% Articulación 2
% Posición de la junta 2 respecto a 1
P(:,:,2) = [l2*cos(th2); %% Solo th2 por que las transformadas de derivación homogéneas locales ya las consideran, así que ya no es necesario colocar la suma de th1
            l2*sin(th2);
                     0]; % Vector de posición indexado por página

% Matriz de rotación de la articulación 2 respecto a 1
R(:,:,2) = [cos(th2) -sin(th2) 0; % Análisis de robot péndulo
            sin(th2)  cos(th2) 0;
            0         0        1];


%% Articulación 3
% Posición de la junta 3 respecto 2
P(:,:,3) = [l3*cos(th3); 
            l3*sin(th3);
                     0]; % Vector de posición indexado por página

% Matriz de rotación de la articulación 3 respecto a 2
R(:,:,3) = [cos(th3) -sin(th3) 0; % Análisis de robot péndulo
            sin(th3)  cos(th3) 0;
            0         0        1];


%% ETAPA II. Transformaciones Homogéneas

% Creamos un vector de ceros
Vector_zeros = zeros(1, 3);

% Inicializamos las matrices de transformación homogeneas locales (la de
% cada junta)
A(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_zeros 1]);

% Inicializamos las matrices de transformación homogeneas globales (con respecto a la referencia)
T(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_zeros 1]);

% Inicializamos los vectores de posición vistos desde el marco de
% referencia inercial
PO(:,:,GDL) = P(:,:,GDL);

% referencia inercial
RO(:,:,GDL) = R(:,:,GDL);

for i = 1:GDL
    i_str = num2str(i);

    % Locales
    disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i) = simplify([R(:,:,i) P(:,:,i); Vector_zeros 1]);
    pretty (A(:,:,i));

    % Globales
    try
        T(:,:,i) = T(:,:,i-1)*A(:,:,i); % Obtienes la global de la parte del sistema a partir de la global del sistema anterior
    catch
        T(:,:,i) = A(:,:,i); % Caso específico cuando i = 1 nos marcaría error en try
    end
    disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i) = simplify(T(:,:,i));
    pretty(T(:,:,i));

 % Obtenemos la matriz de rotación "RO" y el vector de traslación "PO" de
 % la matriz de transformación Homogénea global T(:,:,GDL)
    
    RO(:,:,i) = T(1:3, 1:3, i);
    PO(:,:,i) = T(1:3, 4,   i);
    disp('Matriz de rotación')
    pretty(RO(:,:,i));
    disp('Vector de traslación')
    pretty(PO(:,:,i));
end


%% ETAPA III. Jacobianos

% Derivadas parciales de x respecto a th1, th2 y th3
Jv11 = functionalDerivative(PO(1,1,GDL), th1);
Jv12 = functionalDerivative(PO(1,1,GDL), th2);
Jv13 = functionalDerivative(PO(1,1,GDL), th3);

% Derivadas parciales de y respecto a th1, th2 y th3
Jv21 = functionalDerivative(PO(2,1,GDL), th1);
Jv22 = functionalDerivative(PO(2,1,GDL), th2);
Jv23 = functionalDerivative(PO(2,1,GDL), th3);

% Derivadas parciales de z respecto a th1, th2 y th3
Jv31 = functionalDerivative(PO(3,1,GDL), th1);
Jv32 = functionalDerivative(PO(3,1,GDL), th2);
Jv33 = functionalDerivative(PO(3,1,GDL), th3);

% Creamos la matriz del Jacobiano lineal
jv_d = simplify([Jv11, Jv12, Jv13; ...
                 Jv21, Jv22, Jv23; ...
                 Jv31, Jv32, Jv33]);
pretty(jv_d);

% Calculamos el jacobiano lineal y angular de forma analítica
% Inicializamos jacobianos analíticos (lineal y angular)
Jv_a(:,GDL) = PO(:,:,GDL); % Velocidad lineal
Jw_a(:,GDL) = PO(:,:,GDL); % Velocidad angular


for k = 1:GDL
    if ((RP(k)==0)|(RP(k)==1)) % Casos: articulación rotacional y prismática

        % Para las articulaciones rotacionales
        try
            Jv_a(:,k) = cross(RO(:,3,k-1), PO(:,:,GDL) - PO(:,:,k-1));
            Jw_a(:,k) = RO(:,3,k-1);
        catch % Casos específicos
            Jv_a(:,k) = cross([0,0,1], PO(:,:,GDL));% Matriz de rotación de 0 con respecto a 0 es la matriz de identidad, la posición previa
            Jw_a(:,k) = [0,0,1]; % Si no hay matriz de rotación previa se obtiene la matriz identidad
        end
    else
        % Para las articulaciones prismáticas
        try
            Jv_a(:,k) = RO(:,3,k-1);
        catch % Casos específicos
            Jv_a(:,K) = [0,0,1]; % Si no hay matriz de rotación previa se obtiene la matriz identidad
        end
            Jw_a(:,k) = [0,0,0];
    end
end

Jv_a = simplify(Jv_a);
Jw_a = simplify(Jw_a);
disp('Jacobiano lineal obtenido de forma analítica');
pretty(Jv_a);
disp('Jacobiano ángular obtenido de forma analítica')
pretty(Jw_a)


%% ETAPA IV. Velocidades

disp('Velocidad lineal obtenida mediante el Jacobiano Lineal')
V = simplify(Jv_a*Qp');
pretty(V);

disp('Velocidad angular obtenida mediante el Jacobiano Angular')
W = simplify(Jw_a*Qp');
pretty(W);
