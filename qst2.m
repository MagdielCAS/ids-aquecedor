% Identificação de Sistemas
% Questão 2 do projeto final de Identificação de sistemas
% Autores: Nicollas e Magdiel
close all
clear
clc

%% Definição dos dados para identificação
dados = load('F0407.DAT');
u = dados(:,1);
y = dados(:,2);

%% Definição dos dados para Validação
dados_val = load('F0307.DAT');
u_val = dados_val(:,1);
y_val = dados_val(:,2);

% Amostras
N = length(y);
T = (1:N)';
e = 0.01*randn(N,1);

%% METODO 1 : método MQR

yest = zeros(N,1); % saida estimada
erro = zeros(1,N);  % erro de previsão

% numero de parametros em A e B 
na = 3;  % n e m dimensões para y-na e x-nb
nb = 1;

dim = na+nb+1;

% vetor de parametros
theta = zeros(dim,1);
% vetor de medidas
phi = zeros(dim,1);
% Matriz de covanriancia inicial
p = 1000*eye(dim,dim);

%% Aplica��o do MQR

% Definindo intervalos para as possi��es iniciais de phi()
m = min(na,nb);  % menor qdt de n's  entre A e B
lim = max(na,nb)+1;  % diferen�a de n's entre A e B

for t = 1:lim
    % valor de b0 sempre � u(t)
    phi(na+1) = u(t);
    % verificando os intervalos de phi
    if(t>1 && t<=m+1) % intervalo de a1,b1 at� am,bm
       for i = 1:t-1
           phi(i) = -y(t-i);         % a1 � am
           phi(na+1+i) = u(t-i);    % b1 � bm
       end
    elseif(t>m+1) % intervalo das sobras, onde na>nb ou nb>na 
       if (na>nb)
           for i= 1:t-1
               phi(i) = -y(t-i);     % a1 � a.na
           end
           for i= 1:nb
               phi(na+1+i) = u(t-i);    % b1 � b.nb
           end
       elseif (na<nb)
           for i= 1:t-1
               phi(na+1+i) = u(t-i);     % a1 � a.na
           end
           for i= 1:na
               phi(i) = -y(t-i);     % a1 � a.na
           end 
       end
    end
    
    % Calculando o erro de previsão
    erro(t)= y(t) - phi'*theta;
    % Calculando o ganho
    K = p*phi/(1+phi'*p*phi);
    % Calculando novo valor da matriz de estimadores
    theta = theta+(K*erro(t));
    % Calculando a matriz covariância
    p = p - K*phi'*p;
    % Saida Estimada
    yest(t) = phi'*theta + e(t);
end

% Novo calculo agora para o resto das amostras
for t = lim+1 : N
    % alterando os valores de phi
    phi = [-y(t-(1:na));u(t-(0:nb))];
    % Calculando o erro
    erro(t)= y(t) - phi'*theta;
    % Calculando o ganho
    K = p*phi/(1+phi'*p*phi);
    % Calculando novo valor da matriz de estimadores
    theta = theta+(K*erro(t));
    % Calculando a matriz covariância
    p = p - K*(phi'*p);
    % Saida estimada
    yest(t) = phi'*theta + e(t);
end
% Diferença entre as saidas
err = y-yest;

theta1 = theta;
den = [1 theta1(1:na)'];
num = [theta1(na+1:dim)' 0];
% função discreta
Hs1 = tf(num,den,1,'variable','z^-1');
% função continua
Hsc1 = d2c(Hs1);

%% Plotando os Resultados
% Comparação entre os dados Identificados
figure
plot(y)
title('Estima��o da Resposta de um Aquecedor')
hold on
plot(yest)
legend('Real','Estimado')
grid on

% Comparação entre os dados de Validação
figure
lsim(Hsc1,u_val,T)
hold on
plot(y_val)
title('Dados de saida para valida��o')
grid on

%% MÉTODO 2 - Usando MQ

% Numero de parametros a estimar
na = 2;
nb = 1;

dim = na+nb+1;
m = max(na,nb+1);

% matriz de observação
phi = zeros(N,dim);
for t=m+1:N
    phi(t,:) = [-y(t-(1:na))' u(t-(0:nb))'];
end

% valor da matriz estimada dos parametros
theta2 = phi'*phi\phi'*y;

% obtendo a função de transferencia
den = [1 theta2(1:na)'];
num = [theta2(na+1:dim)' 0];
% função discreta
Hs2 = tf(num,den,1,'variable','z^-1');
% função continua
Hsc2 = d2c(Hs2);

%% Validando os dados
% Comparação entre os dados de Validação
figure
lsim(Hsc2,u_val,T)
hold on
plot(y_val)
title('Dados de saida para valida��o')
grid on

%% METODO 3 - Usando a correlação





