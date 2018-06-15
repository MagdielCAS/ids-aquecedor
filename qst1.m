% Identificação de Sistemas
% Questão 1 do projeto final de Identificação de sistemas
% Autores: Nicollas e Magdiel
close all
clear
clc

%% Definição dos dados
plotar = 3; % 1 - MQ; 2 - MQR; 3 - Os dois; 4 - Nenhum

% Divisão das parcelas para identificação e validação
N = 200;    % numero de amostras.
pN = 0.3;   % parcela de N para identificação (0 < pN < 1).

Ni = N*pN;  % Parcela para identificação
Nv = N-Ni;  % Parcela para validação

% A função de transferencia
num = [1 5];
den = conv([1 1],[1 4]);
G = tf(num,den); % função de transferencia continua

Ts = 1;
Gd = c2d(G,Ts); % função de transferencia discreta

% Obtemos os valores de a1 e a2, b0 e b1 da função acima
[numData, denData] = tfdata(Gd,'v');
grau = length(denData);
% pega os parâmetros da função eliminando os zeros iniciais gerados pelos
% graus diferentes do numerador e denominador
a = denData(find(denData,1,'first')+1:grau); % a1 ... an
b = numData(find(numData,1,'first'):grau); % b0 ... bn

% resposta do sistema ao sinal PRBS e a um rudio branco
u = idinput(Ni,'prbs',[0 0.5],[0,1]);
Ti = (1:Ni)';

y = lsim(G,u,Ti);

% Numero de parametros a estimar
na = 2;
nb = 1;

dim = na+nb+1; % dimensão
m = max(na,nb+1); % n máximo
d=0; % atraso de entrada

% % plotando a resposta real
% figure
% plot(u),title('Entrada PRBS'),
% grid on

%% METODO 1 - MQ



% matriz de observação
phi = zeros(Ni,dim);
for t=m+1:Ni
    phi(t,:) = [-y(t-(1:na))' u(t-(0:nb))'];
end

% valor da matriz estimada dos parametros
theta = phi'*phi\phi'*y;
% função estimada
yest = phi*theta;

% obtendo a função de transferencia
den = [1 theta(1:na)'];
num = [theta(na+1:dim)' 0];
% função discreta
Hs = tf(num,den,1);
% função continua
Hsc = d2c(Hs);

% RESULTADOS
if plotar == 1 || plotar == 3
    % Dados que foram Estimados
    figure
    plot(y)
    hold on
    plot(yest)
    title('Dados de saida na Estimação MQ')
    legend('Real','Estimado')
    grid on

    % Validando os dados
    uv = idinput(Nv,'prbs',[0 0.25],[0,1]);
    Tv = (1:Nv)';

    yv = lsim(G,uv,Tv);
    yest = lsim(Hsc,uv,Tv);

    % Comparação entre os dados de Validação
    figure
    plot(yv)
    hold on
    plot(yest)
    title('Dados de saida na validação MQ')
    legend('Real','Estimado')
    grid on

end

%% Método 2 - MQR


% vetor de parametros
theta = [a(1:na) b(1:nb+1)]';

%valores iniciais
param = zeros(dim,Ni);
erro = zeros(Ni,1);
phi = zeros(dim,1);
yMqrEst = zeros(Ni,1);
p = 1000*eye(dim,dim);
e = zeros(Ni,1);

% MQR
for t = dim : Ni
   phi = [-yMqrEst(t-(1:na)); u(t-d-(0:nb))]; % alterando matriz de observação
   
   erro(t) = y(t)-phi'*theta; % calculando erro
   
   K = p*phi/(1+phi'*p*phi); % calculando ganho
   
   theta = theta+(K*erro(t)); % nova matriz de estimadores
   
   p = p - K*(phi'*p); % calculando matriz de covariância
   
   param(:,t) = theta; % atualizando vetor de parametros
   
   yMqrEst(t) = phi'*theta + e(t); %obtendo valores de saída
end

% obtendo a função de transferencia
denMqr = [1 theta(1:na)'];
numMqr = [theta(na+1:dim)' 0];
% função discreta
HsMqr = tf(numMqr,denMqr,1);
% função continua
HscMqr = d2c(HsMqr);

% RESULTADOS
if plotar == 2 || plotar == 3
    % Dados que foram Estimados
    figure
    plot(1:Ni,y)
    hold on
    plot(1:Ni,yMqrEst)
    title('Dados de saida na Estimação MQR')
    legend('Real','Estimado')
    grid on

    % Validando os dados
    uv = idinput(Nv,'prbs',[0 0.25],[0,1]);
    Tv = (1:Nv)';

    yv = lsim(G,uv,Tv);
    yMqrEst = lsim(HscMqr,uv,Tv);

    % Comparação entre os dados de Validação
    figure
    plot(1:Nv,yv)
    hold on
    plot(1:Nv,yMqrEst)
    title('Dados de saida na validação MQR')
    legend('Real','Estimado')
    grid on
end


