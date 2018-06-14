% Identificação de Sistemas
% Questão 1 do projeto final de Identificação de sistemas
% Autores: Nicollas e Magdiel
close all
clear
clc

%% Definição dos dados
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

% % plotando a resposta real
% figure
% plot(u),title('Entrada PRBS'),
% grid on

%% METODO 1 - MQ

% Numero de parametros a estimar
na = 2;
nb = 1;

dim = na+nb+1;
m = max(na,nb+1);

% matriz de observa��o
phi = zeros(Ni,dim);
for t=m+1:Ni
    phi(t,:) = [-y(t-(1:na))' u(t-(0:nb))'];
end

% valor da matriz estimada dos parametros
theta = phi'*phi\phi'*y;
% fun��o estimada
yest = phi*theta;

% obtendo a fun��o de transferencia
den = [1 theta(1:na)'];
num = [theta(na+1:dim)' 0];
% fun��o discreta
Hs = tf(num,den,1);
% fun��o continua
Hsc = d2c(Hs);

% RESULTADOS
% Dados que foram Estimados
figure
plot(y)
hold on
plot(yest)
title('Dados de saida na Estima��o')
legend('Real','Estimado')
grid on

% Validando os dados
uv = idinput(Nv,'prbs',[0 0.25],[0,1]);
Tv = (1:Nv)';

yv = lsim(G,uv,Tv);
yest = lsim(Hsc,uv,Tv);

% Compara��o entre os dados de Valida��o
figure
plot(yv)
hold on
plot(yest)
title('Dados de saida na valida��o')
legend('Real','Estimado')
grid on

