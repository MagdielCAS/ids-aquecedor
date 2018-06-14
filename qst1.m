% Identifica��o de Sistemas
% Quest�o 1 do projeto final de Identifica��o de sistemas
% Autores: Nicollas e Magdiel
clear
clc

%% Defini��o dos dados
% Divis�o das parcelas para identifica��o e valida��o
N = 200;    % numero de amostras.
pN = 0.3;   % parcela de N para identifica��o (0 < pN < 1).

Ni = N*pN;  % Parcela para identifica��o
Nv = N-Ni;  % Parcela para identifica��o

% A fun��o de transferencia
num = [1 5];
den = conv([1 1],[1 4]);
G = tf(num,den); % fun��o de transferencia continua

Ts = 1;
Gd = c2d(G,Ts); % fun��o de transferencia discreta

% Obtemos os valores de a1 e a2, b0 e b1 da fun��o acima
a = [-0.3862 0.0067];
b = [0.761 0.0147];

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

