%% Métodos Numéricos 2018 - Lista 1 - questões 1 a 4

% Professor: Cézar Santos
% Aluno: Gustavo Bulhões Carvalho da Paz Freire

clear all
close all
clc

%% Questão 1 - Tauchen (1986)

% z(t) = ro*z(t-1) + eps(t), onde eps ~ N(0,sigma^2)

% Parâmetros 

N = 9;               % enunciado pede 9 pontos
ro = 0.95;           % parâmetro, coeficiente do AR(1)
sigma = 0.007;       % parâmetro, desvio padrão do eps ~ N(0,sigma^2)
m = 3;               % scaling parameter (estou escolhendo 3)
z = zeros(N,1);      % pré-definindo vetor coluna Nx1 dos estados
P = zeros (N,N);    % pré-definindo matriz de transição NxN (linhas são estado atual,
                     % colunas são estado futuro.
                     
% Especificando o grid                     
                     
z(N) = m*(sigma/sqrt(1-ro^2));
z(1) = -z(N);
delta_z = (z(N)-z(1))/(N-1);
z = linspace (z(1), z(N), N)';

% Matriz de transição

for j = 1:N
    for k = 1:N
        if k == 1
            P(j,k) = normcdf(((z(1) - ro*z(j) + delta_z/2) / sigma));
        elseif k == N
            P(j,k) = 1 - normcdf(((z(N) -ro*z(j) - delta_z/2) / sigma));
        else
            P(j,k) = normcdf(((z(k) + delta_z/2 - ro*z(j) ) / sigma)) - normcdf(((z(k) - delta_z/2 - ro*z(j)) / sigma));
        end
    end
end


%% Questão 2 - Rouwenhorst

% Especificando o grid

sigma_theta2 = (sigma^2/(1 - ro^2));
theta(N) = sqrt(sigma_theta2)*sqrt(N-1);
theta(1) = -theta(N); 
theta = linspace(theta(1), theta(N), N)';
p = (1+ro)/2;

% Matriz de transição

Pr = [ p  (1-p);         % P2
        (1-p) p];


for t = 3:N              % PN recursivamente
      v = zeros(t-1,1);
      a = [Pr v;v' 0];
      b = [v Pr;0 v'];
      c = [v' 0;Pr v];
      d = [0 v';v Pr];
      Pr_prime = p*a + (1-p)*b + (1-p)*c + p*d;
      S = [ones(1,t);0.5*ones(t-2,t);ones(1,t)];
      Pr = S.*Pr_prime;
end


%% Questão 3 - Simulando para 10000 períodos e comparando

time = 10000;

% AR(1) contínuo 

y = NaN*ones(time, 1);
y(1) = 0;          % ponto inicial
e = ones(1,time);

for t = 1:time     % epslon (normal padrão)
    e(t) = randn;
end


for t = 1:(time-1)             % simulando
    y(t+1) = ro*y(t) + sigma*e(t);
    
end


% AR(1) discretizado (Tauchen)

tauch = ones(time,1);
tauch(1) = 1;
cum   = cumsum(P,2);

for t = 2:time                % simulando
    x = find(normcdf(e(t)) <= cum(tauch(t-1),:));
    tauch(t) = x(1);
end

tauchen = z(tauch);    % recuperando os valores no grid

% AR(1) discretizado (Rouwenhorst)

rouwen = ones(time,1);
rouwen(1) = 1;
cum_r   = cumsum(Pr,2);

for t = 2:time              % simulando
    x = find(normcdf(e(t)) <= cum_r(rouwen(t-1),:));
    rouwen(t) = x(1);
end

rouwenhorst = theta(rouwen);    % recuperando no grid

% Comparativo

tempo = linspace(1,time,time);   % todos os períodos

figure(1)
plot (tempo, y, tempo, tauchen, tempo, rouwenhorst);
xlabel('Tempo');
ylabel('Processo');
title('Comparativo das simulações do AR(1) contínuo e discretizado');
legend('Contínuo', 'Tauchen', 'Rouwenhorst');

tempo2 = linspace(1, 1000, 1000);    % 1000 primeiros períodos
y1000 = y(1:1000,1);
tauchen1000 = tauchen(1:1000,1);
rouwenhorst1000 = rouwenhorst(1:1000,1);

figure(2)
plot (tempo2, y1000, tempo2, tauchen1000, tempo2, rouwenhorst1000);
xlabel('Tempo');
ylabel('Processo');
title('Comparativo das 1000 primeiras simulações do AR(1) contínuo e discretizado');
legend('Contínuo', 'Tauchen', 'Rouwenhorst');

figure(3)
plot (tempo, y, tempo, tauchen);
xlabel('Tempo');
ylabel('Processo');
title('Comparativo das simulações do AR(1) contínuo e discretizado');
legend('Contínuo', 'Tauchen')

figure(4)
plot (tempo, y, tempo, rouwenhorst);
xlabel('Tempo');
ylabel('Processo');
title('Comparativo das simulações do AR(1) contínuo e discretizado');
legend('Contínuo', 'Rouwenhorst')


%% Questão 4 - Estimando AR(1) com base nos dados simulados

% Temos que ro_hat = cov(yt, yt-1)/var(yt) e sigma2_hat = (1-ro^2)var(yt)

ro_hat_tauchen = cov(tauchen(1:time-1,1), tauchen(2:time,1))/var(tauchen(1:time,1));
ro_hat_tauch = ro_hat_tauchen(1,2)

ro_hat_rouwenhorst = cov(rouwenhorst(1:time-1,1), rouwenhorst(2:time,1))/var(rouwenhorst(1:time,1));
ro_hat_rouwen = ro_hat_rouwenhorst(1,2)

sigma2_hat_tauchen = (1 - ro_hat_tauchen(1,2)^2)*var(tauchen(1:time,1));
sigma_hat_tauch = sqrt(sigma2_hat_tauchen)

sigma2_hat_rouwenhorst = (1 - ro_hat_rouwenhorst(1,2)^2)*var(rouwenhorst(1:time,1));
sigma_hat_rouwen = sqrt(sigma2_hat_rouwenhorst)

