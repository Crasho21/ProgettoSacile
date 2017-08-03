clear;
clc;

%intervallo di dicretizzazione delta_t;
delta_t = 0.1;
%orizzonte da 0 a 59 secondi;
orizz = 9;

A = eye(3);
%Consideriamo solo u1 u2 u5 u6 come variabili di controllo, dato che u3 e
%u4 dipendono da u2
B = [1 -3  0  0
     0 -7  0  0
     0 10  1 -1];
C = [1  0  0
     0  0 -1];
D = zeros(2, 4); 
sys = ss(A, B, C, D, delta_t);

%Numero di stati
nx = size(A, 2);
%Numero di uscite
ny = size(C, 1);

%Covarianze del rumore
Qn = [4 2 0; 2 1 0; 0 0 1];
Rn = [1.4 0; 0 3];

R = 1;
ro = 0.1;
R = ro * [1    0    0    0
          0  101    0  -10
          0    0    1    0
          0  -10    0    1];

Q = 0.1 * eye(nx);
F = Q;

QXU = blkdiag(Q, R);
QWV = blkdiag(Qn, Rn);

%Calcolo la matrice di riccati;
P(:, :, orizz - 1) = F;
E = B * inv(R) * B';
V = Q;
I = eye(nx);
for i = orizz - 2 : -1 : 1
   P(:, :, i) = A' * P(:, :, i + 1) * inv(I + E * P(:, :, i + 1)) * A + V;
end

%Calcolo la retroazione ottima di riccati;
for t = 1 : (orizz - 2)
    L(:, :, t) = R^-1 * B' * (A')^-1 * (P(:, :, t) - Q);
end

%Verifico con P infinito;
[L_INF,P_INF,e] = dlqr(A,B,Q,R,0);

%Filtro di Kalman
%Generazione del rumore
media = [0 ,0, 0];
%Per rispoducibilita
rng default

xsi = mvnrnd(media, Qn, orizz - 1)';
eta = mvnrnd(media(1 : ny), Rn, orizz - 1)';
%Media x0 e covarianza x0;
alfa = [3; 1; 0];
sigma_x0 = [1  0    0;
            0  1.5  1;
            0  1    2];

k = zeros(nx, ny, t);
mu = zeros(3, orizz - 1);
x = zeros(3, orizz - 1);
sigma = zeros(nx, nx, orizz - 1);

sigma(:, :, 1) = inv(inv(sigma_x0) + C' * inv(Rn) * C);
k(:, :, 1) = sigma(:, :, 1) * C' * inv(Rn);

y = zeros(2, orizz - 2);

u = zeros(4, orizz - 1);

%Osservo y0
%Valori presi dalle equazioni dell'esercizio, calcolati su x(:,1)
y(:, 1) = [4 -2];
x(:, 1) = [4 7 2]';

%mu è la stima dello stato
%-L * mu è il controllo OTTIMO

mu(:, 1) = alfa + k(:, :, 1) * (y(1, 1) - C * alfa);
for t = 1 : (orizz - 2)
    u(:, t) = -L(:, :, t) * mu(:, t);
    
    for i = 1 : length(u(:, t))
        if(u(i, t) > 1)
            u(i, t) = 1;
        elseif(u(i, t) < -1)
            u(i, t) = -1;
        end
    end
    
    %Parte del sistema... la x non la vedo nel controllo!;
    x(:, t+1) = A * x(:, t) + B * u(:,t) + xsi(:, t);
    y(:, t+1) = C * x(:, t + 1) + eta(:, t + 1);
    %Parte del controllo stimando lo stato;
    sigma(:, :, t + 1) = inv(inv(A * sigma(:, :, t) * A' + Qn) + C' * inv(Rn) * C);
    k(:, :, t + 1) = sigma(:, :, t + 1) * C' * inv(Rn);
    mu(:, t + 1) = A * mu(:, t) + B * u(:, t) + k(:, :, t + 1) * (y(1, t + 1) ...
                   - C *(A * mu(:, t) + B * u(:, t)));
    
    %Per x1 e x3 non e' necessario usare Kalman
    mu(1, t + 1) = y(1, t + 1);
    mu(3, t + 1) = -1 * y(2, t + 1);
    
end

%Eqm(errore quadratico medio) tra x e mu
for t = 1 : (orizz - 2)
    eqm(:, t) = (x(:, t) - mu(:, t))' * (x(:, t) - mu(:, t));
end

Jstate = zeros(1, orizz - 2);
Jcontrol = zeros(1, orizz - 2);
for i = 1 : (orizz - 2)
    Jstate(1, i) = x(:, i)' * Q * x(:, i);
end
 
for i = 1 : orizz - 2
    Jcontrol(1, i) = u(:, i)' * R * u(:, i);
end

Jtot = Jstate + Jcontrol;

%Costo totale per confronto con pid
JFinale = 0;
for i = 1 : orizz - 2
    JFinale = JFinale + Jtot(i);
end

JFinale

subplot(7, 1, 1);
plot(1 : (orizz - 1), x(1, :), 1 : (orizz - 1), mu(1, :));
legend('State x1');
subplot(7, 1, 2);
plot(1 : (orizz - 1), x(2, :), 1 : (orizz - 1), mu(2, :));
legend('State x2');
subplot(7, 1, 3);
plot(1 : (orizz - 1), x(3, :), 1 : (orizz - 1), mu(3, :));
legend('State x3');
subplot(7, 1, 4);
plot(1 : (orizz - 1), u);
legend('Control u');
subplot(7, 1, 5);
plot(1 : (orizz - 2), eqm);
legend('Estimation error');
subplot(7, 1, 6);
plot(1 : (orizz - 2), Jtot);
legend('Cost');