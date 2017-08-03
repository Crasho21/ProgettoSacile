clear;
clc;

%La prima equazione dell'esercizio è errata, infatti rende impossibile
%controllare il livello del liquido, quindi come suggerito dal professore
%la sostituiamo con x2
A = [0  -2    0
     0   0    1
     0   0  -10];
B = [0   0      0
     0   0      0
     0   0   9000];
C = eye(3);
D = zeros(3);

T = 100;
Ts = 1;

sysc = ss(A, B, C, D);
sysd = c2d(sysc, Ts);

Ad = sysd.a;
Bd = sysd.b;
Cd = C;
Dd = D;

Q = [1000       0        0
     0     0.0001        0
     0          0   0.0001];
QF = Q;
%ro = 0.3 => convergenza all'ultimo passo
%ro = 0.0001 => convergenza in un passo
ro = 0.1;
R = ro * eye(3);

x = zeros(3, T + 1);
u = zeros(3, T);

%Definizione controllo ottimo ad orizzonte infinito a tempo discreto
[Kdinf, Pdinf, e] = dlqr(Ad, Bd, Q, R, 0); 

x0 = [2 -4 1]';

%Risolvo le equazioni alle differenze usando Riccati
P(:, :, T + 1) = Q;

for i = T : -1 : 1
    P(:, :, i) = Q + Ad' * P(:, :, i + 1) * Ad - Ad' * P(:, :, i + 1) * ...
     Bd * (inv(R + Bd' * P(:, :, i + 1) * Bd)) * Bd' * P(:, :, i + 1) * Ad;
end


for i = 1 : T
    Kd_fin(:, :, i) = inv(R + Bd' * P(:, :, i + 1) * Bd) * Bd' * P(:, :, i + 1) * Ad;
end


%Traccio il valore desiderato
Z = [10 0 0]';

%Calcolo di G;
G(:, T + 1) = C * Q * Z;
E = Bd * inv(R) * Bd';

W = C' * Q;

for k = T : -1 : 1
    G(:, k) = Ad' * (eye(3) - inv(inv(P(:, :, k + 1)) + E) * E) * ...
              G(:, k + 1) + W * Z;
end

%Calcolo di LG;
for k = 1 : T
    LG(:, :, k) = inv(R + Bd' * P(:, :, k + 1) * Bd) * Bd';
end

%Evoluzione del sistema in control loop
x(:, 1) = x0;
for k = 1 : T
    x(:, k + 1) = (Ad - Bd * Kd_fin(:, :, k)) * x(:, k) + Bd * ...
                  LG(:, :, k) * G(:, k + 1);
    u(:, k) = -Kd_fin(:, :, k) * x(:, k) + LG(:, :, k) * G(:, k + 1);
end

subplot(2, 1, 1);
plot(1 : (T + 1), x(1, :));
legend('State x');
subplot(2, 1, 2);
plot(1 : (T), u(3, :));
legend('Control u');