clear;
clc;

A = eye(3);
%Consideriamo solo u1 u2 u5 u6 come variabili di controllo, dato che u3 e
%u4 dipendono da u2
B = [ 1 -3  0  0
      0 -7  0  0
      0 10  1 -1];
%Sistema completamente osservabile => C = I
C = eye(3);

%Se Q = 0.1 * C' * C si ha una convergenza molto più lenta
Q = C' * C;
QF = Q;
%ro = 0.3 => convergenza all'ultimo passo
%ro = 0.0001 => convergenza in un passo
ro = 0.1;
R = ro * [1    0    0    0
          0  101    0  -10
          0    0    1    0
          0  -10    0    1];

%Orizzonte da raggiungere
T = 7;
%Variabili di stato ed input
x = zeros(3, T + 1);
u = zeros(4, T);
%Costi
Jstate = zeros(1,T + 1);
Jcontrol = zeros(1, T);

%Calcolo matrice P Riccati
%Le dimesnsioni di P sono la dimensione dello stato e l'orizzonte da
%raggiungere +1
P = zeros(3, 3, T + 1);
P(:, :, T + 1) = QF;
for t = T : -1 : 1
    P(:, :, t) = Q + A' * P(:, :, t + 1) * A - A' * P(:, :, t + 1) * B * ...
                 inv(R + B' * P(:, :, t + 1) * B) * B' * P(:, :, t + 1) * A;
end

%Calcolo di K, necessario per il calcolo del controllo
K = zeros(4, 3, T);
for t = 1 : T
    K(:, :, t) = inv(R + B' * P(:, :, t + 1) * B) * B' * P(:, :, t + 1) * A;
end
 
%Calcolo dello stato e del controllo
x(:,1)=[4 7 2]';
for t = 1 : T
    u(:, t) = - K(:, :, t) * x(:, t);
    x(:, t + 1) = A * x(:, t) + B * u(:, t);
end

[K_inf, P_inf, e] = dlqr(A, B, Q, R, 0); 
 
%Grafici sugli andamenti degli stati e dei controlli nell'orizzonte definito
subplot(2, 1, 1);
plot(1 : (T + 1), x);
legend('state x');
subplot(2, 1, 2);
plot(1 : (T), u);
legend('control u');
 
for t = 1 : (T + 1)
    Jstate(1, t) = x(:, t)' * Q * x(:, t);
end
 
for t = 1 : T
    Jcontrol(1, t) = u(:, t)' * R * u(:, t);
end
 
JstateF = sum(Jstate);
JcontrolF = sum(Jcontrol);
J = JstateF + JcontrolF;