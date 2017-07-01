clear;
clc;

A = eye(3);
%Consideriamo solo u1 u2 u5 u6 come variabili di controllo, dato che u3 e
%u4 dipendono da u2
B = [1 -3  0  0
     0 -8  0  0
     0 10  1 -1];
%Sistema completamente osservabile => C = I
C = eye(3);
D = zeros(3, 4);

%Se Q = 0.1 * C' * C si ha una convergenza molto più lenta
Q = [0.01 0 0
     0 0.01 0
     0 0 1000];
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
TS = 1;

x = zeros(3, T + 1);
u = zeros(4, T);

%Definizione sistema a tempo discreto
sysc=ss(A,B,C,D);
sysd=c2d(sysc,TS);
Ad=sysd.a;
Bd=sysd.b;
Cd=C;
Dd=D;

%Definizione controllo ottimo ad orizzonte infinito a tempo discreto
[Kdinf,Pdinf,e] = dlqr(Ad,Bd,Q,R,0); 

x0=[1 -2  3]';

%solve the problem at finite horizon;
%discrete time:
%solve the RICCATI equation at the difference
P(:,:,T+1)=Q;

for i=T:-1:1
P(:,:,i)=Q+Ad'*P(:,:,i+1)*Ad-Ad'*P(:,:,i+1)*Bd*...
     (inv(R+Bd'*P(:,:,i+1)*Bd))*Bd'*P(:,:,i+1)*Ad;
end

for i=1:T
    Kd_fin(:,:,i)=inv(R + Bd'*P(:,:,i+1)*Bd)*...
          Bd'*P(:,:,i+1)*Ad;
    %Kd_fin_vettore(i,:)=reshape(Kd_fin(:,:,i),4,1);
end

%Kd_fin_2simulink=[[1:T]' Kd_fin_vettore];

Z = [0 0 0]';
Z2 = [0 0 10]';

%compute G;
G=zeros(3,T+1);
G(:,T+1)=C*Q*Z;
E=B*inv(R)*B';
W=C'*Q;
for k=T:-1:1
    if(k == 2)
        G(:,k)=A'*(eye(3)-inv(inv(P(:,:,k+1))...
        +E)*E)*G(:,k+1)+W*Z2;
    else
        G(:,k)=A'*(eye(3)-inv(inv(P(:,:,k+1))...
        +E)*E)*G(:,k+1)+W*Z;
    end
end

%compute LG;
for k=1:T
    LG(:,:,k)=inv(R+B'*P(:,:,k+1)*B)...
    *B';
end

%evolution of the system in control loop
x(:,1)=x0;
for k=1:T
    x(:,k+1)=(A-B*Kd_fin(:,:,k))*x(:,k)+...
        B*LG(:,:,k)*G(:,k+1);
    u(:,k) =  -Kd_fin(:,:,k)*x(:,k) + LG(:,:,k)*G(:,k+1);
end
%plot(1:T+1,x);
subplot(2, 1, 1);
plot(1:(T+1),x);
legend('state x');
subplot(2, 1, 2);
plot(1:(T),u);
legend('control u');