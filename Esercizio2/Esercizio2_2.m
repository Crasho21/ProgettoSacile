clear;
clc;

%intervallo di dicretizzazione delta_t;
delta_t=0.1;
%orizzonte da 0 a 59 secondi;
orizz=9;

A = eye(3);
%Consideriamo solo u1 u2 u5 u6 come variabili di controllo, dato che u3 e
%u4 dipendono da u2
B = [ 1 -3  0  0
      0 -7  0  0
      0 10  1 -1];
C = [1  0  0
     0  0 -1];
D = zeros(2,4); 
sys = ss(A,B,C,D,delta_t);

nx = size(A,2);    %Number of states
ny = size(C,1);    %Number of outputs

%Covarianze del rumore
Qn = [4 2 0; 2 1 0; 0 0 1];
Rn = [1.4 0; 0 3];

R = 1;
ro = 0.1;
R = ro * [1    0    0    0
          0  101    0  -10
          0    0    1    0
          0  -10    0    1];

Q=0.1*eye(nx);
F=Q;

QXU = blkdiag(Q,R);
QWV = blkdiag(Qn,Rn);

%------------------------------------------------------
%calcolo la matrice di riccati;
P(:,:,orizz-1)=F;
E=B*inv(R)*B';
V=Q;
I=eye(nx);
for i=orizz-2:-1:1
   P(:,:,i)=A'*P(:,:,i+1)*inv(I+E*P(:,:,i+1))*A+V;
end

%calcolo la retroazione ottima di riccati;
for t=1:(orizz-2)
    L(:,:,t)=R^-1 * B'* (A')^-1*(P(:,:,t)-Q);
end

%serve?
%verifico con P infinito;
[L_INF,P_INF,e] = dlqr(A,B,Q,R,0);

%------------------------------------------------------
% %kalman
% 
%genero il rumore
media = [0 ,0, 0];
rng default  % For reproducibility

xsi = mvnrnd(media,Qn,orizz-1)';
eta = mvnrnd(media(1:ny),Rn, orizz-1)';
%media_x0 e covarianza x0;
alfa=[3; 1; 0];
sigma_x0=[1 0 0;
          0 1.5 1;
          0 1 2];

k=zeros(nx, ny, t);
mu=zeros(3,orizz-1);
x=zeros(3,orizz-1);
sigma=zeros(nx, nx, orizz-1);
 
sigma(:,:,1)=inv(inv(sigma_x0)+C'*inv(Rn)*C);
k(:,:,1)=sigma(:,:,1)*C'*inv(Rn);

y=zeros(2,orizz-2);

u=zeros(4,orizz-1);

%osservo y0;
y(:,1)=[4 -2]; %Valori presi dalle equazioni dell'esercizio, calcolati su x(:,1)
x(:,1)=[4 7 2]';

% MU è la stima dello stato
% L* MU è il controllo OTTIMO

mu(:,1)=alfa+k(:,:,1)*(y(1,1)-C*alfa);
for t=1:(orizz-2)
    
    u(:,t) = L(:,:,t)*mu(:,t);
    
    %parte del sistema... la x non la vedo nel controllo!;
    x(:,t+1)=A*x(:,t)-B*L(:,:,t)*mu(:,t)+ xsi(:,t);
    y(:,t+1)=C*x(:,t+1)+eta(:,t+1);
    %parte del controllo stimando lo stato;
    sigma(:,:,t+1)=inv(inv(A*sigma(:,:,t)*A'+Qn)+C'*inv(Rn)*C);
    k(:,:,t+1)=sigma(:,:,t+1)*C'*inv(Rn);
    mu(:,t+1)=A*mu(:,t)-B*L(:,:,t)*mu(:,t)+ k(:,:,t+1)*(y(1,t+1)-C*(A*mu(:,t)-B*L(:,:,t)*mu(:,t)));
    
    % non serve usare il filtro di kalman
    mu(1,t+1) = y(1,t+1);
    mu(3,t+1) = -1*y(2,t+1);
    
end

%eqm tra x e mu;
% eqm:  errore quadratico medio
for t=1:(orizz-2)
    eqm(:,t)=(x(:,t)-mu(:,t))'*(x(:,t)-mu(:,t));
end

Jstate=zeros(1,orizz-2);
Jcontrol=zeros(1, orizz-2);
for i=1:(orizz-2)
Jstate(1,i)=x(:,i)'*Q*x(:,i);
end
 
for i=1:orizz-2
Jcontrol(1,i)=u(:,i)'*R*u(:,i);
end

Jtot = Jstate + Jcontrol;

%Costo totale per confronto con pid
JFinale = 0;
for i=1:orizz-2
JFinale = JFinale + Jtot(i);
end

JFinale

subplot(7, 1, 1);
plot(1:(orizz-1),x(1,:) , 1:(orizz-1),mu(1,:));
legend('state x1');
subplot(7, 1, 2);
plot(1:(orizz-1),x(2,:) , 1:(orizz-1),mu(2,:));
legend('state x2');
subplot(7, 1, 3);
plot(1:(orizz-1),x(3,:) , 1:(orizz-1),mu(3,:));
legend('state x3');
subplot(7, 1, 4);
plot(1:(orizz-1),u);
legend('control u');
subplot(7, 1, 6);
plot(1:(orizz-2),eqm);
legend('estimation error');
subplot(7, 1, 7);
plot(1:(orizz-2),Jtot);
legend('cost');