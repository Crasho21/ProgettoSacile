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
      0 -8  0  0
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

QWV = blkdiag(Qn,Rn);

%------------------------------------------------------
% %kalman
% 
%genero il rumore
media = [0 ,0, 0];
rng default  % For reproducibility

% R = mvnrnd(MU,SIGMA) returns an n-by-d matrix R of random vectors 
% chosen from the multivariate normal distribution with mean MU, and 
% covariance SIGMA.
xsi = mvnrnd(media,Qn,orizz-1)';
eta = mvnrnd(media(1:2),Rn, orizz-1)';
%media_x0 e covarianza x0;
alfa=[3; 1; 0];
sigma_x0=[1 0 0;
          0 1.5 1;
          0 1 2];

sigma=zeros(nx, nx, orizz-1);
sigma(:,:,1)=inv(inv(sigma_x0)+C'*inv(Rn)*C);
k(:,:,1)=sigma(:,:,1)*C'*inv(Rn);
%osservo y0;
y(:,1)=[4 -2]; %Valori presi dalle equazioni dell'esercizio, calcolati su x(:,1)
x(:,1)=[4 7 2]';

%  MU è la stima dello stato 



%------------------------------------------
% CONTROLLORE PID
%   Abbiamo bisogno di 4 pid, uno per ingresso u


mu(:,1)=alfa+k(:,:,1)*(y(1,1)-C*alfa);

kp_1=[-3 0 0];
ki_1=[0.8 0.2 0.4];
kd_1=[0 0 0.9091];

kp_2=[0 0 0];
ki_2=[0 0 0];
kd_2=[0 0 0];

kp_3=[0 0 -1.10001];
ki_3=[0 0 0];
kd_3=[0 0 0];

kp_4=[0 0 0];
ki_4=[0 0 0];
kd_4=[-2 3 0];


last_errore= [0 0 0]';
integrale= [0 0 0]';


for t=1:(orizz-2)
    stima_errore_pid = mu(:,t)- [0 0 0]';
    integrale=integrale+stima_errore_pid;
    derivata=(stima_errore_pid -last_errore);
    u(1,t)=kp_1*stima_errore_pid + ki_1*integrale + kd_1*derivata;
    u(2,t)=kp_2*stima_errore_pid + ki_2*integrale + kd_2*derivata;
    u(3,t)=kp_3*stima_errore_pid + ki_3*integrale + kd_3*derivata;
    u(4,t)=kp_4*stima_errore_pid + ki_4*integrale + kd_4*derivata;
    
    %parte del sistema... la x non la vedo nel controllo!;
    x(:,t+1)=A*x(:,t)+B*u(:,t)+ xsi(:,t);
    y(:,t+1)=C*x(:,t+1)+eta(:,t+1);
    %parte del controllo stimando lo stato;
    sigma(:,:,t+1)=inv(inv(A*sigma(:,:,t)*A'+Qn)+C'*inv(Rn)*C);
    k(:,:,t+1)=sigma(:,:,t+1)*C'*inv(Rn);
    mu(:,t+1)=A*mu(:,t)+B*u(:,t)+ k(:,:,t+1)*(y(1,t+1)-C*(A*mu(:,t)+B*u(:,t)));
    
    % sto dicendo che x(3) e x(1) è misurata
    % perfettamente da y(2) perciò non serve usare il filtro di kalman
    mu(1,t+1) = y(1,t+1);
    mu(3,t+1) = -1*y(2,t+1);
    
    
    last_errore = stima_errore_pid;
end


Jstate=zeros(1,orizz-1);
Jcontrol=zeros(1, orizz-1);
ro=0.2;
R=ro*eye(4);
Q=eye(3);





for i=1:(orizz-2)
Jstate(1,i)=x(:,i)'*Q*x(:,i);
end
 
for i=1:orizz-2
Jcontrol(1,i)=u(:,i)'*R*u(:,i);
end

Jtot = Jstate + Jcontrol;

JFinale = 0;
for i=1:orizz-2
JFinale = JFinale + Jtot(i);
end

JFinale

% COMPARAZIONE DEI COSTI
% 
%   
%   : costo del controllo ottimo (vedi es 2.2)
%  : costo controllore PID
%
%







subplot(5, 1, 1);
plot(1:(orizz-1),x(1,:) , 1:(orizz-1),mu(1,:));
legend('state x1');
subplot(5, 1, 2);
plot(1:(orizz-1),x(2,:) , 1:(orizz-1),mu(2,:));
legend('state x2');
subplot(5, 1, 3);
plot(1:(orizz-1),x(3,:) , 1:(orizz-1),mu(3,:));
legend('state x3');
subplot(5, 1, 4);
plot(1:(orizz-1),u);
legend('control u');
subplot(5, 1, 5);
plot(1:(orizz-1),Jtot);
legend('cost');
