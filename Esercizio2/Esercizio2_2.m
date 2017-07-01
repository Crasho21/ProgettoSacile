clear;
clc;

A = eye(3);
%Consideriamo solo u1 u2 u5 u6 come variabili di controllo, dato che u3 e
%u4 dipendono da u2
B = [ 1 -3  0  0
      0 -8  0  0
      0 10  1 -1];
C = [1  0  0
     0  0 -1];
 