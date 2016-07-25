clear all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Schranke für Ungleichbedingung
epsilon = 1e-4;

% Anzahl von Play-Operatoren
% N = 51;
N = 51;

% Eingangssignal [V]
load('x.mat')
x = x-3;
x = x';
% Ausgangssignal [V]
load('y.mat')
y = y - 3.2;

% Initialisierung der Startwerte der Operatorausgänge
% y_H0 = mean(y)*ones(N, 1);
y_H0 = 0*ones(N, 1);

% Schwellenwerte der Playoperatoren werden hier inisialisiert
% Nach Klaus Kuhnen Dissertation:
r_H = zeros(N,1);
r_H(1) = 0;
parfor i=2:N
    r_H(i) = i/N * max(abs(x));
end

%STARTWERTE GEWICHTE NACH KUHNEN:
w_H = zeros(N,1);
w_H(1) = 1;

% Matrizen der superpositionierten Operatoren
H_rH = zeros(N,length(x));
H_rH_inv = zeros(N,length(x));

%%
% %HYSTERESEBERECHNUNG
parfor i = 1:N
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% FEHLERMODELL          %
%%%%%%%%%%%%%%%%%%%%%%%%%

%Errorfunktion



w_MK_tilde = w_H;
M_Kr = H_rH;

% E_MK = x + (w_MK_tilde' * M_Kr);

H = zeros(length(w_MK_tilde), length(w_MK_tilde));
f = zeros(length(w_MK_tilde),1);
c = zeros(1,1);
parfor i=1:length(x)
H = H + M_Kr(:,i)*M_Kr(:,i)';
f = f + (-y(i))*M_Kr(:,i);
c = c + y(i) * y(i);
end

%%
% Constrains       
U_H = eye(N,N);

u_H = zeros(N,1);
u_H(1) = epsilon;



U_MK = U_H;

u_MK = u_H;

%Parameter für Inequality constraints A*x <= b
%Vergleichen mit U_MK * w - u_MK >= 0
A = -U_MK;
b = -u_MK;

%%
%Minimierungsverfahren

%Lower bound

options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

[w_MK_tilde_ident,fval,exitflag,output,lambda] = quadprog(H,f,A,b,[],[],[],[],[],options);

fprintf('Quadratischer Fehler:                   %d\n', fval+c/2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Hysterese mit identifizierten Parametern                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%HYSTERESEBERECHNUNG

parfor i = 1:N
    % Trajektorie des Playoperators
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end

H_delta = w_MK_tilde_ident' * H_rH;

figure(1)
plot(x, H_delta, 'r')
hold on
plot(x, y, 'b')
legend('Identifizierte Daten', 'Reelle Daten')

w_H_inv_ident = fun_w_H_inv(w_MK_tilde_ident);
r_H_inv = fun_r_H_inv(r_H, w_MK_tilde_ident);
y_H0_inv = fun_y_H0_inv(w_H_inv_ident, y_H0);
parfor i = 1:N
    % Trajektorie des Playoperators
    H_rH_inv(i,:) = fun_trajectory_PlayOperator(H_delta, r_H_inv(i), y_H0_inv(i));
end

H_delta_inv = w_H_inv_ident'*H_rH_inv;
figure(2)
plot(x,H_delta_inv);