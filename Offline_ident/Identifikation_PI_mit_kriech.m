clear all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Schranke für Ungleichbedingung
epsilon = 1e-4;

% Anzahl von Play-Operatoren
N = 51;

% Anzahl von elementaren Kriechoperatoren
M = 10;

% Abtastzeit (10 kHz)
T_s = 0.0001;

% Eingangssignal [V]
load('x.mat')
x = x-3;
x = x';
% Ausgangssignal [V]
load('y.mat')
y = y - 3.2;

% Initialisierung der Operator-/ Systemausgänge
y_kriech = zeros(length(x),1);

% Initialisierung der Startwerte der Operatorausgänge
% y_H0 = mean(y)*ones(N, 1);
y_H0 = 0*ones(N, 1);
y_K0 = 0*ones(N, 1);

% Schwellenwerte der Playoperatoren werden hier inisialisiert
% Nach Klaus Kuhnen Dissertation:
r_H = zeros(N,1);
r_H(1) = 0;
parfor i=2:N
    r_H(i) = i/N * max(abs(x));
end

r_K = zeros(N,1);
r_K(1) = 0;
for i=2:N
    r_K(i) = i/N * max(abs(x));
end

%STARTWERTE GEWICHTE NACH KUHNEN:
w_H = zeros(N,1);
w_H(1) = 1;

w_K = zeros(N, 1);

% Matrizen der superpositionierten Operatoren
H_rH = zeros(N,length(x));
H_rH_inv = zeros(N,length(x));
K_rK = zeros(N,length(x));

%%
% %HYSTERESEBERECHNUNG
parfor i = 1:N
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end
parfor i=1:M
    a_Kj(i) = 1/(10^(i-1)*T_s);
end

parfor j = 1:N
    for i = 1:M
        y_kriech = fun_trajectory_KriechOperator(x, y_K0(j), a_Kj(i), r_K(j), T_s);
        
        K_rK(j,:) = K_rK(j,:) + y_kriech;
         if (i == M)
             K_rK(j,:) = K_rK(j,:)/M;
         end
    end
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% FEHLERMODELL          %
%%%%%%%%%%%%%%%%%%%%%%%%%

%Errorfunktion



w_MK_tilde = [w_H; w_K];
M_Kr = [H_rH; K_rK];

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
U_K = eye(N,N);

u_H = zeros(N,1);
u_H(1) = epsilon;

u_K = zeros(N,1);

U_MK = zeros(N+N, N+N);
U_MK(1:N,1:N) = U_H;
U_MK(N+1:N+N, N+1:N+N) = U_K;


u_MK = [u_H; u_K];

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

w_H_ident = w_MK_tilde_ident(1:N);
w_K_ident = w_MK_tilde_ident(N+1:N+N);


%%
%HYSTERESEBERECHNUNG

parfor i = 1:N
    % Trajektorie des Playoperators
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end

H_delta = w_H_ident' * H_rH;

parfor i=1:M
    a_Kj(i) = 1/(10^(i-1)*T_s);
end
tic
parfor j = 1:N
    for i = 1:M
        y_kriech = fun_trajectory_KriechOperator(x, y_K0(j), a_Kj(i), r_K(j), T_s);
        
        K_rK(j,:) = K_rK(j,:) + y_kriech;
         if (i == M)
             K_rK(j,:) = K_rK(j,:)/M;
         end
    end
end
toc
K_delta = w_K_ident' * K_rK;


figure(1)
plot(x, (H_delta+K_delta), 'r')
hold on
plot(x, y, 'b')
legend('Identifizierte Daten', 'Reelle Daten')
