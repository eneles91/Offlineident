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

% Halbe Anzahl (+1) von Superpositionsoperatoren
L = 20;

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
y_H0 = 0*ones(N, 1);
y_K0 = 0*ones(N, 1);

% Schwellenwerte der Playoperatoren werden hier inisialisiert
% Nach Klaus Kuhnen Dissertation:
r_H = zeros(N,1);
r_H(1) = 0;
parfor i=2:N
    r_H(i) = i/N * max(abs(x));
end

% Schwellenwerte der Superpositionsoperatoren werden hier inisialisiert
% r_Sl_= linspace(min(x), 0-(abs(min(x)/L)), L);
% r_S0 = 0;
% r_Sl = linspace(0+(max(x)/L), max(x), L);
% r_S = [r_Sl_, r_S0, r_Sl]';
%NACH KUHNEN:
r_S_inv = zeros((2*L+1),1);
parfor i=1:L
    r_S_inv(i) = (i-(L+1)+0.5)/L*min(y);
end
parfor i=(L+2):(2*L+1)
    r_S_inv(i) = (i-(L+1)-0.5)/L*max(y);
end
% Schwellenwerte der Kriechoperatoren werden hier inisialisiert
r_K = zeros(N,1);
r_K(1) = 0;
for i=2:N
    r_K(i) = i/N * max(abs(x));
end

% Identifizierte Gewichte werden hier initialisiert
% w_H = linspace(1,3,N)';
% w_H = ones(N,1);
%STARTWERTE GEWICHTE NACH KUHNEN:
w_H = zeros(N,1);
w_H(1) = 1;

%STARTWERTE GEWICHTE NACH KUHNEN:
w_S_inv = zeros((2*L)+1,1);
w_S_inv(L+1) = 1;

% w_K = 1*linspace(2,3,N)'/60;
% w_K = ones(N, 1);
%STARTWERTE GEWICHTE NACH KUHNEN:
w_K = zeros(N, 1);

% Identifizierte Kriechoperatoren
a_Kj = zeros(length(M), 1);

% Matrizen der superpositionierten Operatoren
H_rH = zeros(N,length(x));
H_rH_inv = zeros(N,length(x));
S_rS = zeros((2*L+1),length(x));
S_rS_inv = zeros((2*L+1),length(x));
K_rK = zeros(N,length(x));
M_rM = zeros((2*L)+1,length(x));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung inverser Parameter                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_H_inv = fun_w_H_inv(w_H);

r_H_inv = fun_r_H_inv(r_H, w_H);

y_H0_inv = fun_y_H0_inv(w_H, y_H0);

% w_S_inv = fun_w_S_inv(w_S, L);
w_S = fun_w_S_inv(w_S_inv, L);

% r_S = fun_r_S_inv(r_S_inv, w_S_inv, L);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrizen für Ungleichungsbedingung                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_H = eye(N,N);
U_K = eye(N,N);
U_S = eye((2*L+1),(2*L+1));
for i=1:L+1
    parfor j=i:L+1
        U_S(i, j) = 1;
    end
end
for i=2*L+1:-1:L+1
    for j=i:-1:L+1
        U_S(i, j) = 1;
    end
end

u_H = zeros(N,1);
u_H(1) = epsilon;
u_K = zeros(N,1);
u_S = epsilon*ones((2*L+1),1);

U_MK = zeros(N+N+(2*L+1), N+N+(2*L+1));
U_MK(1:N, 1:N) = U_H;
U_MK(N+1:N+(2*L+1), N+1:N+(2*L+1)) = U_S;
U_MK(N+(2*L+1)+1:N+(2*L+1)+N, N+(2*L+1)+1:N+(2*L+1)+N) = U_K;

u_MK = [u_H; u_S; u_K];

w_MK = [w_H; w_S; w_K];

if ((U_MK*w_MK - u_MK) >= 0)
    disp('Ungleichungsbedingung erfüllt!')
else
    disp('Ungleichungsbedingung NICHT erfüllt!')
end
   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Ungleichsbedingung                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_H_inv = -1 * eye(N,N);
U_H_inv(1,1) = 1; 

U_MK_inv = zeros(N+N+(2*L+1), N+N+(2*L+1));
U_MK_inv(1:N, 1:N) = U_H_inv;
U_MK_inv(N+1:N+(2*L+1), N+1:N+(2*L+1)) = U_S;
U_MK_inv(N+(2*L+1)+1:N+(2*L+1)+N, N+(2*L+1)+1:N+(2*L+1)+N) = U_K;

w_MK_inv = [w_H_inv; w_S_inv; w_K];

if ((U_MK_inv*w_MK_inv - u_MK) >= 0)
    disp('Inverse Ungleichungsbedingung erfüllt!')
else
    disp('Inverse Ungleichungsbedingung NICHT erfüllt!')
end

%%
% %HYSTERESEBERECHNUNG
parfor i = 1:N
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end

% H_delta = w_H' * H_rH;

% parfor i = 1:(2*L+1)
%     S_rS(i,:) = fun_trajectory_SuperpositionsOperator(x, r_S(i));
% end

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
% K_delta = w_K' * K_rK;

% parfor i = 1:(2*L+1)
%     M_rM(i,:) = fun_trajectory_SuperpositionsOperator((H_delta + K_delta), r_S(i));
% end
% M_Kdelta = w_S' * M_rM;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% FEHLERMODELL          %
%%%%%%%%%%%%%%%%%%%%%%%%%

%Errorfunktion
w_H_tilde = w_H(2:end);
H_rH_tilde = H_rH(2:end,:);

%-S_rS_inv
parfor i = 1:(2*L+1)
    S_rS_inv(i,:) = fun_trajectory_SuperpositionsOperator(y, r_S_inv(i));
end

w_MK_tilde = [w_H_tilde; w_S_inv; w_K];
M_Kr = [H_rH_tilde; -S_rS_inv; K_rK];

% E_MK = x + (w_MK_tilde' * M_Kr);

H = zeros(length(w_MK_tilde), length(w_MK_tilde));
f = zeros(length(w_MK_tilde),1);
c = zeros(1,1);
parfor i=1:length(x)
H = H + M_Kr(:,i)*M_Kr(:,i)';
f = f + x(i)*M_Kr(:,i);
c = c + x(i) * x(i);
end

%%
% Constrains       
U_H = eye(N-1,N-1);

U_K = eye(N,N);
U_S = eye((2*L+1),(2*L+1));
for i=1:L+1
    parfor j=i:L+1
        U_S(i, j) = 1;
    end
end
for i=2*L+1:-1:L+1
    for j=i:-1:L+1
        U_S(i, j) = 1;
    end
end

u_H = zeros(N-1,1);
%u_H(1) = epsilon;
u_K = zeros(N,1);
u_S = epsilon*ones((2*L+1),1);

U_MK = zeros(N-1+N+(2*L+1), N-1+N+(2*L+1));
U_MK(1:N-1, 1:N-1) = U_H;
U_MK(N:N+(2*L+1)-1, N:N+(2*L+1)-1) = U_S;
U_MK(N+(2*L+1):N+(2*L+1)+N-1, N+(2*L+1):N+(2*L+1)+N-1) = U_K;

u_MK = [u_H; u_S; u_K];

%Parameter für Inequality constraints A*x <= b
%Vergleichen mit U_MK * w - u_MK >= 0
A = -U_MK;
b = -u_MK;

%%
%Minimierungsverfahren

%Lower bound

options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

[w_MK_tilde_ident,fval,exitflag,output,lambda] = quadprog(H,f,A,b,[],[],[],[],[],options);

fprintf('\nMaximale Abweichung eines Gewichts:     %d\n', abs(max(w_MK_tilde-w_MK_tilde_ident)))
fprintf('Quadratischer Fehler:                   %d\n', fval + c*0.5)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Hysterese mit identifizierten Parametern                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_H_ident = ones(length(w_H),1);
w_H_ident(2:end) = w_MK_tilde_ident(1:N-1);

w_S_inv_ident = ones(length(w_S),1);
w_S_inv_ident(1:end) = w_MK_tilde_ident(N:N+2*L);

w_K_ident = ones(length(w_K),1);
w_K_ident(1:end) = w_MK_tilde_ident(N+2*L+1:N+2*L+1+N-1);

w_S_ident = fun_w_S_inv(w_S_inv_ident, L);

r_S = fun_r_S_inv(r_S_inv, w_S_inv_ident, L);

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
t = linspace(0,T_s * length(K_rK),length(K_rK));
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

parfor i = 1:(2*L+1)
%     S_rS(i,:) = fun_trajectory_SuperpositionsOperator((H_delta), r_S(i));
    M_rM(i,:) = fun_trajectory_SuperpositionsOperator((H_delta + K_delta), r_S(i));
end
M_Kdelta = w_S_ident'*M_rM;

figure(1)
plot(x, M_Kdelta, 'r')
hold on
plot(x, y, 'b')
legend('Identifizierte Daten', 'Reelle Daten')