clear all
clc
tic
%% PARAMETER

%Schranke für Ungleichbedingung
epsilon = 1e-6;

% Anzahl von Play-Operatoren +1
N = 80
% Halbe Anzahl (+1) von Superpositionsoperatoren
L = 40;

% Anzahl von elementaren Kriechoperatoren
M = 20;

% Abtastzeit (10 kHz)
T_s = 0.0001;

% Eingangssignal [V]
% load('x.mat')
% x = x-3;
% x = x';
x0 = 0;
deltaX = 1e-1;
xMax = 5;
xMin = -5;
x = [x0:deltaX:xMax xMax-deltaX:-deltaX:xMin xMin+deltaX:deltaX:xMax]; 

% Ausgangssignal [V]
% load('y.mat')
% y = y - 3.2;
% y = y(41:end);
load('y_sim.mat');
y = y_sim2;

% Initialisierung der Startwerte der Operatorausgänge
y_H0 = 0*ones(N+1, 1);
y_K0 = 0*ones(N+1, 1);

r_H = zeros(N+1,1);
parfor i=1:(N+1)
    r_H(i) = (i-1)/(N+1) * max(abs(x));
end

%NACH PESOTSKI:
r_S_inv = zeros((2*L+1),1);
parfor i=1:L
    r_S_inv(i) = (-(i-(L+1)))/(L+1)*min(y);
end
parfor i=(L+1):(2*L+1)
    r_S_inv(i) = ((i-(L+1)))/(L+1)*max(y);
end


% Schwellenwerte der Kriechoperatoren werden hier inisialisiert
r_K = zeros(N+1,1);
for i=1:(N+1)
    r_K(i) = (i-1)/(N+1) * max(abs(x));
end

% Matrizen der superpositionierten Operatoren
H_rH = zeros(N+1,length(x));
S_rS_inv = zeros((2*L+1),length(x));
K_rK = zeros(N+1,length(x));
M_rM = zeros((2*L)+1,length(x));

%%
% H_rH
parfor i = 1:N+1
    H_rH(i,:) = fun_trajectory_PlayOperator(x, r_H(i), y_H0(i));
end

% K_rK
a_Kj = zeros(length(M), 1);
parfor i=1:M
    a_Kj(i) = 1/(10^(i-1)*T_s);
end

parfor j = 1:(N+1)
    for i = 1:M
        y_kriech = fun_trajectory_KriechOperator(x, y_K0(j), a_Kj(i), r_K(j), T_s);
        
        K_rK(j,:) = K_rK(j,:) + y_kriech;
         if (i == M)
             K_rK(j,:) = K_rK(j,:)/M;
         end
    end
end

% -S_rS_inv
parfor i = 1:(2*L+1)
    S_rS_inv(i,:) = fun_trajectory_SuperpositionsOperator(y, r_S_inv(i));
end

%% FEHLERMODELL

H_rH_tilde = H_rH(2:end,:);

M_Kr = [H_rH_tilde; -S_rS_inv; K_rK];

H = zeros(length(N+(2*L+1)+N+1), length((N)+(2*L+1)+N+1));
f = zeros(length(N+(2*L+1)+N+1),1);
c = zeros(1,1);
parfor i=1:length(x)
H = H + M_Kr(:,i)*M_Kr(:,i)';
f = f + x(i)*M_Kr(:,i);
c = c + x(i) * x(i);
end

%%
% Constrains       
U_H = eye(N,N);

U_K = eye(N+1,N+1);
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
%u_H(1) = epsilon;
u_K = zeros(N+1,1);
u_S = epsilon*ones((2*L+1),1);

U_MK = zeros(N+N+1+(2*L+1), N+N+1+(2*L+1));
U_MK(1:N, 1:N) = U_H;
U_MK((N+1):N+(2*L+1), (N+1):N+(2*L+1)) = U_S;
U_MK((N+1)+(2*L+1):N+(2*L+1)+N+1, N+1+(2*L+1):N+(2*L+1)+N+1) = U_K;

u_MK = [u_H; u_S; u_K];

%Parameter für Inequality constraints A*x <= b
%Vergleichen mit U_MK * w - u_MK >= 0
A = -U_MK;
b = -u_MK;

%% OPTIMIERUNGSVERFAHREN

options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

[w_MK_tilde_ident,fval,exitflag,output,lambda] = quadprog(H,f,A,b,[],[],[],[],[],options);
if exitflag
   disp('Function converged to the solution') 
end

fprintf('Quadratischer Fehler:                   %d\n', fval + c*0.5)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Hysterese mit identifizierten Parametern                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_H_ident = ones(N+1,1);
w_H_ident(2:end) = w_MK_tilde_ident(1:N);

w_S_inv_ident = w_MK_tilde_ident(N+1:(N+1)+2*L);

w_K_ident = w_MK_tilde_ident((N+1)+2*L+1:(N+1)+2*L+N+1);

w_S_ident = fun_w_S_inv(w_S_inv_ident, L);

r_S = fun_r_S_inv(r_S_inv, w_S_inv_ident, L); %Richtig nach KUHNEN

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


K_delta = w_K_ident' * K_rK;

parfor i = 1:(2*L+1)
    M_rM(i,:) = fun_trajectory_SuperpositionsOperator((H_delta + K_delta), (r_S(i)));
end

M_Kdelta = w_S_ident'*M_rM;
toc
figure(1)
plot(x, M_Kdelta, 'r')
hold on
plot(x, y, 'b')
legend('Identifizierte Daten', 'Reelle Daten')


% S_rs_inv = zeros((2*L)+1,length(x));
% parfor i = 1:(2*L+1)
%     S_rs_inv(i,:) = fun_trajectory_SuperpositionsOperator((y), (r_S_inv(i)));
% end
% 
% S_delta_inv = w_S_inv_ident'*S_rS_inv;
% 
% H_rH_inv = zeros(N, length(x));
% r_H_inv = fun_r_H_inv(r_H, w_H_ident);
% y_H0_inv = fun_y_H0_inv(w_H_ident, y_H0);
% parfor i = 1:N
%     % Trajektorie des Playoperators
%     H_rH_inv(i,:) = fun_trajectory_PlayOperator((S_delta_inv-K_delta), r_H_inv(i), y_H0_inv(i));
% end
% w_H_inv = fun_w_H_inv(w_H_ident);
% 
% H_delta_inv = w_H_inv'*H_rH_inv;
% figure
% plot(x, H_delta_inv);
