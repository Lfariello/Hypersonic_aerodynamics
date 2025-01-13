%% CUNEO - STUDIO GRANDEZZE TERMOFLUIDODINAMICHE
% Confronto campo conico (3D) vs campo generato da un cuneo (2D)
% Studio campo di Taylor-Maccoll all'interno dello strato limite
% supersonico di un cono posto a zero angolo d'attacco in una corrente
% ipersonica. Si confrontano i risultati ottenuti con quelli del cuneo per
% le stesse condizioni asintotiche
clc
close all
clear all

% Condizioni asintotiche (manuale)
M_inf = 8;
gamma = 1.4;    % Valori dell'aria ad una data quota
T_inf = 273;
p_inf = convpres(0.1,'atm','Pa');
a_inf = 331;
rho_inf = 0.129;

% Parametri del flusso
n = 2/(gamma - 1); % (= 5)
R = 287; % J/(kg*K), per l'aria
V_inf = M_inf*a_inf; % (= 2648)
cp = gamma*R/(gamma-1);
h_inf = cp*T_inf;
H = V_inf^2/2 + h_inf;
%V_lim = sqrt(2*H);
T0 = H/cp;
p0_inf = (T0/T_inf)^(gamma/(gamma-1))*p_inf;

% Geometria --> Angolo di semiapertura del cono (e cuneo)
delta_c = deg2rad(15); % Angolo di deflessione imposto alla corrente

% Calcolo il campo di moto (condizioni termofluidodinamiche nello shock
% layer)
[omega,V_r,V_omega,V,M,T,p,rho,beta_c] = TaylorMaccoll(M_inf,T_inf,p_inf,a_inf,rho_inf,gamma,delta_c);
% La funzione restituisce vettori (13x1) che rappresentano le variazioni
% delle grandezze T.F.D muovendosi lungo lo spessore dello shock layer in
% piu il valore dell'angolo d'urto

% Plot campo di moto
omega = convang(omega, 'deg','rad');
fprintf('beta_cono = %f\n',beta_c);
beta_c = convang(beta_c, 'deg','rad'); % (= 17.8 << beta_2D = 20.9)
N = 8; % Riduco i vettori (13x1 --> 8x1)
omega_fitto = linspace(omega(1),omega(end),N);
V_omega = interp1(omega,V_omega, omega_fitto);
V_r = interp1(omega,V_r, omega_fitto);

% Estendo i risultati su tutto il cono
Lc = 1; % Lunghezza
Yc = tan(delta_c)*Lc; % Superficie laterale
% Plot del cono
plot([0,Lc],[0,Yc],'k',[0,Lc],[0,-Yc],'k')
hold on
Yb = tan(beta_c)*Lc;
% Plot onda d'urto
plot([0,Lc],[0,Yb],'r',[0,Lc],[0,-Yb],'r')


x = fliplr(linspace(0.001,Lc,50));
omega_fitto = fliplr(omega_fitto);
V_r = fliplr(V_r);
V_omega = fliplr(V_omega);
y = x.*tan(omega_fitto(:));
X = x;
Y = y;
for i = 1:N-1
    X = [X;x];
end

x = X;
r = flipud(sqrt(y.^2+x.^2));
Vx = V_r(:).*x./r+flipud(abs(V_omega(:)).*y./r);
Vy = flipud(V_r(:).*y./r)-flipud(abs(V_omega(:)).*x./r);
VX = Vx;
VY = Vy;

% pcolor(X,-Y,hypot(VX,-VY))
% shading interp
% Y = flipud(Y);
% X = flipud(X);

quiver(X,Y,VX,VY,'b','AutoScale','on', 'AutoScaleFactor',0.2)
x = linspace(0,1,10);
plot(x,0*x,'--k')
axis equal
ylim([-0.32 0.32])
xlim([0 1])

quiver(X,-Y,VX,-VY,'b','AutoScale','on', 'AutoScaleFactor',0.2)
x = linspace(0,1,10);
plot(x,0*x,'--k')
axis equal
ylim([-0.32 0.32])
xlim([0 1])

% Confronto con cuneo

% Calcolo dei salti di grandezze attraverso l'urto
beta_2D = AngoloUrto(M_inf,delta_c*180/pi,5);
[p2p1,eps,~,~,T2T1,p02p01] = RankineHugoniot(M_inf,beta_2D,5);

p2_cuneo = p2p1*p_inf;
p2_cuneo_v = linspace(p2_cuneo,p2_cuneo,13)';

rho_cuneo = (1./eps)*rho_inf;
rho_cuneo_v = linspace(rho_cuneo,rho_cuneo,13)';

T_cuneo = T2T1*T_inf;
T_cuneo_v = linspace(T_cuneo,T_cuneo,13)';


% Plots
figure
plot(rad2deg(omega),V,'-ok','LineWidth',1);
grid on;
title('Velocita totale vs omega')
ylabel('V (m/s)');
xlabel('\omega (deg)');

figure
plot(rad2deg(omega_fitto),V_r,'-ok','LineWidth',1);
grid on;
title('Velocita radiale vs omega')
ylabel('V_r (m/s)');
xlabel('\omega (deg)');

figure
plot(rad2deg(omega_fitto),V_omega,'-ok','LineWidth',1);
grid on;
title('Velocita angolare vs omega')
ylabel('V_{\omega} (m/s)');
xlabel('\omega (deg)');

figure
plot(rad2deg(omega),p,'-ok','LineWidth',1);
grid on;
hold on;
plot(rad2deg(omega),p2_cuneo_v,'-ob','LineWidth',1)
title('Pressione vs omega')
ylabel('P (Pa)');
xlabel('\omega (deg)');
legend('Cono','Cuneo','-ok','ob');

figure
plot(rad2deg(omega),T,'-ok','LineWidth',1);
grid on;
hold on;
plot(rad2deg(omega),T_cuneo_v,'-ob','LineWidth',1)
title('Temperatura vs omega')
ylabel('T (K)');
xlabel('\omega (deg)');
legend('Cono','Cuneo','-ok','ob');

figure
plot(rad2deg(omega),rho,'-ok','LineWidth',1);
grid on;
hold on;
plot(rad2deg(omega),rho_cuneo_v,'-ob','LineWidth',1)
title('Densita vs omega')
ylabel('\rho (Kg/m^{3})');
xlabel('\omega (deg)');
legend('Cono','Cuneo','-ok','ob');

