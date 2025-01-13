function [omega,V_r,V_omega,V,M,T,p,rho,beta_c] = ...
    TaylorMaccoll(M_inf,T_inf,p_inf,a_inf,rho_inf,gamma,delta_c)

% Funzione per la soluzione del campo ipersonico conico, attraverso
% l'integrazione delle equazioni di Taylor-Maccoll con il metodo diretto.
% La funzione accetta in ingresso le condizioni asintotiche e l'angolo di
% deviazione della corrente in radianti, e restituisce le condizioni
% termofluidodinamiche nello shock layer. I valori di omega in output sono
% espressi in gradi.

% Dichiarazione di variabili globali
global V_lim n V_inf 

% Parametri del flusso
n = 2/(gamma - 1);
R = 287; % J/(kg*K), per l'aria
V_inf = M_inf*a_inf;
cp = gamma*R/(gamma-1);
h_inf = cp*T_inf;
H = V_inf^2/2 + h_inf;
V_lim = sqrt(2*H);
T0 = H/cp;
p0_inf = (T0/T_inf)^(gamma/(gamma-1))*p_inf;

Dbeta = 0.99; % Fattore riduttivo dato che beta_3D < beta_2D

% Calcolo dell'angolo d'urto per il cuneo (beta 2D valore partenza metodo iter.)
% Primo valore di tentativo per beta_c
beta_2D = AngoloUrto(M_inf,delta_c*180/pi,n);
% beta_2D = 14.42;
beta_c = convang(beta_2D,'deg','rad');

rapp_omega = 1;

% Metodo diretto di integrazione (assegno beta_c per ricavare delta_c)
while rapp_omega > 1e-6
    
    % Condizioni a valle dell'urto (condizioni iniziali integrazione)
    eps = 1/(n+1)*(1 + n/(M_inf^2*sin(beta_c)^2));
    V_r_0     = V_inf*cos(beta_c); % Vel. radiale
    V_omega_0 = -eps*V_inf*sin(beta_c); % Vel. tangenziale
    V0 = [V_r_0, V_omega_0];
    
  % Chiamata alla funzione evento per l'arresto dell'integrazione --> StopTM.m
  % Funzione che arresta la procedura di integrazione delle equazioni di
  % Taylor-Maccoll quando la componente di velocità V_omega si annulla.
    options = odeset('Event',@StopTM);
    
    % Integrazione equazioni di Taylor-Maccoll
    [omega,V] = ode45(@EquazioniTaylorMaccoll,[beta_c 0],V0,options);
    
    % Confronto dei valori di omega e delta_c
    rapp_omega = (omega(end) - delta_c)/delta_c;
    
    % Aggiornamento del valore di beta (con un valore ridotto)
    beta_c = beta_c*Dbeta;
end

beta_c = convang(beta_c/Dbeta,'rad','deg');

% Conversione di omega in gradi

omega = convang(omega,'rad','deg');

% Calcolo dei salti di grandezze attraverso l'urto
[p2p1,eps,~,~,~,p02p01] = RankineHugoniot(M_inf,beta_c,n);

% Calcolo delle grandezze immediatamente a valle dell'urto
p2  = p2p1*p_inf;
p0 = p02p01*p0_inf;
rho2 = rho_inf/eps;

% Estrazione vettori soluzione campo cinematico
V_r     = V(:,1);
V_omega = V(:,2);

% Modulo della velocità
V = sqrt(V_r.^2 + V_omega.^2);

% Calcolo grandezze termofluidodinamiche
T = (H - V.^2/2)/cp;
M  = sqrt(n*(T0./T - 1));
p = p0.*(T./T0).^(gamma/(gamma-1));
% p = p0/(1 + (gamma-1)/2.*M).^(gamma/(gamma-1));
rho = (p./p2).^(1/gamma)*rho2;

end


