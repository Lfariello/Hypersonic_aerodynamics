function beta = AngoloUrto(M1,theta,n)

% Funzione per il calcolo dell'angolo d'urto attraverso il metodo di
% Newton-Raphson. Gli angoli in input e in output sono in gradi.

% Conversione di theta in radianti
theta = convang(theta,'deg','rad');

% Valore di tentativo per beta (Tale da fornire un onda d'urto debole)
beta = theta*((n+1)/(2*n) + sqrt((n+1)^2/(4*n^2) + 1/(M1*theta)^2));

rapp = 1;

% Procedura iterativa per il calcolo dell'angolo d'urto
% F(beta) + (dF/dbeta)*D_beta = 0 ----> D_beta
while rapp > 1e-3
    F = tan(beta - theta)/tan(beta) - 1/(n+1)*(1 + n/(M1*sin(beta))^2);
    DF_Dbeta = 1/(cos(beta-theta)^2*tan(beta)) - ... % Sviluppo di Tayolr al 1 ordine
        tan(beta - theta)/sin(beta)^2 + ...
        2*n*cos(beta)/((n+1)*M1^2*sin(beta)^3);
    % Pongo lo sviluppo in serie uguale a zeero e ricavo:
    D_beta = - F/DF_Dbeta;  
    beta_new = beta + D_beta; % Valore di beta aggiornato   
    rapp = (beta_new - beta)/beta;    % Valuta se proseguire l'iterazione
    beta = beta_new; 
end

% Conversione di beta in gradi restituisco al main OndeDurto.m
beta = convang(beta,'rad','deg');

end