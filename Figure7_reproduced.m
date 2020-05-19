% This code models the change in free energy for an adhesion process between a coronavirus particle and a type II alveolar epithelial cell
% Adapted from Bell et al. "Cell Adhesion: Competition Between Nonspecific Repulsion and Specific Binding"
% Written By: Julia Caserto and Alexis Ostwalt
% CHEME 7770 Final Project 

%% 1 = cell, 2 = virus

%%--Glycocalyx Parameters------------------------------------------------------------------------%%
tau = 5; % Thickness coefficient of cell glycocalyx, nm
Ns1 = 0.5e14; % # polymer segments/unit area of glycocalyx (saccaride residues/nm^2),Bell et al. 
Vs = 3e-22; % cm^3, volume of a typical saccaride moeity of the glycocalyx,Value from Bell et al. 
Vw = 3e-23; % cm^3, Bionumbers: 106549 volume of a water molecule
%%------------------------------------------------------------------------------------------------%%

kappa = 0:.01:14; % kappa range,(kg/s^s)
k = 1.380649*10^-5; %Boltzmann constant (nm^2*kg/ s^2*K)
T = 310.15; % temperature (K), this is body temperature
L = 10; % nm, unstrained length of cell-virus bridges
%Separation distance between cell glycocalyx and virus
%  S = .01:0.01:10;  %From 0 to 2tau (nm)
% S_hat = 1/2*(L + ((k*T)./(kappa.*tau)) + sqrt((L+((k*T)./(kappa.*tau))).^2 + (4*k*T)./kappa));
S_hat = 10;
A_set = pi.*((S_hat+tau)./2).^2; % Area of the set for contact region,depends on S
A1 = 3.1 * 10^8; % nm^2, surface area alveolar type II epithelial cells 
A2 = 4.5 * 10^4; % nm^2, surface area virus
Amax = 10^9; % maximum contact area, nm^2

% ------------------------------------------------------------------------------------------------%%

N1t = 10^4; % total receptor # per cell

% Dimenstions of spike proteins 
length = 30; % diameter, nm
width = 5; % thickness, nm 
spike_base = pi*(width/2)^2; % base area of spike 
N2t = A2 / spike_base; % assuming the base is circular, and the width is the diameter 

gamma = kappa*T*(Vs^2/Vw)*Ns1^2./10^7; % compressibility coeff. glycocalyx, kg*nm/s^2
Lambda = ((gamma./S_hat).*exp(-S_hat./tau));  % Repulsive potential per unit area 

Ks = 1e7; % nm^3, from Bell 
Sigma = 10; % nm, from Bell 
KL = Ks./Sigma; % binding constant for formation of unstrained cell-virus bridges
KS = KL.*exp(-0.5.*kappa.*(S_hat-L).^2 ./ (k*T)); 
squiggle = (A1.*A2.*Lambda)./k*T.*KS; % this is squiggle for S_hat


Nb_hat = 0.5.*(N2t-(sqrt((N2t.^2)+(4*squiggle))/1e11));  % take out N1t because virus is smaller, the limiting binding number thing 
N_new = Nb_hat./N2t; % use N2t because it is the smallest and therefore the max # of bridges possible, dimensionless

scatter(gamma,N_new,'m')   % plot fraction of spikes bound vs. gamma
ax = gca;
ax.FontSize = 16;
xlabel('Gamma (kg*nm/s^2)','FontSize',20)
ylabel('Fraction of Spikes Bound to Cell Surface','FontSize',20)
title('Fraction of Spikes Bound to Cell Surface vs. Cell Glycocalyx Compressibility','FontSize',25)
grid on 

