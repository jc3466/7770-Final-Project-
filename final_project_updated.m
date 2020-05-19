% This code models the change in free energy for an adhesion process between a coronavirus particle and a type II alveolar epithelial cell
% Adapted from Bell et al. "Cell Adhesion: Competition Between Nonspecific Repulsion and Specific Binding"
% Written By: Julia Caserto and Alexis Ostwalt
% CHEME 7770 Final Project 
%% 1 = cell, 2 = virus
%Thickness coefficient of cell glycocalyx
tau = 5; %nm
kappa = 4.28*10^-5; %kg/s^2
k = 1.380649*10^-5; %Boltzmann constant (nm^2*kg/ s^2*K)
T = 310.15; %temperature (K)
L = 10; % nm, unstrained length of cell-virus bridges

S_hat = 1/2*(L + ((k*T)/(kappa*tau)) + sqrt((L+((k*T)/(kappa*tau))).^2 + (4*k*T)/kappa));
A_set = pi.*((S_hat+tau)./2).^2; %Area of the set for contact region,depends on S
A1 = 3.1 * 10^8; %nm^2, cell surface area 
A2 = 4.5 * 10^4; %nm^2, virus surface area 
Amax = 10^9; % nm^2, Maximum contact area

% ---------------------------------------------------------------%
N1t = 10^7:-10^4:10^3; % total number of receptors per cell

% Spike proteins size
length = 30; % diameter, nm
width = 5; % thickness, nm 
spike_base = pi*(width/2)^2; % base area of spike 
N2t = A2 / spike_base; % virus, found using bionumbers thickness 
% assuming the base is circular, and the width is the diameter 

%Absolute number of bridges holding the cell and virus together
Nb = 50; 

%Local surface density of unattached receptors(1,2) and cell-virus bridges(b)
n_1 = (N1t - Nb)./A1;
n_2 = (N2t - Nb)./A2;
n_b = Nb./A_set;

%Compressibility coefficient of cell glycocalyx
gamma = 1.99*10^-2; %kg*nm/s^2
%Chemical potentials of free(1,2) or bound receptors(b) at some standard
%surface density
mu_10 = -13165e21; % nm^2/s^2
mu_20 = -13165e21; % nm^2/s^2
mu_b0 = -237.18 + 0.5.*kappa.*(S_hat-L).^2; %Taylors Thrm
mu_1 = mu_10 + k*T.*log(n_1);
mu_2 = mu_20 + k*T.*log(n_2);
mu_b = mu_b0 + k*T.*log(n_b);

Lambda = ((gamma./S_hat).*exp(-S_hat./tau));  % repulsive potential per unit area 

Ks = 1e7; % nm^3, from Bell 
Sigma = 10; % nm, from Bell 
KL = Ks./Sigma; % binding constant for formation of unstrained cell-virus bridges
KS = KL.*exp(-0.5.*kappa.*(S_hat-L).^2 ./ (k*T)); 
 
squiggle = (A1.*A2.*Lambda)./k*T.*KS; % this is squiggle for S_hat
line1 = N2t.*N1t;
plot(N1t,line1)
xlabel('N1t')
ylabel('squiggle')
legend('Region I-II Boundary')


