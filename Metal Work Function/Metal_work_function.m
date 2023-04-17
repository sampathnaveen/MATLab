% ------------------------------- Semiconductor theory II - Assignment_04 ------------------------------
% U88648766
% Malyadri Venkata Ssampath Naveen Padmanabhuni 
clc;
close all;
clear all;
% Inputs
silicon_type = "p";
N = 10^17; % doping density
tox = 2e-7; %oxide thickness (cm)
Vg = 2; % gate voltage (V)
Pox = 3.7*8.85e-14; %SiO2 permitivity
% Constants
eX = 4.05; %eV - electron affinity of slilcon
Q = 1e-8; % trapped charges
phiM = 4.06; %Work fuction of the metal
q = 1.609e-19; % Charge of the electron (Units: Coulombs)
E = 11.7*(8.85e-14); % Permittivity of silicon in free space (Units: Farad/cm)
%Oxide capacitance of given by Ci
Ci = Pox/tox;
Eb = 0.025*log(N/10^10); % Bulk potential
% Phi_s - Work function of the silicon.
if silicon_type == 'n'
    phiS = eX + 1.12/2 - abs(Eb);
else
    phiS = eX + 1.12/2 + abs(Eb);
end
%PhiMs
phiMS = phiM - phiS; %phiMS
vfb = phiMS - Q/Ci; %flatband voltage
Vg = Vg - vfb; %Vg voltage
%Now to calculate the Surface potential we need to calculate the Alpha
%and Beta values
beta = sqrt((E*q*N)/2)./Ci;
alpha = (-beta + sqrt(beta^2 + 4*Vg))/2;
Es = alpha^2; % Surface potential value
Wd = sqrt((2.*E.*abs(Es))./(abs(q).*N)); % width of depletion region
Wd_max = sqrt((4*E*Eb)/(q*N)) ; %Maximum depletion width
 if Wd > Wd_max 
     Wd = Wd_max;
 end
% Calculate Poissons Equation constants
a = (q*N)./(2*E);
b = (-q.*N*Wd)./(E); % b constant Poisson
c = Es;
% set x-axis length
x = linspace(0,2*Wd);
if Wd == 0
    x = linspace(0,0.2e-7);
end
 
Ps = []; %Potential Ps
Ev = []; %Valence band Ev

if silicon_type == 'p'
    for i = 1:length(x)
        if  x(i)<= Wd
            Ps(i) = a.*(x(i)^2) + b.*x(i) + c ; 
            Ev(i) = (-1)*Ps(i) ;  
            last_value = Ps(i);
        else
            Ps(i) = last_value;                     
            Ev(i) = -Ps(i);                  
        end
    end
    Eb = -Eb;
end
if silicon_type == 'n'
    for i = 1:length(x)
        if  x(i)<= Wd
            Ps(i) = -a*(x(i)^2) - b*x(i) - c ; 
            Ev(i) = (-1)*Ps(i) ;  
            last_value = Ps(i);
        else
            Ps(i) = last_value;                     
            Ev(i) = -Ps(i);                  
        end
    end
end

 x = (10000)*x;
    
Ec = Ev + (1.12); % Conduction band + energy gap
Ei = (Ev+Ec)/2; % Intrinsic energy level (Units : eV)
Ef = [];
Ef = Ei + Eb;
Ef_O = []; % constant value of Ef
    
for i = 1:length(Ef)
       Ef_O(i) = Ef(51);    % to make the fermi level constant
end 

xm = linspace(-30*tox,-tox,10);
xm = (10000)*xm;

for i = 1:length(xm)
    if Vg > 0 
        mpotential(i) = Ef(51) - Vg;
    elseif Vg < 0 
        mpotential(i) = Vg + abs(Ef(51));
    else
        mpotential(i) = Ef(51);
    end
end


    % Plot the energy band diagram/ with metal side 
    
         close all
        figure('units','normalized','outerposition',[0 0 1 1])
        title('Energy Band Diagram',"Color","black")
        xlabel('x(μm)','FontWeight','bold')
        ylabel('Energy(eV)','FontWeight','bold')
        hold
        plot(x,Ec,x,Ef_O,'--',x,Ei,x,Ev,xm,mpotential); % plot the curves togather in one plot
        xline(-tox*1e4)
        xline(0)
        hold

% calculate the minority and majority carrier concentration
   
        for i = 1:length(x)
           if silicon_type == 'p'
           Np(i) = (10^10)*exp((- Ei(i)+ Ef_O(i))/0.025);
           Nn(i) = 10^20/Np(i);
           else
           Nn (i) = (10^10)*exp((Ei(i)- Ef_O(i))/0.025);
           Np(i) = 10^20/Nn(i);
            end
        end
       
    % plot the minority and majority carrier concentration distribution
        %figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5])
        figure()
        semilogy(x,Np,x,Nn)
        title("Minority and Majority carriers")
        xlabel('x(μm)','FontWeight','bold')
        ylabel('Carrier concentration(atoms/cm3)')
        grid on
        hold
        legend(["Hole concentration","Electron concentration"],'FontSize',10,'TextColor','black')  