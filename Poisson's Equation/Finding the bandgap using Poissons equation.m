%% ------------------------------- Semiconductor theory II - Assignment_02 ------------------------------
% U88648766
% Malyadri Venkata Ssampath Naveen Padmanabhuni 

clc;
close all;
clear all;

%disp ("Please enter the input values as follows: ");
silicon_type = input ("Please enter the silicon type N or P = ", "s");

%fprintf (['If you have selected N type please provide the surface potential in negative value'
%'If you have selected P type, please provide the surface potential in positive value']);

Es = input ("Please enter the surface potential value 𝜓s(Units: electron Volts - eV) = ");
N = input ("Please enter the doping concentration (Units: cm^-3) = ");

%Simulation block

% ----------------------------------------  Now for the N - type  ----------------------------------------

if silicon_type == "N"

% These values will be the same entire code   
    q = 1.60217663*10e-19; % Charge of the electron (Units: Coulombs)
    E = 11.68*(8.85e-14); % Permittivity of silicon in free space (Units: Farad/cm)
    Wd = sqrt((2*E*Es)/(q*N)); % Width of the depletion region
   % Wd = sqrt((2*E*abs(Es))/(q*N)); % Width of the depletion region - with the ABS in the formulae
    x_max = linspace(0,2*Wd);
    Ef_Ei = (0.025)*log(N/1e10); %Energy difference between Ef and Ei
   

   % To determine if or not to include this block
   % if Es == 0
   %     x_max = linspace(0,0.2e-6);
   % end
    
    %Eg = 1.12; % Energy gap between Ec and Ev (Units: eV)
    %Ei = Eg/2; % Intrinsic energy level (Units : eV)
    
    % Now lets calculate the a,b,c values from the 2nd order differential
    % equation of the Poisson's equation - ax^2 + bx + c
    
    a = (-q*N)/(2*E); % derived from poisson equation
    b = ((q*N*Wd)/E); % derived from poisson equation
    c = -Es; % derived from poisson equation
    
    Px = []; % create an empty vector for the potential value
    Ev = [];
    
    for i = 1:length(x_max)
        if x_max(i) <= Wd
            Px(i) = a*(x_max(i)^2) + b*x_max(i) + c; 
            Ev(i) = (-1)*Px(i); % Calculating energy value
            Ev1 = Ev(i);
        else
            Px(i) = Ev1;
            Ev(i) = -Px(i);
        end
    end
    
    x_max = (10000)*x_max;
    Ec = Ev + (1.12);
    Ei = (Ev+Ec)/2; % Intrinsic energy level (Units : eV)
    
    Ef = [];
    Ef = Ei + Ef_Ei;
    
    for i = 1:length(Ef)
            Ef(i) = Ef(51);
    end
    
% ----------------------------------------  Now for the P - type  ----------------------------------------
    
elseif silicon_type == "P"

    q = 1.60217663*10e-19; % Charge of the electron (Units: Coulombs)
    E = 11.68*(8.85e-14); % Permittivity of silicon in free space (Units: Farad/cm)
    Wd = sqrt((2*E*Es)/(q*N)); % Width of the depletion region
    % Wd = sqrt((2*E*abs(Es))/(q*N)); % Width of the depletion region - with the ABS in the formulae
    x_max = linspace(0,2*Wd);
    Ef_Ei = (0.025)*log(N/1e10); %Energy difference between Ef and Ei
    
   % if Es == 0
   %     x_max = linspace(0,0.2e-6);
   % end
    
    %Eg = 1.12; % Energy gap between Ec and Ev (Units: eV)
    %Ei = Eg/2; % Intrinsic energy level (Units : eV)
    % Now lets calculate the a,b,c values from the 2nd order differential
    % equation of the Poisson's equation - ax^2 + bx + c
    
    a = (q*N)/(2*E); % derived from poisson equation
    b = ((-q*N*Wd)/E); % derived from poisson equation
    c = Es; % derived from poisson equation
    Px = []; % create an empty vector for the potential value
    Ev = [];
    
    for i = 1:length(x_max)
        if x_max(i) <= Wd
            Px(i) = a*(x_max(i)^2) + b*x_max(i) + c; 
            Ev(i) = (-1)*Px(i); % Calculating energy value
            Ev1 = Ev(i);
        else
            Px(i) = Ev1;
            Ev(i) = -Px(i);
        end
    end
    
    x_max = (10000)*x_max;
    Ec = Ev + (1.12);
    Ei = (Ev+Ec)/2; % Intrinsic energy level (Units : eV)
    Ef = [];
    Ef = Ei - Ef_Ei;
    
    for i = 1:length(Ef)
            Ef(i) = Ef(51);
    end
end

%plot(x_max,Ec,x_max,Ei,x_max,Ef,'--',x_max,Ev) % Energy Band Plot

%plot(x_max,Ec,x_max,Ei,x_max,Ef,x_max,Ev) % Energy Band Plot
%legend('Ec','Ei','Ef','Ev')
%xlabel('Distance (um)')
%ylabel('Energy Level (eV)')

% --------------------------------------------------------- Energy Band diagram plot graph ---------------------------------------------------------

figure
    plot(x_max,Ev)
    title('Energy Band Diagram')
    xlabel('Wd(um)')
    ylabel('Energy(eV)')
    hold
    plot(x_max,Ec,x_max,Ei,x_max,Ef,'--' )
    hold