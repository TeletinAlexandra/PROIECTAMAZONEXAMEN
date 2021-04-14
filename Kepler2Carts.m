function [X, Y, Z, Vx, Vy, Vz] = Kepler2Carts(alt, ecc, inc, w, nu, RAAN)
%--------------------------------------------------------------------------------------------------------%
%
% 			Convertire din elemente orbitale in sistem OXYZ
%
% 			Autor: Teletin Alexandra-Ionela
%
% 			DATA: 15.04.2021
%
% 			INPUT:
% 			alt:    altitudine.....................(Km)							
% 			ecc:    excentricitate											    
% 			inc:	inclinatie..................(rad)							
% 			w:	    argumentul perigeului..........(rad)	
% 			nu:	    pozitia satelitului...........(rad)							
% 			RAAN:	Longitudinea nodului ascendent............(rad)							
%
% 			OUTPUT:
%			Componentele vectorului de pozitie: 		
% 			[X Y Z]...(Km)
%
%			Componentele vectorului viteza:
% 			[Vx Vy Vz]...(Km/s) 														
%
%
%%---------------------------- Constante ----------------------------------------------%
mu_earth = 3.986 * 10^5; % Constanta gravitationala a Pamantului
re = 6378.1; % raza Pamantului
%
%%--------------------------------------------------------------------------------------
a = alt + re;
p = a*(1-ecc ^2);
r_0 = p / (1 + ecc * cos(nu));
%
%%--------------- Coordonatele in sistemul de referinta Oxyz -----------------%
%
% coordonatele vectorului de pozitie
x = r_0 * cos(nu);
y = r_0 * sin(nu);
%
%
% coordonatele vectorului viteza
Vx_ = -(mu_earth/p)^(1/2) * sin(nu);
Vy_ = (mu_earth/p)^(1/2) * (ecc + cos(nu));
%
%
%%--------------  OXYZ ---------------------%
%dupa aplicarea celor doua rotatii avem
%  X, Y, and Z
X = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * x + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * y;
Y = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * x + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * y;
Z = (sin(w) * sin(inc)) * x + (cos(w) * sin(inc)) * y;
%  X', Y', and Z'
Vx = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * Vx_ + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * Vy_;
Vy = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * Vx_ + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * Vy_;
Vz = (sin(w) * sin(inc)) * Vx_ + (cos(w) * sin(inc)) * Vy_;
