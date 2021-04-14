function [ R, V ] = OrbitalElements2rvGeo( h, muo, e, i, Omega, w, theta )
% TELETIN ALEXANDRA-IONELA 
% GRUPA MAE1
% 15.04.2021
% INPUTS:
% h              : vectorul momentului unghiular in km^2/s
% muo        :  Constanta Gravitationala
% e              : excentricitatea
% i               : inclinatia
% Omega    : longitudinea nodului ascendent
% w              : argumentul periheliului
% theta        : anomalia adevarata
% OUTPUTS:
% r    : vectorul de pozitie
% v    : vecctorul viteza in km/s
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
%r
R_bar=h^2/muo/(1+e*cosd(theta))*[cosd(theta);sind(theta);0];
% v
V_bar=muo/h*[-sind(theta);e+cosd(theta);0];
QxX=[-sind(Omega)*cosd(i)*sind(w)+cosd(Omega)*cosd(w),-sind(Omega)*cosd(i)*cosd(w)-cosd(Omega)*sind(w),sind(Omega)*sind(i);...
           cosd(Omega)*cosd(i)*sind(w)+sind(Omega)*cosd(w),cosd(Omega)*cosd(i)*cosd(w)-sind(Omega)*sind(w),-cosd(Omega)*sind(i);...
           sind(i)*sind(w),sind(i)*cosd(w),cosd(i)];
R=QxX*R_bar;
V=QxX*V_bar;
end