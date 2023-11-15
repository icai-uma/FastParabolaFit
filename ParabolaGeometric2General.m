function [GeneralForm]=ParabolaGeometric2General(Directrix,Focus)
% Convert the geometric form of a parabola in 2D 
% given by the directrix:
% ax + by + c = 0
% and the focus (u,v), to the general form:
% Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
% Created by Ezequiel López-Rubio and Karl Thurnhofer-Hemsi
% Last modification: 10/06/2018

a=Directrix(1);
b=Directrix(2);
c=Directrix(3);
u=Focus(1);
v=Focus(2);


GeneralForm=[ ((a^2)/(a^2+b^2))-1 ...
    2*a*b/(a^2+b^2) ...
    ((b^2)/(a^2+b^2))-1 ...
    2*u+2*a*c/(a^2+b^2) ...
    2*v+2*b*c/(a^2+b^2) ...
    ((c^2)/(a^2+b^2))-u^2-v^2];

GeneralForm=GeneralForm/norm(GeneralForm);
if GeneralForm(1)<0
    GeneralForm=-GeneralForm;
end



