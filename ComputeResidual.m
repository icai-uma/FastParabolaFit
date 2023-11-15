% Compute the residual error between the true and the estimated parabola
% Created by Ezequiel López-Rubio and Karl Thurnhofer-Hemsi
% Last modification: 10/06/2018
function Residual = ComputeResidual(Focus,Directrix,Samples)

NumSamples = size(Samples,2);

% Compute the distances E_f(x,y) to the true focus
DifferencesFocus=Samples-repmat(Focus,[1 NumSamples]);
DistFocus=sqrt(sum(DifferencesFocus.^2,1));

% Compute the distances E_d(x,y) to the true directrix
LineEval=Directrix(1:2)'*Samples+Directrix(3); % a*x_i + b*y_i + c
DistDirectrix=sqrt((LineEval.^2)/...
    (Directrix(1)^2+Directrix(2)^2));

% Compute the differences of the distances
Residual=mean((DistFocus-DistDirectrix).^2);
