function [BestFocus,BestDirectrix,GeneralForm,Errors,Foci,Directrices,StepSizes,GradientVectors,Penalty]=ParabolaFit(Samples,Lambda,StepSize,MaxStepSize,MinStepSize,NumSteps,T)
% Fit a parabola to data, enhanced version, with further simplifications,
% normalization of the equation of the directrix at each step and in the 
% inicialization, and random restart each 1000 steps.
% Gradient descent with adaptive step size. The best solution is returned
% rather than the last one.
% Positive values of the Lambda parameters are recommended because this way
% the solutions with the focus very close to the directrix are avoided.
% Inputs:
%   Samples=Matrix of size 2 x NumSamples with the training samples
% Outputs:
%   BestFocus=Vector of size 2 x 1 with the focus point of the fitted
%       parabola.
%   BestDirectrix=Vector of size 3 x 1 with the coefficients of the
%       directrix (in general form, Mx + Ny + P = 0) of the fitted parabola.
%   GeneralForm=Vector of size 6 x 1 with the coefficients of the
%       fitted parabola in general form, Ax^2 + Bxy + Cy^2 + Dx + Ey + F =
%       0.
% Created by Ezequiel López-Rubio and Karl Thurnhofer-Hemsi
% Last modification: 10/06/2018

[Dimension, NumSamples]=size(Samples);

% Initialize the focus to the mean of the distribution and the
% directrix to the linear regression line of the distribution
StdDev=mean(std(Samples,0,2));
% BaseFocus=mean(Samples,2)+randn(2,1); % [u v]
BaseFocus=mean(Samples,2)+0.5*std(Samples,0,2).*randn(2,1); % [u v]
%Focus=[    8.1472; 9.0579];
LinPoly=polyfit(Samples(1,:),Samples(2,:),1);
BaseDirectrix=[LinPoly(1) -1 LinPoly(2)]'; % [a b c]
BaseDirectrix=BaseDirectrix/sqrt(BaseDirectrix(1)^2+BaseDirectrix(2)^2);
%BaseDirectrix=[1 0 -10]';
%Directrix=[1.2699; 9.1338; 6.3236];

% Initialize log variables
Errors=zeros(1,NumSteps);
Penalty=zeros(1,NumSteps);
Directrices=zeros(3,NumSteps);
Foci=zeros(2,NumSteps);
StepSizes=zeros(1,NumSteps);
StepSizes(1)=StepSize;
GradientVectors=zeros(Dimension,NumSteps);

% Main loop
BestError=1.0e300;
Gradient=zeros(5,1);
Focus=BaseFocus;
Directrix=BaseDirectrix;
for NdxStep=1:NumSteps
    % Compute the distances E_f(x,y) to the current focus 
    DifferencesFocus=Samples-repmat(Focus,[1 NumSamples]);
    DistFocus=sqrt(sum(DifferencesFocus.^2,1));
    
    % Compute the distances E_d(x,y) to the current directrix
    LineEval=Directrix'*[Samples; ones(1,NumSamples)]; % a*x_i + b*y_i + c
    NormDirectrix=Directrix(1)^2+Directrix(2)^2;
    DistDirectrix=sqrt((LineEval.^2)/...
        NormDirectrix); % a^2+b^2
    EvalDirectrix=Directrix'*[Focus;1]; % a*u+b*v+c
    
    % Find the points inside and outside the ellipsoid
    FocusPoints=DistFocus<DistDirectrix; % R_f
    DirectrixPoints=~FocusPoints; % R_d
%     NumInnerPoints=nnz(FocusPoints);
%     NumOuterPoints=nnz(DirectrixPoints);
    
    % Find the gradient of the parameter vector p
    AuxVec=sign(LineEval).*...
        (Directrix(1)*Directrix(3)+Directrix(1)*Directrix(2)*Samples(2,:) ...
        -Directrix(2)^2*Samples(1,:)) ./ ...
        (NormDirectrix^(3/2));
    PenaltyTerm=sign(EvalDirectrix)*...
        (Directrix(1)*Directrix(3)+Directrix(1)*Directrix(2)*Focus(2) ...
        -Directrix(2)^2*Focus(1)) /...
        ((NormDirectrix^(3/2)));
    Gradient(1)=(sum(AuxVec(DirectrixPoints)) ...
        -sum(AuxVec(FocusPoints)))/NumSamples ...
        +Lambda*PenaltyTerm;
    
    AuxVec=sign(LineEval).*...
        (Directrix(2)*Directrix(3)+Directrix(1)*Directrix(2)*Samples(1,:) ...
        -Directrix(1)^2*Samples(2,:)) ./ ...
        (NormDirectrix^(3/2));
    PenaltyTerm=sign(EvalDirectrix)*...
        (Directrix(2)*Directrix(3)+Directrix(1)*Directrix(2)*Focus(1) ...
        -Directrix(1)^2*Focus(2)) /...
        ((NormDirectrix^(3/2)));
    Gradient(2)=(sum(AuxVec(DirectrixPoints)) ...
        -sum(AuxVec(FocusPoints)))/NumSamples ...
        +Lambda*PenaltyTerm;

    AuxVec=sign(LineEval) ./ ...
        ((Directrix(1)^2+Directrix(2)^2)^(3/2));
    PenaltyTerm=sign(EvalDirectrix) /...
        ((NormDirectrix^(3/2)));
    Gradient(3)=(-sum(AuxVec(DirectrixPoints)) ...
        +sum(AuxVec(FocusPoints)))/NumSamples ...
        -Lambda*PenaltyTerm;
    
    AuxVec=-DifferencesFocus(1,:)./DistFocus;
    PenaltyTerm=Directrix(1)*sign(EvalDirectrix) /...
        ((NormDirectrix^(3/2)));
    Gradient(4)=(sum(AuxVec(DirectrixPoints)) ...
        -sum(AuxVec(FocusPoints)))/NumSamples ...
        -Lambda*PenaltyTerm;
    
    AuxVec=-DifferencesFocus(2,:)./DistFocus;
    PenaltyTerm=Directrix(2)*sign(EvalDirectrix) /...
        ((NormDirectrix^(3/2)));
    Gradient(5)=(sum(AuxVec(DirectrixPoints)) ...
        -sum(AuxVec(FocusPoints)))/NumSamples ...
        -Lambda*PenaltyTerm;
    
    % Update the best solution found so far
    Errors(NdxStep)=(sum(DistDirectrix(FocusPoints)-DistFocus(FocusPoints))+...
        sum(DistFocus(DirectrixPoints)-DistDirectrix(DirectrixPoints)))/NumSamples-...
        Lambda*sqrt(((Directrix'*[Focus;1])^2)/NormDirectrix);
    Penalty(NdxStep)=Lambda*sqrt(((Directrix'*[Focus;1])^2)/NormDirectrix);
    if Errors(NdxStep)<BestError
        BestError=Errors(NdxStep);
        BestFocus=Focus;
        BestDirectrix=Directrix;
    end
       
    % Update the directrix and normalize it
    Directrix=Directrix-StepSize*Gradient(1:3);
    Directrix=Directrix/sqrt(Directrix(1)^2+Directrix(2)^2);
    
    % Update the focus
    Focus=Focus-StepSize*Gradient(4:5);
    
    if mod(NdxStep,T)==0
        Focus=BaseFocus+StdDev*0.5*randn(size(Focus));
        Directrix=BaseDirectrix+StdDev*0.5*randn(size(Directrix));
    end
    
%     if mod(NdxStep,1000)==0
%         Directrix=10*randn(size(Directrix));
%     end
%     
%     if mod(NdxStep,1000)==500
%         Focus=10*randn(size(Focus));
%     end
    
    % Update the logs
    Directrices(:,NdxStep)=Directrix;
    Foci(:,NdxStep)=Focus;
    
    % Update the step size each 10 iterations of the main loop
    if mod(NdxStep,10)==0
        % Check whether the error has grown in a robust way
        if median(Errors((NdxStep-4):NdxStep))<median(Errors((NdxStep-9):(NdxStep-5)))
            % The error is smaller, so we increase the step size provided
            % that it is not too big
            if StepSize<MaxStepSize
                StepSize=1.1*StepSize;
            end
        else
            % The error is bigger, so we decrease the step size provided
            % that it is not too small
            if StepSize>MinStepSize
                StepSize=0.9*StepSize;
            end
        end
    end
    StepSizes(NdxStep)=StepSize;
end

a=BestDirectrix(1);
b=BestDirectrix(2);
c=BestDirectrix(3);
u=BestFocus(1);
v=BestFocus(2);
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


