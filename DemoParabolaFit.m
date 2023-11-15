% Parabola fitting demo
%
% Please, reference this work as:
% 
% Ezequiel López-Rubio, Karl Thurnhofer-Hemsi, Elidia Beatriz Blázquez-Parra, Óscar David de Cózar-Macías, M. Carmen Ladrón-de-Guevara-Muñoz,
% A fast robust geometric fitting method for parabolic curves,
% Pattern Recognition, Volume 84, 2018, Pages 301-316, ISSN 0031-3203,
% https://doi.org/10.1016/j.patcog.2018.07.019.
%
% Last modification: 11/09/2018

clear all
close all

%% Seed inicialization
rng('default');

%% Parameter definition

% Square limits from where the points will be defined
UpperLimit=20;
%Number of candidate samples
NumCandidates=100000;
%Number of samples of the simulated parabola
NumSamples=50;
%Percentage of outliers added
percOut = 0.2;
%Penalization term of the fitting algorithm
Lambda=0.1;
%Step size of the fitting algorithm
StepSize=0.1;
%Limitations for the step size
MinStepSize=1.0e-4;
MaxStepSize=0.1;
%Maximum number of steps of the fitting algorithm
NumSteps=10000;
%Number of steps for the reinicialization of the focus and directrix
T=1000;

%% Generate points by Monte Carlo

%Create focus and directrix of the true parabola
TrueFocus=0.5*UpperLimit*rand(2,1);
TrueDirectrix=0.5*UpperLimit*rand(3,1);
TrueGeneralForm=ParabolaGeometric2General(TrueDirectrix,TrueFocus);

%Generate the sample candidates
SampleCandidates=-UpperLimit+2*UpperLimit*rand(2,NumCandidates);
% Compute the distances E_f(x,y) to the true focus
DifferencesFocus=SampleCandidates-repmat(TrueFocus,[1 NumCandidates]);
DistFocus=sqrt(sum(DifferencesFocus.^2,1));

% Compute the distances E_d(x,y) to the true directrix
LineEval=TrueDirectrix(1:2)'*SampleCandidates+TrueDirectrix(3); % a*x_i + b*y_i + c
DistDirectrix=abs(LineEval)/...
    sqrt(TrueDirectrix(1)^2+TrueDirectrix(2)^2);

% Compute the differences of the distances
DistDiff=(DistFocus-DistDirectrix).^2;

%Select those samples that are near to the true parabola
SamplesParabBig=SampleCandidates(:,DistDiff<2);
SamplesParab = SamplesParabBig(:,randperm(size(SamplesParabBig,2),NumSamples));

%Add outliers (if percOut>0)
SamplesNoParab=SampleCandidates(:,DistDiff>2);      
Samples=[SamplesParab SamplesNoParab(:,randperm(size(SamplesNoParab,2),round(percOut*NumSamples)))];


%% Execute the parabola fitting
tstart = tic;
[BestFocus,BestDirectrix,GeneralForm,Errors,Foci,Directrices,StepSizes,GradientVectors,~]=...
    ParabolaFit(Samples,Lambda,StepSize,MaxStepSize,MinStepSize,NumSteps,T);
CPUtime = toc(tstart);
Residuals = ComputeResidual(BestFocus,BestDirectrix,SamplesParab);


%% Plot results

fprintf('Execution time: %d\n',CPUtime);
fprintf('True general form: [%d %d %d %d %d %d]\n',TrueGeneralForm);
fprintf('Estimated general form: [%d %d %d %d %d %d]\n',GeneralForm);
fprintf('Residual error: %d\n',Residuals);


% Error evolution
figure
plot(Errors)
title('Errors')

%Step size evolution
figure
plot(StepSizes)
title('Step sizes')

% Fitted parabola
figure;
hold on
plot(Samples(1,:),Samples(2,:),'.','Color',[0 0 0],'MarkerSize',10);
plot(Foci(1,:),Foci(2,:),'-m')
Handle=ezplot(sprintf('%f*x^2+%f*x*y+%f*y^2+%f*x+%f*y+%f=0',TrueGeneralForm),[-UpperLimit UpperLimit]);
set(Handle,'Color',[255 200 60]/255,'LineStyle','-','LineWidth',2);
Handle=ezplot(sprintf('%f*x^2+%f*x*y+%f*y^2+%f*x+%f*y+%f=0',GeneralForm),[-UpperLimit UpperLimit]);
set(Handle,'Color',[1 0 0]);
plot(BestFocus(1),BestFocus(2),'sk','Color',[1 0 0])
DirectrixYest=polyval([-BestDirectrix(1)/BestDirectrix(2) -BestDirectrix(3)/BestDirectrix(2)],[-UpperLimit UpperLimit]);
line([-20 20],DirectrixYest,'Color',[1 0 0],'LineStyle','--');
axis([-UpperLimit UpperLimit -UpperLimit UpperLimit]);
title('Fitted parabola')
legend('Samples','Trace of the foci','True parabola','Fitted parabola','Best Focus','Best Directrix')





