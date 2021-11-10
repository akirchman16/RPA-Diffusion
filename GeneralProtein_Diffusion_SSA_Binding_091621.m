clearvars;
close all;

%This code will simply describe a general protein diffusing around on a
%1-dimensional lattice. It will run a SSA calculation for binding and
%unbinding, but with each step a simple Monte-Carlo method will be used to
%test for diffusion in one direction or another if possible. The diffusion
%will essentially be a random walk simulation for the proteins

%Stochastic Binding w/ Deterministic Diffusion

N = 5000;    %length of DNA lattice
n = 3;  %length of each protein

k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 1;      %cooperativity parameter
L_A = 1;  %concentration of free proteins

Diffusion_Rate = 100000;   %rate of diffusion for proteins (diffusion events/time unit)
Left_Prob = 0.5;    %probability of a protein diffusing right

Iterations = 10000;  %number of events which will occur

K = k_on/k_off; %equilibrium constant
Locations = 2:1:N-(n-1)+1;  %list of all possible locations on the lattice

DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros

xB = zeros(1,Iterations+1);
xAB = zeros(1,Iterations+1);
t = zeros(1,Iterations+1);
a_f = zeros(1,Iterations);
a_r = zeros(1,Iterations);
a_0 = zeros(1,Iterations);
dt = zeros(1,Iterations);
Hist = zeros(2,Iterations); %top row is binding event, bottom row is unbinding events
ProteinCount = zeros(1,Iterations+1);
FracCover = zeros(1,Iterations+1);
BindProb = zeros(Iterations,N); %probability of protein binding to each location if a binding event were to occur next
                                %rows are DNA lattices after each iteration
CoopEvents = zeros(3,Iterations);   %top row is isolated events, second row is singly contiguous, third row is doubly contiguous
ProteinTracking = zeros(Iterations+1,N);  %memory allocation to track where proteins are bound over time
DiffusedProteins = zeros(Iterations,N); %memory allocation to record which proteins diffused at which event number (initial location of protein)

t(1) = 0;
BindCounter = 0;
UnbindCounter = 0;
TotalDiffEvent = 0;
TotalLeftDiff = 0;
TotalRightDiff = 0;
MaxOut = 0;

xB(1) = N-(n-1);    %initial values for a free lattice
xAB(1) = 0;

ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice

LeftDiffCounter = 0;
RightDiffCounter = 0;
Event = 0;  %counts number of events which have occured
Equilibrium = 0;    %Is the system in equilibrium? 0 = False, 1 = True

tic
while Equilibrium == 0
    Event = Event+1;
    a_f(Event) = k_on*(L_A)*(xB(Event));   %propensity functions (probability of each event happening)
    a_r(Event) = k_off*(xAB(Event));
    a_0(Event) = a_f(Event)+a_r(Event);     %sum of propensity functions used for determining dt

    dt(Event) = (1/a_0(Event))*log(1/rand); %random time interval for Gillespie method
    
    for j = 2:N-(n-1)+1   %loop to calculate probability of binding at each spot if binding event occurs
        if (DNA(j) ~= 1) && (sum(DNA(j:j+(n-1))) == 0)
            if DNA(j-1) == 0 && DNA(j+n) == 0   %checks for isolated location
                BindProb(Event,j-1) = 1/(xB(Event));
            elseif (DNA(j-1) == 0 && DNA(j+n) == 1) || (DNA(j-1) == 1 && DNA(j+n) == 0) %checks for singly contiguous location
                BindProb(Event,j-1) = w/(xB(Event));
            elseif DNA(j-1) == 1 && DNA(j+n) == 1 %checks for doubly contiguous location
                BindProb(Event,j-1) = (w^2)/(xB(Event));
            end
        end
    end

    if a_f(Event) > rand*a_0(Event) %a forward reaction occurs
        eventB = 0;
        while ~eventB   %repeats until a binding occurs
            SpotB = randsample(Locations,1,true,BindProb(Event,1:N-n+1));  %random location on lattice is chosen
            if DNA(SpotB:SpotB+(n-1)) == 0   %checks if location is free
               DNA(SpotB:SpotB+(n-1)) = 1;   %binds protein to location
               BoundAtSpot(SpotB) = 1;  %stores locatin in BoundAtSpot
               BindCounter = BindCounter+1;    %updates bind counter
               Hist(1,Event) = SpotB;      %stores location of binding event

               FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %function to find number of free spots now available
               xB(Event+1) = FreeSpots;    %updates populations
               xAB(Event+1) = xAB(Event)+1;

               eventB = 1;
            end
        end
        if DNA(SpotB-1) == 0 && DNA(SpotB+n) == 0 %isolated event occured
            CoopEvents(1,Event) = SpotB;
        elseif (DNA(SpotB-1) == 0 && DNA(SpotB+n) == 1) || (DNA(SpotB-1) == 1 && DNA(SpotB+n) == 0) %singly contiguous event occured
            CoopEvents(2,Event) = SpotB;
        elseif DNA(j-1) == 1 && DNA(j+n) == 1   %doubly contiguous binding event occured
            CoopEvents(3,Event) = SpotB;
        end
    else                       %otherwise an unbinding has to occur
        CurrentBound = find(BoundAtSpot == 1);
        pos = randi(length(CurrentBound));
        SpotU = CurrentBound(pos);  %random position for protein to unbind

        DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
        BoundAtSpot(SpotU) = 0; %removes location from BoundAtSpot
        UnbindCounter = UnbindCounter+1;    %updates unbind counter
        Hist(2,Event) = SpotU;  %stores location of unbinding event

        FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %calculates number of free spaces
        xB(Event+1) = FreeSpots;
        xAB(Event+1) = xAB(Event)-1;
    end
%%% Diffusion Process - Deterministic Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CurrentBound = find(BoundAtSpot == 1);  %list of all locations where a protein is bound
    DiffOrder = randperm(numel(CurrentBound)); %random order to test for diffusion
    DiffusionEvents = 0;    %initially no diffusion events have occured in this time step
    CheckCount = 0; %number of checks that have occured
    EventNumCheck = round(Diffusion_Rate*dt(Event));    %number of diffusion events to occur in this time step
    LeftDiffCounter = 0;
    RightDiffCounter = 0;
    MovedProteins = []; %which proteins moved in this time step
    while DiffusionEvents < EventNumCheck & CheckCount ~= numel(CurrentBound)   %completes diffusion events until the proper amount have been done
        CheckCount = CheckCount+1;
        k = DiffOrder(CheckCount);  %checks proteins for diffusion in a random order
        ProteinCheck = CurrentBound(k); %location of bound protein that's being checked right now
        if (DNA(ProteinCheck-1) == 0 & DNA(ProteinCheck+n) ~= 0) | ProteinCheck == N+2-n   %if diffusion is only possible to the left...
            if ProteinCheck == 2    %protein bound at position 2 with protein to the right cannot move
                R = 2;  %no diffusion can occur
            else
                R = rand;
            end
            if R <= (Left_Prob)  %check left diffusion probability
                DNA(ProteinCheck:ProteinCheck+(n-1)) = 0;   %clears protein from current location
                DNA(ProteinCheck-1:ProteinCheck+n-2) = 1; %protein diffuses to the left
                CurrentBound(k) = CurrentBound(k)-1;    %updates CurrentBound list
                BoundAtSpot(ProteinCheck) = 0;  %updates BoundAtSpot record
                BoundAtSpot(ProteinCheck-1) = 1;
                LeftDiffCounter = LeftDiffCounter+1;
                DiffusionEvents = DiffusionEvents+1;
                MovedProteins = [MovedProteins,ProteinCheck];   %records which proteins moved in this time step
            end
        elseif  (DNA(ProteinCheck-1) ~= 0 & DNA(ProteinCheck+n) == 0) | ProteinCheck == 2   %if diffusion is only possible to the right...
            if ProteinCheck == N+2-n    %protein bound at position N+2-n with protein to the left cannot diffuse
                R = 2;  %makes it so no diffusion can occur
            else
                R = rand;
            end
            if R <= (1-Left_Prob)
                DNA(ProteinCheck:ProteinCheck+(n-1)) = 0;   %clears protein from current location
                DNA(ProteinCheck+1:ProteinCheck+n) = 1; %protein diffuses to the right
                CurrentBound(k) = CurrentBound(k)+1;    %updates CurrentBound list
                BoundAtSpot(ProteinCheck) = 0;  %updates BoundAtSpot record
                BoundAtSpot(ProteinCheck+1) = 1;
                RightDiffCounter = RightDiffCounter+1;
                DiffusionEvents = DiffusionEvents+1;
                MovedProteins = [MovedProteins,ProteinCheck];   %records which proteins moved in this time step
            end
        elseif DNA(ProteinCheck-1) == 0 & DNA(ProteinCheck+n) == 0 & ProteinCheck ~= 2 & ProteinCheck ~= N+2-n %if diffusion is possible in either direction...
            if rand <= Left_Prob    %check diffusing to left
                DNA(ProteinCheck:ProteinCheck+(n-1)) = 0;   %clears protein from current location
                DNA(ProteinCheck-1:ProteinCheck+(n-1)-1) = 1;    %protein diffuses to the left
                CurrentBound(k) = CurrentBound(k)-1;    %updates CurrentBound list
                BoundAtSpot(ProteinCheck) = 0;  %updates BoundAtSpot record
                BoundAtSpot(ProteinCheck-1) = 1;
                LeftDiffCounter = LeftDiffCounter+1;
                DiffusionEvents = DiffusionEvents+1;
                MovedProteins = [MovedProteins,ProteinCheck];   %records which proteins moved in this time step
            else
                DNA(ProteinCheck:ProteinCheck+(n-1)) = 0;   %clears protein from current location
                DNA(ProteinCheck+1:ProteinCheck+n) = 1; %protein diffuses to the right
                CurrentBound(k) = CurrentBound(k)+1;    %updates CurrentBound list
                BoundAtSpot(ProteinCheck) = 0;  %updates BoundAtSpot record
                BoundAtSpot(ProteinCheck+1) = 1;
                RightDiffCounter = RightDiffCounter+1;
                DiffusionEvents = DiffusionEvents+1;
                MovedProteins = [MovedProteins,ProteinCheck];   %records which proteins moved in this time step
            end
        end
        if (DiffusionEvents == numel(CurrentBound)) & (DiffusionEvents  < EventNumCheck)
            MaxOut = MaxOut+1;  %counts how many times the diffusions maxed out before reaching the theoretical number of diffusions
        end
    end
    TotalDiffEvent = TotalDiffEvent+DiffusionEvents;   %counts total number of diffusion events that occur
    TotalLeftDiff = TotalLeftDiff+LeftDiffCounter;  %counts how many left diffusion events occur throughout simulation
    TotalRightDiff = TotalRightDiff+RightDiffCounter;   %counts how many right diffusion events occur throughout simulation
    DiffusedProteins(Event,1:numel(MovedProteins)) = MovedProteins;  %records which proteins diffused in which time step event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t(Event+1) = t(Event)+dt(Event);    %advances time

    ProteinCount(Event+1) = sum(DNA)/n;   %all of these values should be integers
    FracCover(Event+1) = sum(DNA)/N;    %fractional coverage of the DNA lattice
    
    ProteinTracking(Event+1,:) = DNA(2:N+1);
    
    if Event >= 100 %only tests for equilibrium after 100 events
        t_Equilibrium_Test = t(Event+1-round(0.25*(Event+1)):end);  %set of time we're looking for and testing for equilibrium
        Sat_Equilibrium_Test = FracCover(Event+1-round(0.25*(Event+1)):end);    %set of FracCover values to test for equilibrium in
        Avg_Saturation = sum(Sat_Equilibrium_Test)/numel(Sat_Equilibrium_Test); %average saturation level
        CurveFit = polyfit(t_Equilibrium_Test,Sat_Equilibrium_Test,1);  %linear fit to the last quarter of events
        Y_Int_Error(Event) = abs(Avg_Saturation-CurveFit(2))/Avg_Saturation;   %percent error in y-intercept fit
        if abs(CurveFit(1)) < 0.01 & (Y_Int_Error < 0.05 | isnan(Y_Int_Error))  %if fitted slope is basically zero and the fitted y-int is basically the average saturation...
            Equilibrium = 1;    %...then the system is considered to be at equilibirum
        else
            Equilibrium = 0;
        end
    end
end
toc
LeftDiffusionOcc = TotalLeftDiff/TotalDiffEvent; %occurence rate of a left diffusion
RightDiffusionOcc = TotalRightDiff/TotalDiffEvent;   %occurence rate of a right diffusion
Left_P_Error = (abs(LeftDiffusionOcc-Left_Prob)/Left_Prob)*100; %percent errors
Right_P_Error = (abs(RightDiffusionOcc-(1-Left_Prob))/(1-Left_Prob))*100;

disp(['Left Diffusion: ', num2str(round(LeftDiffusionOcc,2)), ' (', num2str(Left_Prob), ') (', num2str(round(Left_P_Error,3)), '% Error)']); %comparing directional diffusion probability to the actual results
disp(['Right Diffusion: ', num2str(round(RightDiffusionOcc,2)), ' (', num2str(1-Left_Prob), ') (', num2str(round(Right_P_Error,3)), '% Error)']);

figure(1);
scatter(t,FracCover,1,'r','filled');    %fractional coverage vs. dynamic time
hold on;
xlabel('Time, t');
xlim([0 max(t)]);
ylabel('Saturation Level');
ylim([0 1]);
title('ssDNA Saturation');
box on;

[X_DNA,Y_Time] = meshgrid(1:N,t);
CustMap = [1 1 1; 0 0 1];  %custom color range for white = uncovered, green = covered nt
colormap(figure(2),CustMap);

figure(2);
surf(X_DNA,Y_Time,ProteinTracking,'EdgeColor','none');
set(gca,'Ydir','reverse');  %reverses time axis so beginning is top of figure
view(2);
hold on;
xlabel('ssDNA Location');
xlim([1 N]);
ylabel('Time, t (Inverse)');
ylim([0 max(t)]);
box on;
title('Protein Diffusion');

TheoreticalDiffEvents = t(end)*Diffusion_Rate;  %theoretical number of diffusions that should have occured
DiffusionEventsError = (abs(TheoreticalDiffEvents-TotalDiffEvent)/TheoreticalDiffEvents)*100;   %percent error of diffusion events
disp(['Diffusion Events: ', num2str(round(DiffusionEventsError,2)), '% Error']);

figure(3);
scatter(1:Event,Y_Int_Error,5,'k','filled');
hold on;
xlabel('Event Number'); xlim([0 Event]);
ylabel('Y Intercept Error'); ylim([0 1]);
box on;