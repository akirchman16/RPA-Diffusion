clearvars;
close all;

%This code will simply describe a general protein diffusing around on a
%1-dimensional lattice. It will run a SSA calculation for binding and
%unbinding, but with each step a simple Monte-Carlo method will be used to
%test for diffusion in one direction or another if possible. The diffusion
%will essentially be a random walk simulation for the proteins

N = 1000;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 1;      %cooperativity parameter
L_A = 1;  %concentration of free proteins

Diffusion_Prob = 0.5;   %probability of any protein diffusing on lattice
Left_Prob = 0.5;    %probability of a protein diffusing right

Iterations = 1000;  %number of events which will occur

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

t(1) = 0;
BindCounter = 0;
UnbindCounter = 0;

xB(1) = N-(n-1);    %initial values for a free lattice
xAB(1) = 0;

ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice


for i = 1:Iterations
    a_f(i) = k_on*(L_A)*(xB(i));   %propensity functions (probability of each event happening)
    a_r(i) = k_off*(xAB(i));
    a_0(i) = a_f(i)+a_r(i);     %sum of propensity functions used for determining dt

    dt(i) = (1/a_0(i))*log(1/rand); %random time interval for Gillespie method
    
    for j = 2:N-(n-1)+1   %loop to calculate probability of binding at each spot if binding event occurs
        if (DNA(j) ~= 1) && (sum(DNA(j:j+(n-1))) == 0)
            if DNA(j-1) == 0 && DNA(j+n) == 0   %checks for isolated location
                BindProb(i,j-1) = 1/(xB(i));
            elseif (DNA(j-1) == 0 && DNA(j+n) == 1) || (DNA(j-1) == 1 && DNA(j+n) == 0) %checks for singly contiguous location
                BindProb(i,j-1) = w/(xB(i));
            elseif DNA(j-1) == 1 && DNA(j+n) == 1 %checks for doubly contiguous location
                BindProb(i,j-1) = (w^2)/(xB(i));
            end
        end
    end

    if a_f(i) > rand*a_0(i) %a forward reaction occurs
        eventB = 0;
        while ~eventB   %repeats until a binding occurs
            SpotB = randsample(Locations,1,true,BindProb(i,1:N-n+1))+1;  %random location on lattice is chosen
            if DNA(SpotB:SpotB+(n-1)) == 0   %checks if location is free
               DNA(SpotB:SpotB+(n-1)) = 1;   %binds protein to location
               BoundAtSpot(SpotB) = 1;  %stores locatin in BoundAtSpot
               BindCounter = BindCounter+1;    %updates bind counter
               Hist(1,i) = SpotB;      %stores location of binding event

               FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %function to find number of free spots now available
               xB(i+1) = FreeSpots;    %updates populations
               xAB(i+1) = xAB(i)+1;

               eventB = 1;
            end
        end
        if DNA(SpotB-1) == 0 && DNA(SpotB+n) == 0 %isolated event occured
            CoopEvents(1,i) = SpotB;
        elseif (DNA(SpotB-1) == 0 && DNA(SpotB+n) == 1) || (DNA(SpotB-1) == 1 && DNA(SpotB+n) == 0) %singly contiguous event occured
            CoopEvents(2,i) = SpotB;
        elseif DNA(j-1) == 1 && DNA(j+n) == 1   %doubly contiguous binding event occured
            CoopEvents(3,i) = SpotB;
        end
    else                       %otherwise an unbinding has to occur
        CurrentBound = find(BoundAtSpot == 1);
        pos = randi(length(CurrentBound));
        SpotU = CurrentBound(pos);  %random position for protein to unbind

        DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
        BoundAtSpot(SpotU) = 0; %removes location from BoundAtSpot
        UnbindCounter = UnbindCounter+1;    %updates unbind counter
        Hist(2,i) = SpotU;  %stores location of unbinding event

        FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %calculates number of free spaces
        xB(i+1) = FreeSpots;
        xAB(i+1) = xAB(i)-1;
    end
    CurrentBound = find(BoundAtSpot == 1);
    for k = 1:numel(find(BoundAtSpot == 1)) %check each bound protein for diffusion possibilities
        BoundProtein = CurrentBound(k);
        if (DNA(BoundProtein-1) == 0 & DNA(BoundProtein+n:BoundProtein+n) ~= 0) | BoundProtein == N+2-n    %if diffusion is only possible to the left...
            if BoundProtein == 2    %protein bound at position 2 with protein to the right cannot move
                R = 1;  %makes it so no diffusion can occur
            else
                R = rand;
            end
            if R <= Diffusion_Prob   %check diffusion probability
                DNA(BoundProtein:BoundProtein+(n-1)) = 0;   %clears protein from current location
                DNA(BoundProtein-1:BoundProtein+(n-1)-1) = 1;   %protein diffuses to the left
                CurrentBound(k) = CurrentBound(k)-1;    %updates CurrentBound list
                BoundAtSpot(BoundProtein) = 0;  %updates BoundAtSpot record
                BoundAtSpot(BoundProtein-1) = 1;
            end
        elseif (DNA(BoundProtein-1) ~= 0 & DNA(BoundProtein+n:BoundProtein+n) == 0) | BoundProtein == 2   %if diffusion is only possible to the right...
            if BoundProtein == N+2-n    %protein bound at position N+2-n with protein to the left cannot move
                R = 1;  %makes it so no diffusion can occur
            else
                R = rand;
            end
            if rand <= Diffusion_Prob   %check diffusion probability
                DNA(BoundProtein:BoundProtein+(n-1)) = 0;   %clears protein from current location
                DNA(BoundProtein+1:BoundProtein+n) = 1; %protein diffuses to the right
                CurrentBound(k) = CurrentBound(k)+1;    %updates CurrentBound list
                BoundAtSpot(BoundProtein) = 0;  %updates BoundAtSpot record
                BoundAtSpot(BoundProtein+1) = 1;
            end
        elseif DNA(BoundProtein-1) == 0 & DNA(BoundProtein+n:BoundProtein+n) == 0 & BoundProtein ~= 2 & BoundProtein ~= N+2-n   %if diffusion is possible in both directions...
            if rand <= Diffusion_Prob  %check diffusion probability
                if rand <= Left_Prob    %check diffusing to left
                    DNA(BoundProtein:BoundProtein+(n-1)) = 0;   %clears protein from current location
                    DNA(BoundProtein-1:BoundProtein+(n-1)-1) = 1;    %protein diffuses to the left
                    CurrentBound(k) = CurrentBound(k)-1;    %updates CurrentBound list
                    BoundAtSpot(BoundProtein) = 0;  %updates BoundAtSpot record
                    BoundAtSpot(BoundProtein-1) = 1;
                else
                    DNA(BoundProtein:BoundProtein+(n-1)) = 0;   %clears protein from current location
                    DNA(BoundProtein+1:BoundProtein+n) = 1; %protein diffuses to the right
                    CurrentBound(k) = CurrentBound(k)+1;    %updates CurrentBound list
                    BoundAtSpot(BoundProtein) = 0;  %updates BoundAtSpot record
                    BoundAtSpot(BoundProtein+1) = 1;
                end
            end
        end
    end
    
    t(i+1) = t(i)+dt(i);    %advances time

    ProteinCount(i+1) = sum(DNA)/n;   %all of these values should be integers
    FracCover(i+1) = sum(DNA)/N;    %fractional coverage of the DNA lattice
end

figure();
scatter(t,FracCover,1,'r','filled');    %fractional coverage vs. dynamic time
xlabel('Time, t');
xlim([0 max(t)]);
ylabel('Saturation Level');
ylim([0 1]);
title('ssDNA Saturation');
box on;