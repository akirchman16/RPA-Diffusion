clear;
close all;

% This code will be used to model a generic protein binding and unbinding
% to a ssDNA lattice. The protein will be able to diffuse according to D
% (diffusion rate) and L_Prob/R_Prob (probability of left/right diffusion,
% respectively). The distance traveled by each protein will be recorded and
% will be compared to traditional random walk statistics.

N = 100;    %length of DNA lattice
n = 3;  %length of protein
Iterations = 1000;

FreeProteins(1) = 50;  %number of free proteins initially

k_on = 0.1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding

D = 1e3;    %Diffusion Rate
L_Prob = 0.5;
R_Prob = 1-L_Prob;

DNA = zeros(1,N);   %ssDNA lattice
FracCover = zeros(1,Iterations);
t = zeros(1,Iterations);
BoundAtSpot = zeros(1,N);   %tracks where proteins are currently bound
TotalDiffusionEvents = 0;
LeftDiffEvents = 0;
RightDiffEvents = 0;
for i = 1:Iterations
    % Lattice Search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LeftOpen = find(diff([1 DNA 1]) == -1);   %left edge of open spaces on lattice
    RightOpen = find(diff([1 DNA 1]) == 1)-1;   %right positions of each open space
    SpaceSizes = RightOpen-LeftOpen+1;  %size of each open space
    
    L_Available4Binding = LeftOpen(SpaceSizes >= n);  %left positions of open spaces that are large enough for binding
    SpacesFree = SpaceSizes(SpaceSizes >= n);
    Free4Binding = [];
    for a = 1:numel(L_Available4Binding)
        Free4Binding = [Free4Binding,L_Available4Binding(a):1:(L_Available4Binding(a)+(SpacesFree(a)-n))];
    end
    X_Free(i) = numel(Free4Binding);
    X_Bound(i) = numel(find(BoundAtSpot == 1));
    
    a = [k_on*X_Free(i)*FreeProteins(i),k_off*X_Bound(i)];  %propensity functions
    a_0 = sum(a);
    
    if a(1) >= rand*a_0    %binding reaction
        j = 1;
        BindSpot = Free4Binding(randi(numel(Free4Binding)));    %random spot for binding to occur
        DNA(BindSpot:BindSpot+(n-1)) = 1;   %bind protein to the lattice
        BoundAtSpot(BindSpot) = 1;
        FreeProteins(i+1) = FreeProteins(i)-1;
    else   %unbinding reaction
        j = 2;
        BoundProteins = find(BoundAtSpot == 1); %all locations where proteins are bound
        UnbindSpot = BoundProteins(randi(numel(BoundProteins)));    %random spot for unbinding to occur
        DNA(UnbindSpot:UnbindSpot+(n-1)) = 0;   %unbind protein from lattice
        BoundAtSpot(UnbindSpot) = 0;
        FreeProteins(i+1) = FreeProteins(i)+1;
    end
    
    dt(i) = (1/a_0)*log(1/rand);   %time step
    t(i+1) = t(i)+dt(i);
    
    FracCover(i+1) = sum(DNA)/N;  %saturation level
    
    DiffusionEventsCheck = round(dt(i)*D); %number of Diffusion events that need to occur
    DiffusionCounter = 0;
    while (DiffusionCounter < DiffusionEventsCheck) && (ismember(1,abs(diff(DNA)))) && (sum(DNA) ~= 0)
        BoundLocations = find(BoundAtSpot == 1); %list of all locations where proteins are bound
        DiffProteinCheck = BoundLocations(randi(numel(BoundLocations)));   %random protein to test for diffusing
        if DiffProteinCheck == 1    %left-most position chosen
            if DNA(DiffProteinCheck+n) == 0 %diffusion to right is possible
                DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
                DNA(DiffProteinCheck+1:DiffProteinCheck+n) = 1; %right diffusion
                BoundAtSpot(1,DiffProteinCheck:DiffProteinCheck+1) = [0,1]; %updates BoundAtSpot
                RightDiffEvents = RightDiffEvents+1;
                DiffusionCounter = DiffusionCounter+1;
                TotalDiffusionEvents = TotalDiffusionEvents+1;
            end
        elseif DiffProteinCheck == N-(n-1)  %right-most position chosen
            if DNA(DiffProteinCheck-1) == 0 %diffusion to the left is possible
                DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
                DNA(DiffProteinCheck-1:DiffProteinCheck+(n-1)-1) = 1; %left diffusion
                BoundAtSpot(1,DiffProteinCheck-1:DiffProteinCheck) = [1,0]; %updates BoundAtSpot
                LeftDiffEvents = LeftDiffEvents+1;
                DiffusionCounter = DiffusionCounter+1;
                TotalDiffusionEvents = TotalDiffusionEvents+1;
            end
        elseif DNA(DiffProteinCheck-1) ~= 0 && DNA(DiffProteinCheck+n) == 0    %right diffusion only
            DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
            DNA(DiffProteinCheck+1:DiffProteinCheck+n) = 1; %right diffusion
            BoundAtSpot(1,DiffProteinCheck:DiffProteinCheck+1) = [0,1]; %updates BoundAtSpot
            RightDiffEvents = RightDiffEvents+1;
            DiffusionCounter = DiffusionCounter+1;
            TotalDiffusionEvents = TotalDiffusionEvents+1;
        elseif DNA(DiffProteinCheck-1) == 0 && DNA(DiffProteinCheck+n) ~= 0    %left diffusion only
            DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
            DNA(DiffProteinCheck-1:DiffProteinCheck+(n-1)-1) = 1; %left diffusion
            BoundAtSpot(1,DiffProteinCheck-1:DiffProteinCheck) = [1,0]; %updates BoundAtSpot
            LeftDiffEvents = LeftDiffEvents+1;
            DiffusionCounter = DiffusionCounter+1;
            TotalDiffusionEvents = TotalDiffusionEvents+1;
        else   %either direction is possible
            if rand <= L_Prob    %check left probability
                DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
                DNA(DiffProteinCheck-1:DiffProteinCheck+(n-1)-1) = 1; %left diffusion
                BoundAtSpot(1,DiffProteinCheck-1:DiffProteinCheck) = [1,0]; %updates BoundAtSpot
                LeftDiffEvents = LeftDiffEvents+1;
                DiffusionCounter = DiffusionCounter+1;
                TotalDiffusionEvents = TotalDiffusionEvents+1;
            else
                DNA(DiffProteinCheck:DiffProteinCheck+(n-1)) = 0;
                DNA(DiffProteinCheck+1:DiffProteinCheck+n) = 1; %right diffusion
                BoundAtSpot(1,DiffProteinCheck:DiffProteinCheck+1) = [0,1]; %updates BoundAtSpot
                RightDiffEvents = RightDiffEvents+1;
                DiffusionCounter = DiffusionCounter+1;
                TotalDiffusionEvents = TotalDiffusionEvents+1;
            end
        end
    end
end

figure();
scatter(t,FracCover,3,'r','filled');
xlabel('Time, t');  ylabel('Saturation Level');
box on; xlim([0 max(t)]);   ylim([0 1]);
title('Saturation of ssDNA (Single Protein)');