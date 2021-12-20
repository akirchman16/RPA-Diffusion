clearvars;
close all;

% This code will be a competition model between RPA and RAD51. It will run
% in the same manner as the codes which generated heatmaps but without
% multiple runs. It will only test a single set of paramters. RPA molecules
% will have the possibility to diffuse along the lattice if they're able
% to. This will be set by the RPA_DiffusionProb.

N = 5000;   %ssDNA length
DNA = zeros(2,N);   %represents ssDNA lattice (2nd row is real lattice)

minIterations = 10000;

%RAD51 Properties/Parameters
RAD51 = 51; %how RAD51 will be represented on the lattice
n_RAD51 = 3;    %size of RAD51 protein
TotalCount_RAD51 = 1000; %number of total RAD51 proteins (monomers)
w_RAD51 = 1;    %cooperativity constant for RAD51
k_on_RAD51 = 1; %kinetic rate constant for RAD51 binding
k_off_RAD51 = 1;    %kinetic rate constant for RAD51 unbinding

%RPA Properties/Parameters
RPA_A = 1; %represent RPA-A on lattice
RPA_D = 3;  %represent RPA-D on lattice
n_A = 10;   %size of RPA-A
n_D = 10;   %size of RPA-D
n_RPA = n_A+n_D;
TotalCount_RPA = 125;    %total number of RPA proteins that exist
w_RPA = 1;  %cooperativity of RPA (IDK if this is fully included in the model currently)
k_on_RPA_A = 100;    %kinetic rate constant for RPA-A binding
k_off_RPA_A = 1;    %kinetic rate constant for RPA-A unbinding
k_on_RPA_D = 50;    %kinetic rate consant for RPA-D binding
k_off_RPA_D = 1;    %kinetic rate constant for RPA-D unbinding

DiffusionRate = 100;    %RPA Diffusion Rate constant (events/time interval)

%Memory Allocation
t = zeros(1,minIterations);
xRAD51_M = [TotalCount_RAD51,zeros(1,minIterations-1)];
xRAD51_D = zeros(1,minIterations);
xRPA = [TotalCount_RPA,zeros(1,minIterations-1)];
Free_Proteins = zeros(3,minIterations);
a_Prop = zeros(15,minIterations);
a_0 = zeros(1,minIterations);
dt = zeros(1,minIterations);
LocHist = zeros(15,minIterations);
FracCover_RAD51 = zeros(1,minIterations);
FracCover_RPA_A = zeros(1,minIterations);
FracCover_RPA_D = zeros(1,minIterations);
FracCover_RPA = zeros(1,minIterations);
FracCover_Total = zeros(1,minIterations);
RAD51_Mon_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Monomers are bound
RAD51_Dim_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Dimers are bound
RPA_A_BoundAtSpot = zeros(1,N); %array to record where RPA-A is actively bound
RPA_D_BoundAtSpot = zeros(1,N); %array to record where RPA-D is actively bound
RPA_D_HingedOpen = zeros(1,N);  %array to record where RPA-D is microscopically dissociated from lattice
RPA_A_HingedOpen = zeros(1,N);  %array to record where RPA_A is microscopically dissociated from the lattice
j = zeros(1,minIterations);

Free_Proteins(:,1) = [xRAD51_M(1) ; xRAD51_D(1) ; xRPA(1)];    %starts matrix to count free protein populations

Equilibrium_RPA = 0;
Equilibrium_RAD51 = 0;
Equilibrium = double(Equilibrium_RPA && Equilibrium_RAD51);
Event = 0;
% while Equilibrium ~= 1
for i = 1:minIterations
    Event = Event+1;
    
%%% Lattice Search Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gap_Left = find(diff([1 DNA(2,:) 1])<0 & diff([1 DNA(2,:) 1]) ~= RPA_A-RPA_D & diff([1 DNA(2,:) 1]) ~=RPA_D-RAD51 & diff([1 DNA(2,:) 1]) ~= RPA_A-RAD51);    %left most available location of all gaps on lattice
    Gap_Right = find(diff([1 DNA(2,:) 1])>0 & diff([1 DNA(2,:) 1]) ~= RPA_D-RPA_A & diff([1 DNA(2,:) 1]) ~= RAD51-RPA_D & diff([1 DNA(2,:) 1]) ~= RAD51-RPA_A)-1; %right most available location of all gaps on lattice
    Gap_Size = (Gap_Right-Gap_Left)+1;    %calculates the size of each gap
    Gap_Edges = [Gap_Left;Gap_Right];   %lists all gap edges (1st row = left edge, 2nd row = right edge)
    
    RAD51_M_Available_Gap_Edges = [Gap_Left(Gap_Size >= n_RAD51);Gap_Right(Gap_Size >= n_RAD51)];   %lists the left (1st row) and right (2nd row) edges
                                                                                                    %of gaps large enough for RAD51 monomers
    RAD51_D_Available_Gap_Edges = [Gap_Left(Gap_Size >= 2*n_RAD51);Gap_Right(Gap_Size >= 2*n_RAD51)];   %lists the left (1st row) and right (2nd row) edges
                                                                                                        %of gaps large enough for RAD51 dimers
    RPA_Available_Gap_Edges = [Gap_Left(Gap_Size >= n_RPA);Gap_Right(Gap_Size >= n_RPA)];   %lists left (1st row) and right (2nd row) edges of gaps
                                                                                            %large enough for an RPA protein to bind
    %Doubly Contiguous Search - Easiest
    RAD51_Mon_DC = Gap_Left(Gap_Size == n_RAD51 & Gap_Left > 1 & Gap_Left < N-(n_RAD51-1));   %all available RAD51 Mon. DC sites
    RAD51_Dim_DC = Gap_Left(Gap_Size == 2*n_RAD51 & Gap_Left > 1 & Gap_Left < N-(2*n_RAD51-1)); %all available RAD51 Dim. DC sites
    RPA_DC = Gap_Left(Gap_Size == n_RPA & Gap_Left > 1 & Gap_Left < N-(n_RPA-1));   %all available RPA DC sites
    
    %Singly Contiguous Search
    RAD51_Mon_SC = unique([RAD51_M_Available_Gap_Edges(1,:),RAD51_M_Available_Gap_Edges(2,:)-(n_RAD51-1)]); %position of all locations at the edges of gaps available for RAD51 Mon.
    RAD51_Mon_SC(ismember(RAD51_Mon_SC,RAD51_Mon_DC)) = [];     %clears location if it's already been recorded as a DC location
    RAD51_Mon_SC(ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == RPA_A)-n_RAD51) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == RPA_D)-n_RAD51) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == -RPA_A)) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == -RPA_D)) == 1) = []; %clears any position that neighbors a protein other than RAD51 
    if DNA(2,n_RAD51+1) ~= RAD51    %if the first location isn't SC...
        RAD51_Mon_SC(RAD51_Mon_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-n_RAD51) ~= RAD51       %if the final available location on lattice isn't SC...
        RAD51_Mon_SC(RAD51_Mon_SC == N-(n_RAD51-1)) = [];   %...clear it
    end
    RAD51_Dim_SC = unique([RAD51_D_Available_Gap_Edges(1,:),RAD51_D_Available_Gap_Edges(2,:)-(2*n_RAD51-1)]);   %positions of all locations at the edges of gaps available for RAD51 Dim.
    RAD51_Dim_SC(ismember(RAD51_Dim_SC,RAD51_Dim_DC)) = [];     %clears locations that have already been counted as a DC location
    RAD51_Dim_SC(ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == RPA_A)-(2*n_RAD51)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == RPA_D)-(2*n_RAD51)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == -RPA_A)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == -RPA_D)) == 1) = []; %clears all position that neighbor a protein other than RAD51
    if DNA(2,2*n_RAD51+1) ~= RAD51    %if the first location isn't SC...
        RAD51_Dim_SC(RAD51_Dim_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-(2*n_RAD51)) ~= RAD51       %if the final available location on lattice isn't SC...
        RAD51_Dim_SC(RAD51_Dim_SC == N-(2*n_RAD51-1)) = [];   %...clear it
    end
    RPA_SC = unique([RPA_Available_Gap_Edges(1,:),RPA_Available_Gap_Edges(2,:)-(n_RPA-1)]);   %positions of all locations at the edges of gaps available for RPA
    RPA_SC(ismember(RPA_SC,RPA_DC)) = [];     %clears locations that have already been counted as a DC location
    RPA_SC(ismember(RPA_SC,find(diff([0 DNA(2,:) 0]) == RAD51)-n_RPA) == 1 | ismember(RPA_SC,find(diff([0 DNA(2,:) 0]) == -RAD51)) == 1) = []; %clears all locations that don't neighbor an RPA protein (A or D)
    if DNA(2,n_RPA+1) ~= RPA_A   %if the first location isn't SC...
        RPA_SC(RPA_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-n_RPA) ~= RPA_D & DNA(2,N-n_RPA) ~= RPA_A     %if the final available location on lattice isn't SC...
        RPA_SC(RPA_SC == N-(n_RPA-1)) = [];   %...clear it
    end

    %Isolated Search
    RAD51_Mon_I = [];   %initializes an empty array for RAD51 monomer isolated available sites
    Gap_Size_RAD51_M_I = Gap_Size(Gap_Size > n_RAD51+1);    %gap sizes of gaps large enough for a RAD51 monomer to bind in an isolated position
    RAD51_M_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RAD51_M_I) == 1);  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for a = 1:length(RAD51_M_Available_I_Gap_Edges(1,:))
        RAD51_Mon_I = [RAD51_Mon_I, RAD51_M_Available_I_Gap_Edges(1,a)+1:1:RAD51_M_Available_I_Gap_Edges(2,a)-n_RAD51]; %adds sequential locations that are isolated locations available for RAD51 Mon
    end
    if DNA(2,n_RAD51) ~= RAD51 & DNA(2,1:n_RAD51) == 0   %if left end of lattice is available and isn't SC...
        RAD51_Mon_I = [1,RAD51_Mon_I];  %...add the first location to the list
    end
    if DNA(2,N-n_RAD51) ~= RAD51 & DNA(2,N-(n_RAD51-1):end) == 0   %if the right end of the lattice is available and isn't SC...
        RAD51_Mon_I = [RAD51_Mon_I,N-(n_RAD51-1)];  %...add the last location to the list
    end
    RAD51_Dim_I = [];   %initializes an empty array for RAD51 dimer isolated available sites
    Gap_Size_RAD51_D_I = Gap_Size(Gap_Size > 2*n_RAD51+1);    %gap sizes large enough for RAD51 Dimers
    RAD51_D_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RAD51_D_I));  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for b = 1:length(RAD51_D_Available_I_Gap_Edges(1,:))
        RAD51_Dim_I = [RAD51_Dim_I, RAD51_D_Available_I_Gap_Edges(1,b)+1:1:RAD51_D_Available_I_Gap_Edges(2,b)-(2*n_RAD51)]; %adds sequential locations that are isolated locations available for RAD51 Dim
    end
    if DNA(2,1:2*n_RAD51) == 0 & DNA(2,(2*n_RAD51)+1) ~= RAD51   %if left end of lattice is available and isn't SC...
        RAD51_Dim_I = [1,RAD51_Dim_I];  %...add the first location to the list
    end
    if DNA(2,N-(2*n_RAD51-1):end) == 0 & DNA(2,N-(2*n_RAD51)) ~= RAD51   %if the right end of the lattice isn't SC...
        RAD51_Dim_I = [RAD51_Dim_I,N-(2*n_RAD51-1)];  %...add the last location to the list
    end
    RPA_I = [];   %initializes an empty array for RPA isolated available sites
    Gap_Size_RPA = Gap_Size(Gap_Size > n_RPA);  %gap sizes which are large enough for RPA isolated binding
    RPA_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RPA));  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for c = 1:length(RPA_Available_I_Gap_Edges(1,:))
        RPA_I = [RPA_I, RPA_Available_I_Gap_Edges(1,c)+1:1:RPA_Available_I_Gap_Edges(2,c)-n_RPA]; %adds sequential locations that are isolated locations available for RAD51 Mon
    end
    if DNA(2,1:n_RPA) == 0 & DNA(2,n_RPA+1) ~= RPA_A   %if left end of lattice is available and isn't SC...
        RPA_I = [1,RPA_I];  %...add the first location to the list
    end
    if DNA(2,N-(n_RPA-1):end) == 0 & DNA(2,N-(n_RPA+1)) ~= RPA_A & DNA(2,N-(n_RPA+1)) ~= RPA_D  %if the right end of the lattice is available and isn't SC...
        RPA_I = [RPA_I,N-(n_RPA-1)];  %...add the last location to the list
    end

    %RPA-D Availble Locations
    Free_for_RPA_D = [];    %initializes an array for where hinged open RPA-D can rebind to DNA
    for d = find(RPA_D_HingedOpen == 1)   %check each location where RPA-D is hinged open
        if DNA(2,d:d+(n_D-1)) == 0     %if no proteins are bound below hinged open RPA-D...
            Free_for_RPA_D = [Free_for_RPA_D,d];    %...adds location to store where RPA-D can rebind this event
        end
    end
    Free_for_RPA_A = [];    %initializes an array for where hinged open RPA-A can rebind to DNA
    for f = find(RPA_A_HingedOpen == 1)   %check each location where RPA-A is hinged open
        if DNA(2,f:f+(n_A-1)) == 0     %if no proteins are bound below hinged open RPA-A...
            Free_for_RPA_A = [Free_for_RPA_A,f];    %...adds location to store where RPA-A can rebind this event
        end
    end
    Available_HingeClosed = [numel(Free_for_RPA_A),numel(Free_for_RPA_D)];  %stores number of available locations where RPA can microscopically rebind
    %Population of Bound Proteins
    x_Bound_RAD51_M = numel(find(DNA(2,:) == RAD51))/n_RAD51;    %calculates how many RAD51 Monomers are actively bound to DNA
    x_Bound_RPA_A = numel(find(DNA(2,:) == RPA_A))/n_A;    %calculates the number of bound RPA-A
    x_Bound_RPA_D = numel(find(DNA(2,:) == RPA_D))/n_D;    %calculates the number of actively bound RPA-D
    %Search for RAD51 Dimers (Uses Filament Lengths)
    RAD51_Dim_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Dimers are bound on the lattice
    RAD51_Filament_Edges = [find(diff([0 DNA(2,:) 0]) == 51 | diff([0 DNA(2,:) 0]) == 50 | diff([0 DNA(2,:) 0]) == 48);find(diff([0 DNA(2,:) 0]) == -51 | diff([0 DNA(2,:) 0]) == -50 | diff([0 DNA(2,:) 0]) == -48)-1];  %edges of RAD51 filaments (left = 1st row; right = 2nd row)
    RAD51_Filament_Lengths = RAD51_Filament_Edges(2,:)-RAD51_Filament_Edges(1,:)+1; %length of each RAD51 Filament on the lattice
    RAD51_D_Filament_Locations = RAD51_Filament_Edges(1,RAD51_Filament_Lengths > n_RAD51);  %location of all filaments which contain a dimer
    RAD51_D_Filament_Lengths = RAD51_Filament_Lengths(RAD51_Filament_Lengths > n_RAD51);   %length of RAD51 filaments which contain a dimer
    Monomers_per_Dimer_Filament = RAD51_D_Filament_Lengths./n_RAD51;    %number of monomers in each filament that contains a dimer
    Left_RAD51_Dimer_Filament = []; %initializes array to record locations of dimers in filaments
    if numel(RAD51_D_Filament_Locations) > 0   %if there are any filaments containing dimers...
        for e = 1:numel(RAD51_D_Filament_Locations)     %...check each one to add it to the list of dimer locations
            if RAD51_D_Filament_Lengths(e) == 2*n_RAD51    %if there's only one dimer in the filament...
                Left_RAD51_Dimer_Filament = [Left_RAD51_Dimer_Filament,RAD51_D_Filament_Locations(e)];  %...then record the location of that dimer   %...then record that there is a dimer bound at corresponding location
            else   %...otherwise we have to check each monomer location
                for f = 1:Monomers_per_Dimer_Filament(e)-1     %record that there's a possible dimer at each monomer location (except the last one)
                    Left_RAD51_Dimer_Filament = [Left_RAD51_Dimer_Filament,RAD51_D_Filament_Locations(e)+((f-1)*n_RAD51)];  %record location of the beginning of each dimer within the filament
                end
            end
        end
    end
    RAD51_Dim_BoundAtSpot(Left_RAD51_Dimer_Filament) = 1;   %records where all possible dimers are located
    x_Bound_RAD51_D = numel(find(RAD51_Dim_BoundAtSpot == 1)); %number of RAD51 Dimers bound to lattice
    
    Available_Counts = [numel(RAD51_Mon_I),numel(RAD51_Mon_SC),numel(RAD51_Mon_DC),numel(RAD51_Dim_I),numel(RAD51_Dim_SC),numel(RAD51_Dim_DC),numel(RPA_I),numel(RPA_SC),numel(RPA_DC)];
                %populations of the available locations listed in order
                % 1 - RAD51 Monomer Isolated
                % 2 - RAD51 Monomer Singly Contiguous
                % 3 - RAD51 Monomer Doubly Contiguous
                % 4 - RAD51 Dimer Isolated
                % 5 - RAD51 Dimer Singly Contiguous
                % 6 - RAD51 Dimer Doubly Contiguous
                % 7 - RPA Isolated
                % 8 - RPA Singly Contiguous
                % 9 - RPA Doubly Contiguous
     Bound_Counts = [x_Bound_RAD51_M,x_Bound_RAD51_D,x_Bound_RPA_A,x_Bound_RPA_D];  %counts of how many bound proteins of each type
                % 1 - RAD51 Monomer
                % 2 - RAD51 Dimer
                % 3 - RPA-A
                % 4 - RPA-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_RAD51_Mon = [k_on_RAD51*Available_Counts(1)*Free_Proteins(1,Event);k_on_RAD51*Available_Counts(2)*Free_Proteins(1,Event)*w_RAD51;k_on_RAD51*Available_Counts(3)*Free_Proteins(1,Event)*(w_RAD51^2)];  %propensity functions for RAD51 monomer binding
    a_RAD51_Dim = [k_on_RAD51*Available_Counts(4)*Free_Proteins(2,Event);k_on_RAD51*Available_Counts(5)*Free_Proteins(2,Event)*w_RAD51;k_on_RAD51*Available_Counts(6)*Free_Proteins(2,Event)*(w_RAD51^2)];  %propensity functions for RAD51 Dimer binding
    a_RPA_Macro = [k_on_RPA_A*Available_Counts(7)*Free_Proteins(3,Event);k_on_RPA_A*Available_Counts(8)*Free_Proteins(3,Event)*w_RPA;k_on_RPA_A*Available_Counts(9)*Free_Proteins(3,Event)*(w_RPA^2)];  %propensity functions for RPA Macro binding
    a_RAD51_Unbind = [k_off_RAD51*Bound_Counts(1);k_off_RAD51*Bound_Counts(2)]; %propensity functions for RAD51 unbinding
    a_RPA_Micro = [k_on_RPA_A*Available_HingeClosed(1)*numel(find(RPA_A_HingedOpen == 1));k_on_RPA_D*Available_HingeClosed(2)*numel(find(RPA_D_HingedOpen == 1));k_off_RPA_A*Bound_Counts(3);k_off_RPA_D*Bound_Counts(4)];    %propensity functions for RPA Microscopic binding and unbinding
    a_Prop(:,Event) = [a_RAD51_Mon;a_RAD51_Dim;a_RPA_Macro;a_RAD51_Unbind;a_RPA_Micro];  %all propensity functions together
    a_0(Event) = sum(a_Prop(:,Event));  %sum of propensity functions (for event selection)

    Randoms = [rand,rand];    %random numbers for Monte Carlo steps
    dt(Event) = (1/a_0(Event))*log(1/Randoms(1)); %time until next reaction occurs
    
    if a_Prop(1,Event) >= Randoms(2)*a_0(Event)             %RAD51 Monomer Isolated Binding
        j(Event) = 1;
        RAD51_Mon_I_Bind_Loc = RAD51_Mon_I(randi(numel(RAD51_Mon_I)));  %random location for binding
        DNA(2,RAD51_Mon_I_Bind_Loc:RAD51_Mon_I_Bind_Loc+(n_RAD51-1)) = RAD51; %binds protein to the lattice
        RAD51_Mon_BoundAtSpot(RAD51_Mon_I_Bind_Loc) = 1; %records where protein is now bound to
        LocHist(1,Event) = RAD51_Mon_I_Bind_Loc;    %records where each event occured
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[-1;0;0];  %updates free protein counter
    elseif sum(a_Prop(1:2,Event)) >= Randoms(2)*a_0(Event)  %RAD51 Monomer Singly Contiguous Binding
        j(Event) = 2;
        RAD51_Mon_SC_Bind_Loc = RAD51_Mon_SC(randi(numel(RAD51_Mon_SC)));   %random SC location for RAD51 binding
        DNA(2,RAD51_Mon_SC_Bind_Loc:RAD51_Mon_SC_Bind_Loc+(n_RAD51-1)) = RAD51; %binds protein to lattice
        RAD51_Mon_BoundAtSpot(RAD51_Mon_SC_Bind_Loc) = 1;   %records where a protein is currently bound
        LocHist(2,Event) = RAD51_Mon_I_Bind_Loc;    %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[-1;0;0];  %updates free protein counter
    elseif sum(a_Prop(1:3,Event)) >= Randoms(2)*a_0(Event)  %RAD51 Monomer Doubly Contiguous Binding
        j(Event) = 3;
        RAD51_Mon_DC_Bind_Loc = RAD51_Mon_DC(randi(numel(RAD51_Mon_DC)));   %random DC location for RAD51 binding
        DNA(2,RAD51_Mon_DC_Bind_Loc:RAD51_Mon_DC_Bind_Loc+(n_RAD51-1)) = RAD51; %bind protein to lattice
        RAD51_Mon_BoundAtSpot(RAD51_Mon_DC_Bind_Loc) = 1;   %where protein is currently bound
        LocHist(3,Event) = RAD51_Mon_DC_Bind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[-1;0;0];  %updates free protein counter
    elseif sum(a_Prop(1:4,Event)) >= Randoms(2)*a_0(Event)  %RAD51 Dimer Isolated Binding
        j(Event) = 4;
        RAD51_Dim_I_Bind_Loc = RAD51_Dim_I(randi(numel(RAD51_Dim_I)));  %random I location for RAD51 Dimer binding
        DNA(2,RAD51_Dim_I_Bind_Loc:RAD51_Dim_I_Bind_Loc+(2*n_RAD51-1)) = RAD51; %binds RAD51 protein
        RAD51_Mon_BoundAtSpot([RAD51_Dim_I_Bind_Loc,RAD51_Dim_I_Bind_Loc+n_RAD51]) = 1; %records where each monomer is bound
        LocHist(4,Event) = RAD51_Dim_I_Bind_Loc;    %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;-1;0];  %updates free protein counter
    elseif sum(a_Prop(1:5,Event)) >= Randoms(2)*a_0(Event)  %RAD51 Dimer Singly Contiguous Binding
        j(Event) = 5;
        RAD51_Dim_SC_Bind_Loc = RAD51_Dim_SC(randi(numel(RAD51_Dim_SC)));   %random SC location for RAD51 Dimer binding
        DNA(2,RAD51_Dim_SC_Bind_Loc:RAD51_Dim_SC_Bind_Loc+(2*n_RAD51-1)) = RAD51;    %binds RAD51 dimer
        RAD51_Mon_BoundAtSpot([RAD51_Dim_SC_Bind_Loc,RAD51_Dim_SC_Bind_Loc+n_RAD51]) = 1;   %records where each monomer is bound
        LocHist(5,Event) = RAD51_Dim_SC_Bind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;-1;0];  %updates free protein counter
    elseif sum(a_Prop(1:6,Event)) >= Randoms(2)*a_0(Event)  %RAD51 Dimer Doubly Contiguous Binding
        j(Event) = 6;
        RAD51_Dim_DC_Bind_Loc = RAD51_Dim_DC(randi(numel(RAD51_Dim_DC)));   %random DC location for RAD51 Dimer binding
        DNA(2,RAD51_Dim_DC_Bind_Loc:RAD51_Dim_DC_Bind_Loc+(2*n_RAD51-1)) = RAD51;    %binds RAD51 dimer
        RAD51_Mon_BoundAtSpot([RAD51_Dim_DC_Bind_Loc,RAD51_Dim_DC_Bind_Loc+n_RAD51]) = 1;   %records where each monomer is bound
        LocHist(6,Event) = RAD51_Dim_DC_Bind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;-1;0];  %updates free protein counter
    elseif sum(a_Prop(1:7,Event)) >= Randoms(2)*a_0(Event)  %RPA Macro Isolated Binding
        j(Event) = 7;
        RPA_Macro_I_Bind_Loc = RPA_I(randi(numel(RPA_I)));   %random I location for RPA Marco binding
        DNA(2,RPA_Macro_I_Bind_Loc:RPA_Macro_I_Bind_Loc+(n_A-1)) = RPA_A;   %binds RPA-A to lattice
        DNA(2,RPA_Macro_I_Bind_Loc+n_A:RPA_Macro_I_Bind_Loc+(n_RPA-1)) = RPA_D; %binds RPA-D to lattice
        RPA_A_BoundAtSpot(RPA_Macro_I_Bind_Loc) = 1;    %records where RPA-A is currently bound
        RPA_D_BoundAtSpot(RPA_Macro_I_Bind_Loc+n_A) = 1;    %records where RPA-D is currently bound
        LocHist(7,Event) = RPA_Macro_I_Bind_Loc;    %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;0;-1];  %updates free protein counter
    elseif sum(a_Prop(1:8,Event)) >= Randoms(2)*a_0(Event)  %RPA Macro Singly Contiguous Binding
        j(Event) = 8;
        RPA_Macro_SC_Bind_Loc = RPA_SC(randi(numel(RPA_SC)));   %random SC location for RPA Marco binding
        DNA(2,RPA_Macro_SC_Bind_Loc:RPA_Macro_SC_Bind_Loc+(n_A-1)) = RPA_A;   %binds RPA-A to lattice
        DNA(2,RPA_Macro_SC_Bind_Loc+n_A:RPA_Macro_SC_Bind_Loc+(n_RPA-1)) = RPA_D; %binds RPA-D to lattice
        RPA_A_BoundAtSpot(RPA_Macro_SC_Bind_Loc) = 1;    %records where RPA-A is currently bound
        RPA_D_BoundAtSpot(RPA_Macro_SC_Bind_Loc+n_A) = 1;    %records where RPA-D is currently bound
        LocHist(8,Event) = RPA_Macro_SC_Bind_Loc;    %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;0;-1];  %updates free protein counter
    elseif sum(a_Prop(1:9,Event)) >= Randoms(2)*a_0(Event)  %RPA Macro Doubly Contiguous Binding
        j(Event) = 9;
        RPA_Macro_DC_Bind_Loc = RPA_DC(randi(numel(RPA_DC)));   %random DC location for RPA Marco binding
        DNA(2,RPA_Macro_DC_Bind_Loc:RPA_Macro_DC_Bind_Loc+(n_A-1)) = RPA_A;   %binds RPA-A to lattice
        DNA(2,RPA_Macro_DC_Bind_Loc+n_A:RPA_Macro_DC_Bind_Loc+(n_RPA-1)) = RPA_D; %binds RPA-D to lattice
        RPA_A_BoundAtSpot(RPA_Macro_DC_Bind_Loc) = 1;    %records where RPA-A is currently bound
        RPA_D_BoundAtSpot(RPA_Macro_DC_Bind_Loc+n_A) = 1;    %records where RPA-D is currently bound
        LocHist(9,Event) = RPA_Macro_DC_Bind_Loc;    %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;0;-1];  %updates free protein counter
    elseif sum(a_Prop(1:10,Event)) >= Randoms(2)*a_0(Event) %RAD51 Monomer Unbinding
        j(Event) = 10;
        RAD51_Bound_Monomers = find(RAD51_Mon_BoundAtSpot == 1);    %list of all locations where RAD51 monomers are bound
        RAD51_Mon_Unbind_Loc = RAD51_Bound_Monomers(randi(numel(RAD51_Bound_Monomers)));    %random RAD51 monomer to unbind
        DNA(2,RAD51_Mon_Unbind_Loc:RAD51_Mon_Unbind_Loc+(n_RAD51-1)) = 0;   %unbinds RAD51 monomer
        RAD51_Mon_BoundAtSpot(RAD51_Mon_Unbind_Loc) = 0;    %records that protein is no longer currently bound here
        LocHist(10,Event) = RAD51_Mon_Unbind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[1;0;0];    %updates free protein counter
    elseif sum(a_Prop(1:11,Event)) >= Randoms(2)*a_0(Event) %RAD51 Dimer Unbinding
        j(Event) = 11;
        RAD51_Bound_Dimers = Left_RAD51_Dimer_Filament; %list of all locations where a dimer is bound
        RAD51_Dim_Unbind_Loc = RAD51_Bound_Dimers(randi(numel(RAD51_Bound_Dimers)));    %random RAD51 dimer to unbind
        DNA(2,RAD51_Dim_Unbind_Loc:RAD51_Dim_Unbind_Loc+(2*n_RAD51-1)) = 0; %unbinds RAD51 dimer
        RAD51_Mon_BoundAtSpot([RAD51_Dim_Unbind_Loc,RAD51_Dim_Unbind_Loc+n_RAD51]) = 0; %proteins are no longer currently bound here
        LocHist(11,Event) = RAD51_Dim_Unbind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;1;0];    %updates free protein counter
    elseif sum(a_Prop(1:12,Event)) >= Randoms(2)*a_0(Event) %RPA-A Microscopic Binding
        j(Event) = 12;
        RPA_A_Micro_Bind_Loc = Free_for_RPA_A(randi(numel(Free_for_RPA_A)));    %random location for RPA_A to rebind
        DNA(2,RPA_A_Micro_Bind_Loc:RPA_A_Micro_Bind_Loc+(n_A-1)) = RPA_A;   %binds RPA-A back to lattice
        DNA(1,RPA_A_Micro_Bind_Loc:RPA_A_Micro_Bind_Loc+(n_A-1)) = 0;   %hinges RPA-A closed
        RPA_A_BoundAtSpot(RPA_A_Micro_Bind_Loc) = 1;    %shows RPA-A is now currently bound here
        RPA_A_HingedOpen(RPA_A_Micro_Bind_Loc) = 0; %RPA-A is no longer currently hinged open
        LocHist(12,Event) = RPA_A_Micro_Bind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event);  %free protein counter remains the same
    elseif sum(a_Prop(1:13,Event)) >= Randoms(2)*a_0(Event) %RPA-D Microscopic Binding
        j(Event) = 13;
        RPA_D_Micro_Bind_Loc = Free_for_RPA_D(randi(numel(Free_for_RPA_D)));    %random location for RPA-D to rebind
        DNA(2,RPA_D_Micro_Bind_Loc:RPA_D_Micro_Bind_Loc+(n_D-1)) = RPA_D;   %binds RPA-D back to lattice
        DNA(1,RPA_D_Micro_Bind_Loc:RPA_D_Micro_Bind_Loc+(n_D-1)) = 0;   %hinges RPA-D closed
        RPA_D_BoundAtSpot(RPA_D_Micro_Bind_Loc) = 1;    %shows RPA-D is now currently bound here
        RPA_D_HingedOpen(RPA_D_Micro_Bind_Loc) = 0; %RPA-D is no longer currently hinged open
        LocHist(13,Event) = RPA_D_Micro_Bind_Loc;   %records where each event occurs
        Free_Proteins(:,Event+1) = Free_Proteins(:,Event);  %free protein counter remains the same
    elseif sum(a_Prop(1:14,Event)) >= Randoms(2)*a_0(Event) %RPA-D Microscopic Unbinding
        j(Event) = 14;
        RPA_A_Bound_Proteins = find(RPA_A_BoundAtSpot == 1);    %list of all locations where RPA-A is bound
        RPA_A_Micro_Unbind_Loc = RPA_A_Bound_Proteins(randi(numel(RPA_A_Bound_Proteins)));  %random RPA-A to micro dissociate
        DNA(2,RPA_A_Micro_Unbind_Loc:RPA_A_Micro_Unbind_Loc+(n_A-1)) = 0;   %unbinds RPA-A
        DNA(1,RPA_A_Micro_Unbind_Loc:RPA_A_Micro_Unbind_Loc+(n_A-1)) = 1;   %hinges RPA-A open
        RPA_A_BoundAtSpot(RPA_A_Micro_Unbind_Loc) = 0;  %RPA-A is no longer bound here currently
        LocHist(14,Event) = RPA_A_Micro_Unbind_Loc; %records where each event occurs
        if DNA(1,RPA_A_Micro_Unbind_Loc+n_A) == RPA_D    %if both RPA-A and RPA-D are hinged open here...
            DNA(1,RPA_A_Micro_Unbind_Loc:RPA_A_Micro_Unbind_Loc+(n_RPA-1)) = 0;    %RPA protein is now free
            Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;0;1];  %updates free protein counter
        else
            Free_Proteins(:,Event+1) = Free_Proteins(:,Event);  %free protein counter remains the same
        end
    elseif sum(a_Prop(1:15,Event)) >= Randoms(2)*a_0(Event) %RPA-D Microscopic Unbinding
        j(Event) = 15;
        RPA_D_Bound_Proteins = find(RPA_D_BoundAtSpot == 1);    %list of all locations where RPA-D is bound
        RPA_D_Micro_Unbind_Loc = RPA_D_Bound_Proteins(randi(numel(RPA_D_Bound_Proteins)));  %random RPA-D to micro dissociate
        DNA(2,RPA_D_Micro_Unbind_Loc:RPA_D_Micro_Unbind_Loc+(n_D-1)) = 0;   %unbinds RPA-D
        DNA(1,RPA_D_Micro_Unbind_Loc:RPA_D_Micro_Unbind_Loc+(n_D-1)) = 1;   %hinges RPA-D open
        RPA_D_BoundAtSpot(RPA_D_Micro_Unbind_Loc) = 0;  %RPA-D is no longer bound here currently
        LocHist(15,Event) = RPA_D_Micro_Unbind_Loc; %records where each event occurs
        if DNA(1,RPA_D_Micro_Unbind_Loc-1) == RPA_A    %if both RPA-A and RPA-D are hinged open here...
            DNA(1,RPA_D_Micro_Unbind_Loc:RPA_D_Micro_Unbind_Loc+(n_D-1)) = 0;    %RPA protein is now free
            DNA(1,RPA_D_Micro_Unbind_Loc-n_A:RPA_D_Micro_Unbind_Loc-1) = 0;
            Free_Proteins(:,Event+1) = Free_Proteins(:,Event)+[0;0;1];  %updates free protein counter
        else
            Free_Proteins(:,Event+1) = Free_Proteins(:,Event);  %free protein counter remains the same
        end
    end
    
    t(Event+1) = t(Event)+dt(Event);    %advances time
    FracCover_RAD51(Event+1) = numel(find(DNA(2,:) == RAD51))/N;    %RAD51 Saturation
    FracCover_RPA_A(Event+1) = numel(find(DNA(2,:) == RPA_A))/N;    %RPA-A Saturation
    FracCover_RPA_D(Event+1) = numel(find(DNA(2,:) == RPA_D))/N;    %RPA-D Saturation
    FracCover_RPA(Event+1) = FracCover_RPA_A(Event+1)+FracCover_RPA_D(Event+1); %RPA Saturation
    FracCover_Total(Event+1) = FracCover_RPA(Event+1)+FracCover_RAD51(Event+1); %Total Saturation of ssDNA
end

Max_RAD51_Sat = (Free_Proteins(1,1)*n_RAD51)/N; %maximum saturation for RAD51
Max_RPA_A_Sat = (Free_Proteins(3,1)*n_A)/N; %maximum saturation for RPA-A
Max_RPA_D_Sat = (Free_Proteins(3,1)*n_D)/N; %maximum saturation for RPA-D
Max_RPA_Sat = (Free_Proteins(3,1)*n_RPA)/N; %maximum saturation for RPA
Max_Sat = ((Free_Proteins(1,1)*n_RAD51)+(Free_Proteins(3,1)*n_RPA))/N;   %maximum total saturation

figure(1);  %Saturation Plot
P_RAD51 = scatter(t,FracCover_RAD51,1,'red','filled');    hold on;
yline(Max_RAD51_Sat,'--red');
P_RPA_A = scatter(t,FracCover_RPA_A,1,'cyan','filled'); yline(Max_RPA_A_Sat,'--cyan');
P_RPA_D = scatter(t,FracCover_RPA_D,1,'blue','filled'); yline(Max_RPA_D_Sat,'--blue');
P_RPA = scatter(t,FracCover_RPA,1,'magenta','filled');  yline(Max_RPA_Sat,'--magenta');
P_Total = scatter(t,FracCover_Total,1,'k','filled');    yline(Max_Sat,'--k');
xlabel('Time, t'); xlim([0 max(t)]);
ylabel('Saturation'); ylim([0 1]);
title('RAD51/RPA Competition Saturation');
legend([P_RAD51,P_RPA_A,P_RPA_D,P_RPA,P_Total],'RAD51','RPA-A','RPA-D','All RPA','Total','location','southoutside','orientation','horizontal');
box on;

figure(2);    %Free Protein Count Plot
P_FreeRAD51_M = scatter(t,Free_Proteins(1,:),1,'r','filled');   hold on;
P_FreeRAD51_D = scatter(t,Free_Proteins(2,:),1,'b','filled');
P_Free_RPA = scatter(t,Free_Proteins(3,:),1,'g','filled');
xlabel('Time, t'); xlim([0 max(t)]);
ylabel('Population');   ylim([0 sum(Free_Proteins(1))]);
legend([P_FreeRAD51_M,P_FreeRAD51_D,P_Free_RPA],'RAD51 Mon.','RAD51 Dim.','RPA');
box on;