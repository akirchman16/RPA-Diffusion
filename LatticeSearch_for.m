function [Counts,Locs] = LatticeSearch_for(DNA,n_A,n_D,n_RAD51)
    RPA_A = 1;   RPA_D = 3;  RAD51 = 51;
    ProteinDiffs = [RPA_D-RPA_A,RAD51-RPA_D,RAD51-RPA_A];

    N = length(DNA);
    n_RPA = n_A+n_D;
    
    %Find the edges of each gap along the lattice
    GapLeft = find(diff([1 DNA(2,:) 1]) < 0 & ~ismember(abs(diff([1 DNA(2,:) 1])),ProteinDiffs));   %left edge of all gaps
    GapRight = find(diff([1 DNA(2,:) 1]) > 0 & ~ismember(abs(diff([1 DNA(2,:) 1])),ProteinDiffs))-1; %right edge of all gaps
    GapSizes = GapRight-GapLeft+1;  %length of each gap
    
    %Search for RPA binding sites
    RPA_GapEdges = [GapLeft(GapSizes >= n_RPA) GapRight(GapSizes >= n_RPA)];
    RPA_I_Locs = [];    RPA_SC_Locs = [];   RPA_DC_Locs = [];
    RPA_AllLocs = [];
    for i = 1:numel(RPA_GapEdges(1,:))
        RPA_AllLocs = sort([RPA_AllLocs,(RPA_GapEdges(i,:):1:RPA_GapEdges(i,2)-(n_RPA-1))]);
    end
    for j = RPA_AllLocs
        if j == 1
            if DNA(2,1+n_RPA) == RPA_A | DNA(2,1+n_RPA) == RPA_D
                %Singly Contiguous
                RPA_SC_Locs = [RPA_SC_Locs,j];
            else
                %Isolated
                RPA_I_Locs = [RPA_I_Locs,j];
            end
        elseif j == N-(n_RPA-1)
            if DNA(2,j-1) == RPA_A | DNA(2,j-1) == RPA_D
                %Singly Contiguous
                RPA_SC_Locs = [RPA_SC_Locs,j];
            else
                %Isolated
                RPA_I_Locs = [RPA_I_Locs,j];
            end
        else
            if DNA(2,j-1) == RPA_A | DNA(2,j-1) == RPA_D
                %Doubly Contiguous
                RPA_DC_Locs = [RPA_DC_Locs,j];
            elseif xor((DNA(2,j-1) == RPA_A | DNA(2,j-1) == RPA_D),(DNA(2,j+n_RPA) == RPA_A | DNA(2,j+n_RPA) == RPA_D))
                %Singly Contiguous
                RPA_SC_Locs = [RPA_SC_Locs,j];
            else
                %Isolated
                RPA_I_Locs = [RPA_I_Locs,j];
            end
        end
    end
    X_RPA_I = numel(RPA_I_SC); X_RPA_SC = numel(RPA_SC_Locs); X_RPA_DC = numel(RPA_DC_Locs);
    
    %Search for RAD51 Monomer binding sites
    RAD51_Mon_GapEdges = [GapLeft(GapSizes >= n_RAD51) GapRight(GapSizes >= n_RAD51)];
    RAD51_Mon_GapSizes = GapSizes(GapSizes >= n_RAD51);
    
    %Search for RAD51 Dimer binding sites
    RAD51_Dim_GapEdges = [GapLeft(GapSizes >= 2*n_RAD51) GapRight(GapSizes >= 2*n_RAD51)];
    RAD51_Dim_GapSizes = GapSizes(GapSizes >= 2*n_RAD51);
    
end