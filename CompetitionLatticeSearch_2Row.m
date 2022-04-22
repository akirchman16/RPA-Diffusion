% This function will take the current state of the ssDNA lattice and
% calculate the population counts of each type of protein and binding
% location. These will include: bound RAD51 monomer, bound RAD51 dimer, 
% bound RPA-A, bound RPA-D, hinged open RPA-A, hinged open RPA-D, open
% RAD51 monomer location, open RAD51 dimer location, open RPA location, and
% free proteins of each type (all of these can include cooperative binding
% sites as necessary)

function [Counts,Locations] = CompetitionLatticeSearch_2Row(DNA,N,n_RAD51,n_A,n_D)
    n_RPA = n_A+n_D;
% First to find the gaps along the entire lattice and only look at the ones
% large enough to hold any of the proteins. This will give the counts (and
% locations) of places on the lattice where proteins can bind to.
    ProteinDiffs = [2,48,50];    %values of all possible values of differences between protein values
                                 %prevents them from being counted in gaps
                                 
    %Find all the gaps in the bound section of the lattice
    GapLeft = find(diff([1 DNA(2,:) 1]) < 0 & ~ismember(abs(diff([1 DNA(2,:) 1])),ProteinDiffs));   %left edge of all gaps
    GapRight = find(diff([1 DNA(2,:) 1]) > 0 & ~ismember(abs(diff([1 DNA(2,:) 1])),ProteinDiffs))-1; %right edge of all gaps
    
    GapSizes = GapRight-GapLeft+1;  %length of each gap
    
    %Left edge of gaps that are large enough to hold each type of protein
    RPA_GapEdges = [GapLeft(GapSizes >= n_RPA);GapRight(GapSizes >= n_RPA)];
    RAD51_Mon_GapEdges = [GapLeft(GapSizes >= n_RAD51);GapRight(GapSizes >= n_RAD51)];
    RAD51_Dim_GapEdges = [GapLeft(GapSizes >= 2*n_RAD51);GapRight(GapSizes >= 2*n_RAD51)];
    %List every possible binding location for each protein. These will be
    %removed from this array as we record which types they are.
    RPA_AllLocs = [];   RAD51_Mon_AllLocs = []; RAD51_Dim_AllLocs = [];
    for i = 1:numel(RPA_GapEdges(1,:))
        RPA_AllLocs = [RPA_AllLocs,RPA_GapEdges(1,i):1:RPA_GapEdges(2,i)-(n_RPA-1)];
    end
    for j = 1:numel(RAD51_Mon_GapEdges(1,:))
        RAD51_Mon_AllLocs = [RAD51_Mon_AllLocs,RAD51_Mon_GapEdges(1,j):1:RAD51_Mon_GapEdges(2,j)-(n_RAD51-1)];
    end
    for k = 1:numel(RAD51_Dim_GapEdges(2,:))
        RAD51_Dim_AllLocs = [RAD51_Dim_AllLocs,RAD51_Dim_GapEdges(1,k):1:RAD51_Dim_GapEdges(2,k)-(2*n_RAD51-1)];
    end
    
    %Find gaps which hold doubly contiguous binding locations - This is the
    %easiest type of binding location to find (edges aren't allowed)
    RPA_DC_Locs = RPA_GapEdges(RPA_GapEdges(2,:)-RPA_GapEdges(1,:)+1 == n_RPA);
    RAD51_Mon_DC_Locs = RAD51_Mon_GapEdges(RAD51_Mon_GapEdges(2,:)-RAD51_Mon_GapEdges(1,:)+1 == n_RAD51);
    RAD51_Dim_DC_Locs = RAD51_Dim_GapEdges(RAD51_Dim_GapEdges(2,:)-RAD51_Dim_GapEdges(1,:)+1 == 2*n_RAD51);
    %Each of these DC locations needs to be bounded by a similar type of
    %protein (ex: RPA_DC needs to be bounded by RPA on left and right)
    RPA_DC_Locs = RPA_DC_Locs((DNA(2,RPA_DC_Locs-1) == 1 | DNA(2,RPA_DC_Locs-1) == 3) & (DNA(2,RPA_DC_Locs+n_RPA) == 1 | DNA(2,RPA_DC_Locs+n_RPA) == 3));
    RAD51_Mon_DC_Locs = RAD51_Mon_DC_Locs(DNA(2,RAD51_Mon_DC_Locs-1) == 51 & DNA(2,RAD51_Mon_DC_Locs+n_RAD51) == 51);
    RAD51_Dim_DC_Locs = RAD51_Dim_DC_Locs(DNA(2,RAD51_Dim_DC_Locs-1) == 51 & DNA(2,RAD51_Dim_DC_Locs+(2*n_RAD51)) == 51);
    
    RPA_GapEdges(ismember(RPA_DC_Locs,RPA_GapEdges(1,:))) = []; %Clear DC locations from GapEdges arrays
    RAD51_Mon_GapEdges(ismember(RAD51_Mon_DC_Locs,RAD51_Mon_GapEdges(1,:))) = [];
    RAD51_Dim_GapEdges(ismember(RAD51_Dim_DC_Locs,RAD51_Dim_GapEdges(1,:))) = [];
    RPA_AllLocs(ismember(RPA_DC_Locs,RPA_AllLocs)) = [];    %Clears DC locations from AllLocs arrays
    RAD51_Mon_AllLocs(ismember(RAD51_Mon_DC_Locs,RAD51_Mon_AllLocs)) = [];
    RAD51_Dim_AllLocs(ismember(RAD51_Dim_DC_Locs,RAD51_Dim_AllLocs)) = [];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Counts = [numel(RPA_DC_Locs);numel(RAD51_Mon_DC_Locs);numel(RAD51_Dim_DC_Locs)];
    Locations = [RPA_DC_Locs;RAD51_Mon_DC_Locs;RAD51_Dim_DC_Locs];
end