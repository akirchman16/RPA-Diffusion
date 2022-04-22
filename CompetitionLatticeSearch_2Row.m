% This function will take the current state of the ssDNA lattice and
% calculate the population counts of each type of protein and binding
% location. These will include: bound RAD51 monomer, bound RAD51 dimer, 
% bound RPA-A, bound RPA-D, hinged open RPA-A, hinged open RPA-D, open
% RAD51 monomer location, open RAD51 dimer location, open RPA location, and
% free proteins of each type (all of these can include cooperative binding
% sites as necessary)

function [RPA_GapEdges,RAD51_Mon_GapEdges,RAD51_Dim_GapEdges] = CompetitionLatticeSearch_2Row(DNA,N,n_RAD51,n_A,n_D)
    n_RPA = n_A+n_D;
% First to find the gaps along the entire lattice and only look at the ones
% large enough to hold any of the proteins
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
    
    %Find gaps which hold doubly contiguous binding locations - This is the
    %easiest type of binding location to find
%     RPA_DC_Locs = RPA_GapEdges(RPA_GapEdges(:,2)-RPA_GapEdges(:,1)+1 == n_RPA);
%     RAD51_Mon_DC_Locs = RAD51_Mon_GapEdges(RAD51_Mon_GapEdges(:,2)-RAD51_Mon_GapEdges(:,1)+1 == n_RAD51);
%     RAD51_Dim_DC_Locs = RAD51_Dim_GapEdges(RAD51_Dim_GapEdges(:,2)-RAD51_Dim_GapEdges(:,1)+1 == 2*n_RAD51);
    
    %Find left edges of gaps which hold doubly contiguous binding locations
    % - This is the easiest type of binding location to find
    RPA_DC_Locs = GapLeft(GapRight(GapSizes >= n_RPA)-GapLeft(GapSizes >= n_RPA)+1 == n_RPA);
    RAD51_Mon_DC_Locs = GapLeft(GapRight(GapSizes >= n_RAD51)-GapLeft(GapSizes >= n_RAD51)+1 == n_RAD51);
    RAD51_Dim_DC_Locs = GapLeft(GapRight(GapSizes >= 2*n_RAD51)-GapLeft(GapSizes >= 2*n_RAD51)+1 == 2*n_RAD51);
end