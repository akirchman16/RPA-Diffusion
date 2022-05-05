% This function will take the current state of the ssDNA lattice and
% calculate the population counts of each type of protein and binding
% location. These will include: bound RAD51 monomer, bound RAD51 dimer, 
% bound RPA-A, bound RPA-D, hinged open RPA-A, hinged open RPA-D, open
% RAD51 monomer location, open RAD51 dimer location, open RPA location, and
% free proteins of each type (all of these can include cooperative binding
% sites as necessary). This will be done using clustering techniques with
% bwlabel and property functions.

% Columns of 'Count' output and what the count corresponds to:
%     1: RPA Isolated Sites
%     2: RPA Singly Contiguous Sites
%     3: RPA Doubly Contiguous Sites
%     4: RAD51 Monomer Isolated Sites
%     5: RAD51 Monomer Singly Congiguous Sites
%     6: RAD51 Monomer Doubly Contiguous Sites
%     7: RAD51 Dimer Isolated Sites
%     8: RAD51 Dimer Singly Contiguous Sites
%     9: RAD51 Dimer Doubly Contiguous Sites

function [Counts,RPA_I,RPA_SC,RPA_DC,RAD51_Mon_I,RAD51_Mon_SC,RAD51_Mon_DC,RAD51_Dim_I,RAD51_Dim_SC,RAD51_Dim_DC] = LatticeSearch_Cluster(DNA,n_RAD51,n_A,n_D)
    n_RPA = n_A+n_D;
    
%Find all the gaps in the lattice and record their edges and sizes,
%respectively.
    Blank_Cluster = logical(DNA(2,:) == 0);
        %BlankCluster: array showing all the gaps in the DNA lattice
        %GapCount: number of gaps across the whole lattice  
    bb_Empty = regionprops(Blank_Cluster,'BoundingBox');    bb_Empty = cat(1,bb_Empty.BoundingBox);
    GapSize = bb_Empty(:,3);
    GapEdges = [bb_Empty(:,1)+0.5 bb_Empty(:,1)+bb_Empty(:,3)-0.5]; %edges of all gaps [left right]
    
%Find all the gaps which are large enough to fit each type of protein
    % 1-RPA
    RPA_GapEdges = GapEdges(GapSize >= n_RPA,:);
    % 2-RAD51 Monomer
    RAD51_Mon_GapEdges = GapEdges(GapSize >= n_RAD51,:);
    % 3-RAD51 Dimer
    RAD51_Dim_GapEdges = GapEdges(GapSize >= 2*n_RAD51,:);
    %All possible binding locations for each protein. These will be cleared
    %as we find contiguous sites so that the leftovers will be isolated
    %sites.
    RPA_AllLocs = [];   RAD51_Mon_AllLocs = []; RAD51_Dim_AllLocs = [];
    for i = 1:numel(RPA_GapEdges(:,1))
        RPA_AllLocs = [RPA_AllLocs,RPA_GapEdges(i,1):1:(RPA_GapEdges(i,2)-(n_RPA-1))];
    end
    for j = 1:numel(RAD51_Mon_GapEdges(:,1))
        RAD51_Mon_AllLocs = [RAD51_Mon_AllLocs,RAD51_Mon_GapEdges(j,1):1:(RAD51_Mon_GapEdges(j,2)-(n_RAD51-1))];
    end
    for k = 1:numel(RAD51_Dim_GapEdges(:,1))
        RAD51_Dim_AllLocs = [RAD51_Dim_AllLocs,RAD51_Dim_GapEdges(k,1):1:(RAD51_Dim_GapEdges(k,2)-(2*n_RAD51-1))];
    end
    
%Find all clusters of each type of protein and record their left and
%right edges. These values will be helpful when checking for
%cooperative binding sites.
    % 1-RPA
    RPA_Cluster = logical(DNA(2,:) == 1 | DNA(2,:) == 3);
    bb_RPA = regionprops(RPA_Cluster,'BoundingBox'); bb_RPA = cat(1,bb_RPA.BoundingBox);
    if ~isempty(bb_RPA)
        RPA_Edges = [bb_RPA(:,1)+0.5 bb_RPA(:,1)+bb_RPA(:,3)-0.5];
    else
        RPA_Edges = [];
    end
    % 2-RAD51
    RAD51_Cluster = logical(DNA(2,:) == 51);
    bb_RAD51 = regionprops(RAD51_Cluster,'BoundingBox'); bb_RAD51 = cat(1,bb_RAD51.BoundingBox);
    if ~isempty(bb_RPA)
        RAD51_Edges = [bb_RAD51(:,1)+0.5 bb_RAD51(:,1)+bb_RAD51(:,3)-0.5];
    else
        RAD51_Edges = [];
    end
    
% Record Doubly Contiguous Sites!
    %A site is only Doubly Contiguous if the gap is the same size as the
    %corresponding protein and is doubly bodered by that same protein (this
    %is where my previous bug arose in the old system).
    %First off, record the edges of the gaps which are possible places for
    %Doubly Contiguous binding sites.
    RPA_DC_Gaps = RPA_GapEdges((ismember(RPA_GapEdges(:,1)-1,RPA_Edges)) & (ismember(RPA_GapEdges(:,2)+1,RPA_Edges)) & (RPA_GapEdges(:,2)-RPA_GapEdges(:,1)+1 == n_RPA),:);
    RAD51_Mon_DC_Gaps = RAD51_Mon_GapEdges((ismember(RAD51_Mon_GapEdges(:,1)-1,RAD51_Edges)) & (ismember(RAD51_Mon_GapEdges(:,2)+1,RAD51_Edges)) & (RAD51_Mon_GapEdges(:,2)-RAD51_Mon_GapEdges(:,1)+1 == n_RAD51),:);
    RAD51_Dim_DC_Gaps = RAD51_Dim_GapEdges((ismember(RAD51_Dim_GapEdges(:,1)-1,RAD51_Edges)) & (ismember(RAD51_Dim_GapEdges(:,2)+1,RAD51_Edges)) & (RAD51_Dim_GapEdges(:,2)-RAD51_Dim_GapEdges(:,1)+1 == 2*n_RAD51),:);
    %Now record the actual binding sites (left edge of the gap)
    RPA_DC_Sites = RPA_DC_Gaps(:,1);    RAD51_Mon_DC_Sites = RAD51_Mon_DC_Gaps(:,1);    RAD51_Dim_DC_Sites = RAD51_Dim_DC_Gaps(:,1);
    %Clear DC sites from the array of all possible locations
    RPA_AllLocs(ismember(RPA_AllLocs,RPA_DC_Sites)) = [];   RAD51_Mon_AllLocs(ismember(RAD51_Mon_AllLocs,RAD51_Mon_DC_Sites)) = []; RAD51_Dim_AllLocs(ismember(RAD51_Dim_AllLocs,RAD51_Dim_DC_Sites)) = [];
    
% Now look into Singly Contiguous Sites!
    %A site is Singly Contiguous if the gap size is greater than or equal
    %to the corresponding protein and is singly bounded by the same protein
    %(only one edge has the same protein next to it).
    RPA_SC_Gaps = RPA_GapEdges(xor(ismember(RPA_GapEdges(:,1)-1,RPA_Edges),ismember(RPA_GapEdges(:,2)+1,RPA_Edges)) & (RPA_GapEdges(:,2)-RPA_GapEdges(:,1)+1 >= n_RPA),:);
    RAD51_Mon_SC_Gaps = RAD51_Mon_GapEdges(xor(ismember(RAD51_Mon_GapEdges(:,1)-1,RAD51_Edges),ismember(RAD51_Mon_GapEdges(:,2)+1,RAD51_Edges)) & (RAD51_Mon_GapEdges(:,2)-RAD51_Mon_GapEdges(:,1)+1 >= n_RAD51),:);
    RAD51_Dim_SC_Gaps = RAD51_Dim_GapEdges(xor(ismember(RAD51_Dim_GapEdges(:,1)-1,RAD51_Edges),ismember(RAD51_Dim_GapEdges(:,2)+1,RAD51_Edges)) & (RAD51_Dim_GapEdges(:,2)-RAD51_Dim_GapEdges(:,1)+1 >= 2*n_RAD51),:);
    %Now record the actual binding sites
    RPA_SC_Sites = RPA_SC_Gaps(ismember(RPA_SC_Gaps(:,1)-1,RPA_Edges),1);   RPA_SC_Sites = sort([RPA_SC_Sites;RPA_SC_Gaps(ismember(RPA_SC_Gaps(:,2)+1,RPA_Edges),2)-(n_RPA-1)]);
    RAD51_Mon_SC_Sites = RAD51_Mon_SC_Gaps(ismember(RAD51_Mon_SC_Gaps(:,1)-1,RAD51_Edges),1);   RAD51_Mon_SC_Sites = sort([RAD51_Mon_SC_Sites;RAD51_Mon_SC_Gaps(ismember(RAD51_Mon_SC_Gaps(:,2)+1,RAD51_Edges),2)-(n_RAD51-1)]);
    RAD51_Dim_SC_Sites = RAD51_Dim_SC_Gaps(ismember(RAD51_Dim_SC_Gaps(:,1)-1,RAD51_Edges),1);   RAD51_Dim_SC_Sites = sort([RAD51_Dim_SC_Sites;RAD51_Dim_SC_Gaps(ismember(RAD51_Dim_SC_Gaps(:,2)+1,RAD51_Edges),2)-(2*n_RAD51-1)]);
    %Clear SC sites from the array of all possible locations
    RPA_AllLocs(ismember(RPA_AllLocs,RPA_SC_Sites)) = [];   RAD51_Mon_AllLocs(ismember(RAD51_Mon_AllLocs,RAD51_Mon_SC_Sites)) = []; RAD51_Dim_AllLocs(ismember(RAD51_Dim_AllLocs,RAD51_Dim_SC_Sites)) = [];
    
% Finally, we have to look at Isolated Sites!
    %Isolated sites are all other possible binding sites that aren't
    %considered Doubly Contiguous or Singly Contiguous. Therefore, the gaps
    %are all of the leftover gaps that haven't been recorded as singly or
    %doubly contiguous gaps.
    RPA_I_Sites = RPA_AllLocs.';    RAD51_Mon_I_Sites = RAD51_Mon_AllLocs.';    RAD51_Dim_I_Sites = RAD51_Dim_AllLocs.';
    
% RECORD ALL RESULTS AND COUNTS (recorded in ways that are useful in the
% full model)
    RPA_I = RPA_I_Sites.';  RPA_SC = RPA_SC_Sites.';    RPA_DC = RPA_DC_Sites.';
    RAD51_Mon_I = RAD51_Mon_I_Sites.'; RAD51_Mon_SC = RAD51_Mon_SC_Sites.'; RAD51_Mon_DC = RAD51_Mon_DC_Sites.';
    RAD51_Dim_I = RAD51_Dim_I_Sites.'; RAD51_Dim_SC = RAD51_Dim_SC_Sites.'; RAD51_Dim_DC = RAD51_Dim_DC_Sites.';
    Counts = [numel(RPA_I),numel(RPA_SC),numel(RPA_DC),numel(RAD51_Mon_I),numel(RAD51_Mon_SC),numel(RAD51_Mon_DC),numel(RAD51_Dim_I),numel(RAD51_Dim_SC),numel(RAD51_Dim_DC)];
end