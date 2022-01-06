Graph_DNA(Event+1,:) = DNA(1,:); %HAS TO BE IN FOR LOOP FOR EACH ITERATION

Graph_DNA(Graph_DNA == 0) = 1;      %Empty locations on DNA lattice (white)
Graph_DNA(Graph_DNA == RPA_A) = 2;  %Locations where RPA-A is bound (cyan)
Graph_DNA(Graph_DNA == RPA_D) = 3;  %Locations where RPA-D is bound (blue)
Graph_DNA(Graph_DNA == RAD51) = 4;  %Locations where RAD51 is bound (red)

[X_DNA,Y_Time] = meshgrid(1:N,t);
CustMap = [1, 1, 1; 0, 1, 1; 0, 0, 1; 1, 0, 0];  %custom color range for corresponding proteins (colors labeled above)
colormap(figure(3),CustMap);

figure(3);  %shows proteins moving over time
surf(X_DNA,Y_Time,Graph_DNA,'EdgeColor','none');
set(gca,'Ydir','reverse');  %reverses time axis so beginning is top of figure
view(2);
hold on;
xlabel('ssDNA Location');
xlim([1 N]);
ylabel('Time, t (Inverse)');
ylim([0 max(t)]);
box on;
title('Protein Diffusion');