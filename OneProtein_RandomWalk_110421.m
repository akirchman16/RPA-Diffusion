clearvars;
close all;

% This code consists of a single RPA molecule binding to a short (<200 nt)
% ssDNA lattice and then diffusing left and right along the lattice. The
% diffusion will be modeled as stochastic, selecting a time interval from
% an exponential distribution with mean and std dev equal to a diffusion
% coefficient (k_D). Left vs Right diffusion willbe determined by 
% probability values.

N = 100;    %length of ssDNA lattice
n = 3;  %length of RPA protein

k_D = 0.005;   %diffusion coefficient (used for time selection)
L_Prob = 0.5;   %probability of diffusion to the left
R_Prob = 1-L_Prob;  %probability of diffusion to the right
MaxTime = 100;  %maximum time for simulation

InitialBind = randi([1,N-n+1]); %random location for protein to bind to initially
t(1) = 0;

Location = InitialBind;
Event = 0;  %counts events
while t(end) <= MaxTime
    Event = Event+1;    %increases event counter
    dt(Event) = random('Exponential',k_D);  %selects random time interval
    r = rand;   %random value for Monte Carlo step of direction
    if r <= L_Prob
        Location = [Location,Location(Event)-1];    %diffusion to the left
        if Location(end) < 1
            Location(end) = randi([1,N-n+1]);   %new protein binds at new location if protein diffuses off end
        end
    else
        Location = [Location,Location(Event)+1];    %diffusion to the right
        if Location(end) > N-n+1
            Location(end) = randi([1,N-n+1]);    %new potein binds at new location if old protein diffuses off end of lattice
        end
    end
    t(Event+1) = t(Event)+dt(Event);    %advances timer
end

figure(1);
scatter(t,Location,10,'b','filled');
hold on;
% yline(InitialBind,':k');
xlabel('Time, t (s)');
ylabel('Location');
xlim([0 t(end)]);
ylim([0 N]);
title(['1 Protein Random Walk Location (k_D = ', num2str(k_D), ')']);
box on;