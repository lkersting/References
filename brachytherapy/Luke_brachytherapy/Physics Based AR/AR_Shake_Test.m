% Luke Kersting
% Physics Based Atractiion-Repulsiion Model Test
% 

close all
%clear all

%% LOAD DOSE KERNEL 

% load DoseMatrix
% Seed position in dose rate matrix (a,b,c)
a= 70; b= 70; c= 15;

%% SET UP GRID 
xpts = 60; ypts = 55; zpts = 11;
dx= 0.1; dy= 0.1; dz= 0.5;  % Distance between grid points [cm]

%% PLACE DOSE SINKS 
Dsink=zeros(xpts,ypts,zpts);    %sink dose matrix
% Number of Sinks
Nsinks = 5;

% Sink Positions
sink = zeros(Nsinks,3);
sink(1,:) = [30,20,7];
sink(2,:) = [30,30,5];
sink(3,:) = [20,20,6];
sink(4,:) = [30,30,6];
sink(5,:) = [40,20,6];

for n= 1:Nsinks
    % Shift dose matrix to center on seed
    delta_i= a-sink(n,1);   % index shift in x   
    delta_j= b-sink(n,2);   % index shift in y
    delta_k= c-sink(n,3);   % index shift in z
    % Dose matrix centered on seed location
    centered= circshift(Dose, [-delta_i,-delta_j,-delta_k]);
    % Dose [Gy]
    Dsink= Dsink+centered(1:xpts,1:ypts,1:zpts);
end

%% PLACE SEEDS
Nseed = 5;  % Number of Seeds

% Seed Positions
seed = zeros(Nseed,3);
seed(1,:) = [30,30,7];
seed(2,:) = [30,25,6];
seed(3,:) = [25,25,7];
seed(4,:) = [30,25,6];
seed(5,:) = [35,25,5];


%% ITERATE THROUGH OPTIMIZATION SEQUENCE
for m = 1:30
    %% CALCULATE DOSE
    Dseed=zeros(xpts,ypts,zpts);    %seed dose matrix
    for n= 1:Nseed
        % Shift dose matrix to center on seed
        delta_i= a-seed(n,1);   % index shift in x   
        delta_j= b-seed(n,2);   % index shift in y
        delta_k= c-seed(n,3);   % index shift in z
        % Dose matrix centered on seed location
        centered= circshift(Dose, [-delta_i,-delta_j,-delta_k]);
        % Dose [Gy]
        Dseed= Dseed+centered(1:xpts,1:ypts,1:zpts);
    end
 

    %% CALCULATE AR CHARGE (Q) 
    Q = Dsink - Dseed;
    if Q == 0
        break
    end

    %% CALCULATE ATTRACTION-REPULSION POTENTIAL (V) and FIELD (E)
    V= zeros(Nseed,3);          % Potential Vector Matrix
    E= zeros(Nseed,3);          % Field Vector Matrix 
    for k= 1:zpts     % z voxels
    for j= 1:ypts     % y voxels
    for i= 1:xpts 	  % x voxels

       % Compute distance between voxel(i,j,k) and the seed
       r= [(seed(:,1)-i)*dx, (seed(:,2)-j)*dy, (seed(:,3)-k)*dz];
       % 1/r between seed and voxel(i,j,k)
       U= diag(1./sqrt(sum(r.^2,2))).^2*r;
       % 1/r^2 between seed and voxel(i,j,k)
       U2= diag(1./sqrt(sum(r.^2,2))).^3*r;

       % Voxel located at seed's position
       ind= r==0;
       U(ind)= 0;
       U2(ind)= 0;

       % AR Potential at seed
       V= V + Q(i,j,k).*U;
       % AR Field at seed
       E= E + Q(i,j,k).*U2;
    end
    end
    end
    %ind = abs(V) <= .1;
    %V(ind) = 0;
    %E(ind) = 0;

    %% SELECT SEED TO MOVE
    % select index of seed with largest V magnitude
    maxV =  max(sqrt(sum(V.^2,2))) == sqrt(sum(V.^2,2));
    % select index of seed with largest E magnitude
    maxE =  max(sqrt(sum(E.^2,2))) == sqrt(sum(E.^2,2));
    %maxE = [1,1,0]
    if sum(maxE) > 1
        same= find(maxE==1);
        maxE(same(1))=0;
    end
    
    %% SELECT DIRECTION TO MOVE SEED
    dirV = V(maxV,:)/norm(V(maxV,:));
    dirV2 = V(maxV,:)/max(abs(V(maxV,:)));
    dirE = E(maxE,:)/norm(E(maxE,:));
    dirE2 = E(maxE,:)/max(abs(E(maxE,:)));

    %% MOVE SEED WITHIN TEMPLATE
    g = 1;                          % Gain (integer value only)
    dr = -g*round(dirE);            % Change in seed temlplate position
    dr(1:2) = dr(1:2)*1;            % Distance between template location is 5 voxels
    seed(maxE,:) = seed(maxE,:) + dr% New seed poisition
end




