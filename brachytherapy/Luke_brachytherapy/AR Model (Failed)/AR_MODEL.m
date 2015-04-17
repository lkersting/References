close all
clear all
rng(511)


%%   LOAD DATA   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
fid= fopen('organ_mesh-SUA');    %Opens Sua's organ data file
dimensions= textscan(fid, '%s %n',6, 'delimiter', ':');    %read dimensions in header
%---Number of mesh points [x y z]
xpts= dimensions{2}(1);  %Number of x mesh points
ypts= dimensions{2}(3);  %Number of y mesh points
zpts= dimensions{2}(5);  %Number of z mesh points
Npts= xpts*ypts*zpts;    %Total number of mesh points
%---Distance between mesh points
dx= dimensions{2}(2);    %distance between x mesh points
dy= dimensions{2}(4);    %distance between y mesh points
dz= dimensions{2}(6);    %distance between z mesh points

%---Target Data
Ntarget= textscan(fid, '%*s %n',1);   %mesh points that make up target
target= zeros(xpts,ypts,zpts);        %[0 1] matrix of target organ
%target_c= cell(zpts,1);
%---Read target data from organ_mesh-SUA file
for i= 1:zpts
   data= textscan(fid,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   target(:,:,i)= permute(data{1},[2 1 3]);
   %target_c{i}= data{1};
end

%---Urethra Data
Nurethra= textscan(fid, '%*s %n',1);
urethra= zeros(xpts,ypts,zpts);
%urethra_c= cell(zpts,1);
%---Read urethra data from organ_mesh-SUA file
for i= 1:zpts
   data= textscan(fid,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   urethra(:,:,i)= permute(data{1},[2 1 3]);
   %urethra_c{i}= data{1};
end

%---Rectum Data
Nrectum= textscan(fid, '%*s %n',1);
rectum= zeros(xpts,ypts,zpts);
%rectum_c= cell(zpts,1);
%---Read rectum data from organ_mesh-SUA file
for i= 1:zpts
   data= textscan(fid,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   rectum(:,:,i)= permute(data{1},[2 1 3]);
%   rectum_c{i}= data{1};
end

%---Normal Data
Nnormal= textscan(fid, '%*s %n',1);
normal= zeros(xpts,ypts,zpts);
%normal_c= cell(zpts,1);

for i= 1:zpts
   data= textscan(fid,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   normal(:,:,i)= permute(data{1},[2 1 3]);
%   normal_c{i}= data{1};
end

%---Margin Data
Nmargin= textscan(fid, '%*s %n',1);
margin= zeros(xpts,ypts,zpts);
%margin_c= cell(zpts,1);
%---Read margin data from organ_mesh-SUA file
for i= 1:zpts
   data= textscan(fid,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   margin(:,:,i)= permute(data{1},[2 1 3]);
%   margin_c{i}= data{1};
end

%pretty picture
% % target_con= target(:,:,1);
% % for i= 2:zpts
% %     target_con= target_con+target(:,:,i);
% % end
% % 
% % surf(target_con)

fid2= fopen('needle_template-SUA');    %Opens Sua's needle template file
dim= textscan(fid2, '%s %n',6, 'delimiter', ':');    %read dimensions in header

%---Template Data
template= zeros(xpts,ypts,zpts);        %[0 1] matrix of possible seed locations
%template_c= cell(zpts,1);

%---Read template data from organ_mesh-SUA file
for i= 1:zpts
   data= textscan(fid2,repmat('%n',1,xpts),ypts,'CollectOutput',1);
   template(:,:,i)= permute(data{1},[2 1 3]);
   %template_c{i}= data{1};
end
%---Possible x and y compnents of needle locations
Ntemplate= sum(sum(template(:,:,5)));
[u,v]= find(template(:,:,5));


%%   PLACE SEEDS   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Random needle placement
Nneedles= 13;    %number of initial needles
rnd1= round(rand(Nneedles,1)*Ntemplate);
needles= [u(rnd1) v(rnd1)];
newneedles= zeros(Nneedles,2);

%---Find the number of initial seeds
Nseeds= 0;       
for i= 1:Nneedles
    [o]= find(template(needles(i,1),needles(i,2),:));
    Nseeds= Nseeds+length(o);
end

%---Initial seed placement
%NOTE: there are mutiple seeds for each x and y poisiton of a needle
seeds= ones(Nseeds,3);
newseeds= ones(Nseeds,3);
num_ndl=zeros(Nneedles,1);
mm= 0;
for i= 1:Nneedles
    [o]= find(template(needles(i,1),needles(i,2),:));    %z position along needle i
    mn= mm+length(o);                            
    seeds(mm+1:mn,1)= seeds(mm+1:mn,1)*needles(i,1);     %x position
    seeds(mm+1:mn,2)= seeds(mm+1:mn,2)*needles(i,2);     %y position
    seeds(mm+1:mn,3)= seeds(mm+1:mn,3).*o;               %z poistion
    num_ndl(i)=length(o);                                %number of seeds in the needle
    mm= mn;
end

%%   OPTIMIZE SEEDS   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=zeros(xpts,ypts,zpts);    %dose-rate matrix


%---Target Objective Function, f_T
f_T= zeros(xpts,ypts,zpts); %Target objective function
L_T= 180;            % Left dose boundary
Gra_L_T= -40;        % Slope to left boundary
R_T= 255;            % Right dose boundary
Gra_R_T= -8;         % Slope to right boundary

%--- Urethra Objective Function, f_U
f_U= zeros(xpts,ypts,zpts); %Urethra objective function
L_U= 95;            % Left dose boundary
Gra_L_U= -40;        % Slope to left boundary
R_U= 100;            % Right dose boundary
Gra_R_U= -40;         % Slope to right boundary

%--- Rectum Objective Function, f_R
f_R= zeros(xpts,ypts,zpts); %Rectum objective function
R_R= 80;            % Right dose boundary
Gra_R_R= -40;        % Slope to right boundary

%--- Normal Tissue Objective Function, f_N
f_N= zeros(xpts,ypts,zpts);
R_N= 145;            % Right dose boundary
Gra_R_N= -40;        % Slope to right boundary

%Calculate dose

%---Use SeedDose.m to calculate dose by how many voxels away 
%   from the seed in the x y z direction 
load DoseMatrix
% Seed position in dose rate matrix (a,b,c)
a= 70; b= 70; c= 15;

optimized=1;
iter=0;
iter_end=10;
while optimized
for s= 1:Nseeds
    % Shift dose matrix to center on seed
    delta_i= a-seeds(1,1);   % index shift in x   
    delta_j= b-seeds(1,2);   % index shift in y
    delta_k= c-seeds(1,3);   % index shift in z
    % Dose matrix centered on seed location
    centered= circshift(Dose, [-delta_i,-delta_j,-delta_k]);
    % Dose [Gy]
    D= D+centered(1:xpts,1:ypts,1:zpts);
end

% f_T for underdose
TT= (D < L_T) & (target > 0.5);
f_T(TT)= Gra_L_T*(D(TT)-L_T);

% f_T for overdose
TT= (D > R_T) & (target > 0.5);
f_T(TT)= Gra_R_T*(D(TT)-R_T);


% f_U for underdose
U= (D < L_U) & (urethra > 0.5);
f_U(U)= Gra_L_U*(D(U)-L_U);

% f_U for overdose
U= (D > R_U) & (urethra > 0.5);
f_U(U)= Gra_R_U*(D(U)-R_U);


% f_R for overdose
R= (D > R_R) & (rectum > 0.5);
f_R(R)= Gra_R_R*(D(R)-R_R);


% f_N for overdose
N= (D > R_N) & (normal > 0.5);
f_N(N)= Gra_R_N*(D(N)-R_N);

%--- Total Objective Function, f
f= f_T + f_U + f_R + f_N; 

if f==0
    optimized=0;
    disp('all perameters met')
end

%--- Displacement Vector, T
T= zeros(Nseeds,3);         %Displacement Vector
%dx= 0.1; dy= 0.1; dz= 0.5;  % Distance between grid points [cm]

for k= 1:zpts     % z voxels
for j= 1:ypts     % y voxels
for i= 1:xpts 	  % x voxels

   % Compute distance (voxels) between voxel(i,j,k) and the seeds
   r= [(seeds(:,1)-i), (seeds(:,2)-j), (seeds(:,3)-k)];
   % Voxel located at seeds' position
   ind=r==0;
   if ~isempty(ind)
       r(ind)=0.1;
   end
   % Unit vectors between seeds and voxel(i,j,k)
   U= r/norm(r);
   T= T + f(i,j,k).*U;
end
end
end

%--- Move needles and seeds
mm= 0;
k= 1E-6;     % Gain factor

for i= 1:Nneedles
    mn= mm+num_ndl(i);
%     %--- Move needles in x direction
%     dx= k*sum(T(mm+1:mn,1))/num_ndl(i);
%     if dx > 5/2     % Half the distance between possible needle i locations 
%         newseeds(mm+1:mn,1,1)= seeds(mm+1:mn,1,1)+round(dx);
%         needles(i,1)= needles(i,1)+round(dx);
%     elseif dx < -5/2
%         %seeds(mm+1:mn,1,1)= seeds(mm+1:mn,1,1)-5;
%         %needles(i,1)= needles(i,1)-5;
%     end
%     %--- Move needles in y direction
%     dy= k*sum(T(mm+1:mn,2))/num_ndl(i);
%     if dy > 5/2     % Half the distance between possible needle j locations
%         newseeds(mm+1:mn,2,1)= seeds(mm+1:mn,2,1)+round(dy);
%         needles(i,2)= needles(i,2)+round(dy);
%     elseif dy < -5/2
%         seeds(mm+1:mn,2,1)= seeds(mm+1:mn,2,1)-5;
%         needles(i,2)= needles(i,2)-5;
%     end
%     %--- Move seeds to new z locations
%     [o]= find(template(needles(i,1),needles(i,2),:));
%     for ii= mm+1:mn
%         dz= k*T(ii,3);
%         if dz > 1/2
%             seeds(ii,3)= seeds(ii,3)+round(dz);
%         elseif dz < -1/2
%             seeds(ii,3)= seeds(ii,3)-1;            
%         end
%     end  
%--- Move needles in x direction
    dx(i,1)= k*sum(T(mm+1:mn,1))/num_ndl(i); 
    newseeds(mm+1:mn,1,1)= seeds(mm+1:mn,1,1)+round(dx(i,1));
    newneedles(i,1)= needles(i,1)+round(dx(i,1));
%--- Move needles in y direction    
    dy(i,1)= k*sum(T(mm+1:mn,2))/num_ndl(i,1);
    newseeds(mm+1:mn,2,1)= seeds(mm+1:mn,2,1)+round(dy(i,1));
    newneedles(i,2)= needles(i,2)+round(dy(i,1));
    
    mm= mn;
end     

%--- Move needles in z direction
dz= k*T(:,3);
newseeds(:,3)= seeds(:,3)+round(dz);

%--- If seeds don't move inbetween iterations then end loop
if newseeds==seeds
    optimized= 0;
    disp('seeds dont move')
end

% Make sure needle positions don't overlap
%n_overlap= zeros(Nneedles,1);
%newneedles(2,:)=newneedles(7,:);
%newneedles(3,:)=newneedles(2,:);
for i= 3
    ind= find((newneedles(i,1)==newneedles(:,1)) & (newneedles(i,2)==newneedles(:,2)))
    n_overlap(1:length(ind),1)= ind
    same=n_overlap(:)==i;
    n_overlap(same)=0;
    if any(n_overlap)
        n_overlap(same)=i;
        XX= abs(abs(newneedles(n_overlap,1))-abs(needles(n_overlap,1)+dx(n_overlap)));             
        YY= abs(abs(newneedles(n_overlap,2))-abs(needles(n_overlap,2)+dy(n_overlap)));  
        RR= XX.^2 + YY.^2;
        n_overlap(RR==min(RR))= 0;
        XX(RR==min(RR))= [];     YY(RR==min(RR))= [];
        while ~isempty(n_overlap)
        if min(XX) >= min(YY)
            ind2= n_overlap(XX==min(XX))
            ind= (newneedles(ind2,1)-1==newneedles(:,1)) & (newneedles(ind2,2)==newneedles(:,2));
            if isempty(ind)
                newneedles(ind2,1)= newneedles(ind2,1)-1;
            else
                ind= (newneedles(ind2,1)+1==newneedles(:,1)) & (newneedles(ind2,2)==newneedles(:,2));
                if isempty(ind)
                    newneedles(ind2,1)= newneedles(ind2,1)+1;
                else 
                    newneedles(ind2,1)= newneedles(ind2,1)-2;
                end
            end
            n_overlap(XX==min(XX))=[];
            XX(XX==min(XX))=[]; YY(XX==min(XX))=[];
            
        else
            ind2= n_overlap(find(YY==min(YY),1));
            open= 1;
            while open > 0.5
                ind3= (newneedles(ind2,1)==newneedles(:,1)) & (newneedles(ind2,2)-1==newneedles(:,2));
                if isempty(ind3)
                    newneedles(ind2,2)= newneedles(ind2,2)-open;
                    open= 0;
                else
                    ind4= (newneedles(ind2,1)==newneedles(:,1)) & (newneedles(ind2,2)+1==newneedles(:,2));
                    if isempty(ind4)
                        newneedles(ind2,2)= newneedles(ind2,2)+open;
                        open=0;
                    elseif isempty(n_overlap)
                        open=0;
                    else
                        open=open+1;
                    end
                end
                n_overlap(YY==min(YY))=[];
                XX(YY==min(YY))=[]; YY(YY==min(YY))=[];
            end
        end
        end
    end
end



%seed=seeds;
% Set seed positions to the new positions
needles=newneedles;
seeds=newseeds;
% Seeds (needle) leaves grid in x direction
outx= (seeds(:,1) > xpts) | (seeds(:,1) < 0);
if ~isempty(outx)
seeds(outx,:)= [];
Nseeds= length(seeds);
end
outx_n= (needles(:,1) > xpts) | (needles(:,1) < 0);
if ~isempty(outx_n)
needles(outx_n,:)= [];
Nneedles= length(needles);
num_ndl= zeros(Nneedles,1);
end

% Seeds (needle) leaves grid in y direction
outy= (seeds(:,2) > ypts) | (seeds(:,2) < 0);
if ~isempty(outy)
seeds(outy,:)= [];
Nseeds= length(seeds);
end
outy_n= (needles(:,2) > ypts) | (needles(:,2) < 0);
if ~isempty(outy_n)
needles(outy_n,:)= [];
Nneedles= length(needles);
num_ndl= zeros(Nneedles,1);
end

% Seeds leaves grid in z direction
outz= (seeds(:,3) > zpts) | (seeds(:,3) < 0);
if ~isempty(outz)
seeds(outz,:)= [];
Nseeds= length(seeds);
end

%lng= (ismember(seeds(:,1),needles(:,1))) & (ismember(seeds(:,1),needles(:,1)));

iter= iter+1;
if iter==1
    optimized= 0;
end
% Display # of iterations
say=['number of iterations = ',num2str(iter)];
disp(say)
end

TDose= 145;             % [Gy] target dose
ind= target > 0.5;
Dtarget= D(ind)/TDose;  % Percent of target dose to the target