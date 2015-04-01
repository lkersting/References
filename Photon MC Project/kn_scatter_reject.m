% This script replicates the angular distribution expected for photons
% scattering off electrons (Klein-Nishina angular differential scattering
% cross-section) using the rejection technique.

clear; clf; clc;
ntrials = 1000000; counts = zeros(ntrials,5); xc = zeros(200,5); nc = zeros(200,5);

for i = 1:200                         % this defines 200 bin centers in increments of 0.01
    xc(i,:) = -0.995 + 0.01*(i-1);  
end

% initialize alpha; the values of alpha specified below are for
% illustration; they could be changed to anything, but are physically equal
% to the ratio of initial photon energy to electron rest mass energy
alpha = zeros(5,1); alpha(1) = 0; alpha(2) = 0.5 ; alpha(3) = 1;
alpha(4) = 5; alpha(5) = 10;

for j = 1:5
   i_counter = 0;                     % initalize the number of kept values
   for i=1:ntrials
       xi_1 = rand(1);                    
       xi_2 = rand(1);
       mu = -1 + 2*xi_1;              % this converts the random number, which
                                      % lies in the interval 0 < xi < 1, to a range
                                      % appropriate for the cosine of the scattering angle
       trm1 = (1 + mu^2)/2;
       trm2 = 1/(1 + alpha(j)*(1 - mu));
       trm3 = 2*alpha(j)^2*(1 - mu)^2*trm2/trm1;
       ang_diff_scat = trm1*trm2^2*(1 + trm3);
       if xi_2 < ang_diff_scat        % if this criterion is satisfied, increment the
           i_counter = i_counter + 1; % number of kept values and count mu
           counts(i_counter,j) = mu;
       end
   end
   [nc(:,j),xout(:,j)] = hist(counts(1:i_counter,j),xc(:,j));
end

% compare to the analytical form:

for j = 1:5
   for i = 1:201
       mu_x(i) = -1 + 0.01*(i-1);
       trm1 = (1 + mu_x(i)^2)/2;
       trm2 = 1/(1 + alpha(j)*(1 - mu_x(i)));
       trm3 = 2*alpha(j)^2*(1 - mu_x(i))^2*trm2/trm1;
       f_ana(i,j) = (ntrials/200)*trm1*trm2^2*(1 + trm3);
   end
end

plot(xout(:,1),nc(:,1),'rx',xout(:,2),nc(:,2),'b+',xout(:,3),nc(:,3),'ko',...
    xout(:,4),nc(:,4),'gx',xout(:,5),nc(:,5),'m+')
hold on
plot(mu_x,f_ana(:,1),'k',mu_x,f_ana(:,2),'c',mu_x,f_ana(:,3),'r',mu_x,f_ana(:,4),'b',...
    mu_x,f_ana(:,5),'k')
axis([-1 1 0 ntrials/200])
xlabel('\mu = cos(\theta)')
ylabel('Counts')

% end of program
