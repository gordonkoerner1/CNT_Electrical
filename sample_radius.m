%% This samples ro from a normal and ri from a log-normal distribution
%
function [ro,ri]=sample_radius(numberBeams,mu_ro,sigma_ro,shape_ri,scale_ri)

% mu_ro = 12.34;   % Corrected mean value
% sigma_ro = 3.46;  % Corrected standard deviation
% 
% shape_ri = 0.25;  % Corrected shape value
% scale_ri = 6.37;  % Corrected scale value

% Define the range for t = ro - ri
t_min = 0.90;  % Minimum allowed value for t
t_max = 12.2;  % Maximum allowed value for t

% Initialize storage for valid ro and ri samples
ro = zeros(numberBeams,1,'double');
ri = zeros(numberBeams,1,'double');

% Loop until we collect enough valid samples
count=1;
while count < numberBeams + 1
    % Sample ro and ri
    ro_sample = normrnd(mu_ro, sigma_ro, 1, 1);
    ri_sample = lognrnd(log(scale_ri), shape_ri, 1, 1);

    % Compute t = ro - ri
    t = ro_sample - ri_sample;

    % Check if t is within the desired range
    if t >= t_min && t <= t_max
        ro(count) = ro_sample*1e-9;
        ri(count) = ri_sample*1e-9;
        count=count+1;
    end
end
