function [W] = mosh(N, panickedFraction, roomSize, alpha, sigma)

%alpha = 0.35;
epsilon = 25;
my = 1;
%sigma = 1;

mass = 1;
v0 = 1;
r0 = 1;

dt = 0.1;
T = 20000;

L = roomSize;

p = L*rand(N, 2);
v = zeros(N, 2);

numPanicked = floor(N * panickedFraction);
isPanicked = zeros(N, 1);

isPanicked(1:numPanicked) = 1;

for t = 0:T
    % Update position vector according to boundary conditions
    p = p + dt * v;
    
    v(p > L) = 0; %-v(p > L);
    p(p > L) = L;
    v(p < 0) = 0; %-v(p < 0);
    p(p < 0) = 0;
    
    % Center of mass for active MASHers
    MC = sum(p(isPanicked == 1, :)) / numPanicked;
    
    F = zeros(size(v));
    for i = 1:N
        % Repulsive force
        r_i = [p(:, 1) - p(i, 1) p(:, 2) - p(i, 2)];%bsxfun(@minus, p, p(i, :));   % Relative position vector
        ri = sqrt(sum(r_i.^2, 2)); 
        
        flockIndex = find((ri < 4 * r0) & (ri > 1e-9));
        repIndex = flockIndex(ri(flockIndex) < 2 * r0);
        
        Frep = [0 0];
        if (~isempty(repIndex))
            Frep = bsxfun(@rdivide, r_i(repIndex, :), ri(repIndex));
            Frep = -epsilon / (2 * r0)^1.5 * ((2 * r0 - ri(repIndex)').^1.5) * Frep;
        end
        
        if (isPanicked(i))
            % Propulsion force for Active MASHer
            vi = sqrt(sum(v(i, :).^2));
            Fprop = my * (v0 - vi) * v(i, :) / (vi + (vi == 0));
            
            % Flocking force
            Fflock = [0 0];
            if (~isempty(flockIndex))
                Fflock = Fflock + sum(v(flockIndex, :));
                Fflock = alpha * Fflock / (sqrt(sum(Fflock.^2)) + isequal(Fflock, [0 0]));
            end
            
            % Noise
            theta = pi * rand();
            Fnoise = sigma * randn() * [cos(theta), sin(theta)];
            
            % Total force
            F(i, :) = Frep + Fprop + Fflock + Fnoise;
        else
            % Propulsion force for passive MASHer
            vi = sqrt(sum(v(i, :).^2));
            Fprop = my * (-vi) * v(i, :) / (vi + (vi == 0));
            
            % Total force
            F(i, :) = (Frep + Fprop);
        end
    end
    % Update velocity
    v = v + dt * F / mass;
    
    % Graphics
    
    if (~mod(t, 1))
        hold off
        plot(p(isPanicked == 0, 1) ,p(isPanicked == 0, 2), 'k.', 'markersize', 25);
        hold on
        plot(p(isPanicked == 1, 1) ,p(isPanicked == 1, 2), 'r.', 'markersize', 25);
        plot(MC(1), MC(2), 'b.', 'markersize', 25);
        axis([0 L 0 L]);
        quiver(p(isPanicked == 1, 1), p(isPanicked == 1, 2), v(isPanicked == 1, 1), v(isPanicked == 1, 2), 0);
        pause(0.01);
    end
end
W = 1;
return