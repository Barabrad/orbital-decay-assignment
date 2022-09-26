% ASTE 331a, Orbital Motion
% Brad Barakat
% 22 September 2022

% The script was put in a function to allow for potential scalability
function bbOrbitMotion()

    year0 = year(date);
    
    % If the data was in a file, a for loop would've worked here

    % LEO Orbit Data
    j = 1;
    orbit(j).p0 = [6778, 0, 0]; % Cartesian position (km)
    orbit(j).v0 = [0, 7.66855, 0]; % Cartesian velocity (km/s)
    orbit(j).T = 5554.46; % period (s)
    orbit(j).n_max = 5; % max number of orbits
    orbit(j).BC = 1; % ballistic coefficient (kg/m^2)
    orbit(j).rf = 6478; % final radius (km)
    orbit(j).type = "LEO";

    % Elliptic Orbit Data
    j = 2;
    orbit(j).p0 = [0, -3388.8918, -5869.7327];
    orbit(j).v0 = [8.90236, 0, 0];
    orbit(j).T = 29806.9016;
    orbit(j).n_max = 100;
    orbit(j).BC = 1;
    orbit(j).rf = 6478;
    orbit(j).type = "Elliptic";

    % Main loop
    for i = 1:j
        figure
        sgtitle("Orbit #" + i + " (" + orbit(i).type + ")")
        % solve ODE for each orbit, with and without drag
        for hasDrag = 0:1
            orbit(i).n_max = 1000*hasDrag + orbit(i).n_max;
            [decayTime, times, radii] = orbitaldecay(orbit(i), hasDrag, ...
                year0);
            if (hasDrag == 0)
                graphWord = "No";
            else
                graphWord = "With";
                disp("The decay time is " + decayTime/(60*60*24) + ...
                    " days for an orbit (with drag) characterized by ")
                disp(orbit(i))
            end

            if (orbit(i).type == "LEO")
                % 4-tile subplot
                subplot(2,2,(2*hasDrag)+1)
                plot(radii(:,1), radii(:,2))
                title(newline + "x-y Projection (" + graphWord + " Drag)")
                xlabel("x-coordinate (km)")
                ylabel("y-coordinate (km)")

                subplot(2,2,(2*hasDrag)+2)
                r = sqrt(radii(:,1).^2 + radii(:,2).^2 + radii(:,3).^2);
                plot(times, r)
                title("Orbital Radius (" + graphWord + " Drag)")
                xlabel("Time (s)")
                ylabel("Orbital Radius (km)")
            else
                % 2-tile subplot
                subplot(2,1,hasDrag+1)
                plot3(radii(:,1), radii(:,2), radii(:,3))
                title("3D Orbit (" + graphWord + " Drag)")
                xlabel("x-coordinate (km)")
                ylabel("y-coordinate (km)")
                zlabel("z-coordinate (km)")
            end
        end
    end

end

% Function that calculates the decay
function [decayTime, times, radii] = orbitaldecay(orbit, hasDrag, year0)
    % Constants
    R0 = 6378; % Earth's radius (km)
    mu = 398500; % Earth's gravitational constant (km^3/s^2)
    secPerYear = 60*60*24*365.25;
    
    % Calculations
    t_max = orbit.T*orbit.n_max;
    x0 = [orbit.p0, orbit.v0]';
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Events', ...
        @(t,x)stopAtFinal(t, x, orbit.rf));
    [times, states] = ode45(@(t,x)motionODEs(t, x, orbit.BC, R0, mu, ...
        hasDrag, secPerYear, year0), [0, t_max], x0, options);
    decayTime = times(end);
    radii = [states(:,1), states(:,2), states(:,3)];
end

% ODE function
function dxdt = motionODEs(t, x, BC, R0, mu, hasDrag, secPerYear, year0)
    % x = [px, py, pz, vx, vy, vz]'
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    if (r <= R0)
        r = R0;
    end
    acc_g = -mu/(r^3);

    dxdt = zeros(size(x));
    dxdt(1) = x(4);
    dxdt(2) = x(5);
    dxdt(3) = x(6);
    dxdt(4) = acc_g*x(1);
    dxdt(5) = acc_g*x(2);
    dxdt(6) = acc_g*x(3);

    if (hasDrag)
        rho = atmosphericDensity(r-R0, year0+t/secPerYear);
        v = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
        acc_drag = -(rho*v^2/(2*BC))*1000; % 1000 is a conversion factor
        dxdt(4) = dxdt(4) + acc_drag*x(4)/v;
        dxdt(5) = dxdt(5) + acc_drag*x(5)/v;
        dxdt(6) = dxdt(6) + acc_drag*x(6)/v;
    end
end

% ODE event
function [position,isterminal,direction] = stopAtFinal(~, x, rf)
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    position = r-rf; % The value that we want to be zero
    isterminal = 1; % Halt integration 
    direction = 0; % The zero can be approached from either direction
end
