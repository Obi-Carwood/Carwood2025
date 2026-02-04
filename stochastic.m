function [Fdist] = stochastic(R,D,k,c0,P,T,Nts,Nl,Np,trials,d,coating)

% The function "Stochastic" performs the stochastic process and simulation
% for the governing mechanistic model.
%
% Inputs:
% (R) - device radius
% (D) - diffusivity function
% (k) - reaction-rate function
% (c0) - initial concentration function
% (P) - mass transfer coefficient
% (T) - runtime
% (Nts) - number of time-steps
% (Nl) - number of lattice sites
% (Np) - number of particles
% (trials) - number of simulation trials
% (d) - device geometry dimension
% (coating) - permeability type
%
% Output:
% [Fdist] - distribution array of the stochastic release profile

% Radius vector
r = linspace(0,R,Nl);

% West and east vectors
west = [r(1) (r(1:Nl-1)+r(2:Nl))/2];
east = [(r(1:Nl-1)+r(2:Nl))/2 r(Nl)];

% Finite difference
h = R / (Nl-1);

% Finite volumes
V = (east.^d - west.^d)/d;

% Functions for simplicty
FunW = D(west).*west.^(d-1);
FunE = D(east).*east.^(d-1);

% Conditions on time-step
cond = [V(1)*h/FunE(1), V(2:Nl-1)*h./(FunW(2:Nl-1)+FunE(2:Nl-1)), 1./k(r(1:Nl-1))];
if isequal(coating,'semi-permeable')
    cond = [cond, V(Nl)*h/FunW(Nl), 1/k(r(Nl)), V(Nl)/(P*R^(d-1))];
end

% Check conditions
maxtau = min(cond);
if Nts < ceil(T/maxtau)
    error(['Nts must be greater than or equal to ', num2str(ceil(T/maxtau),'%i')])
end

% Set time-step
tau = T/Nts;

% Create matrix As
As = zeros(Nl,Nl);

% Node i = 1
As(1,1) = -FunE(1)/(V(1)*h);
As(1,2) = FunE(1)/(V(1)*h);

% Nodes i = 2,...,Nl-1
for i = 2:Nl-1
    As(i,i-1) = FunW(i)/(V(i)*h);
    As(i,i) = -(FunW(i) + FunE(i))/(V(i)*h);
    As(i,i+1) = FunE(i)/(V(i)*h);
end

% Node i = Nl
As(Nl,Nl-1) = FunW(Nl)/(V(Nl)*h);
As(Nl,Nl) = -FunW(Nl)/(V(Nl)*h);

% Create diagonal matrix of V and its inverse
Vol = diag(V);
Volinv = diag(1./V);

% Create Markov Matrix Ps
Ps = Vol * (eye(Nl) + tau*As) * Volinv;
Ps = Ps';

% Create Markov Matrix Pr
if isequal(coating,'semi-permeable')
    Pr = tau*P*R^(d-1)/V(Nl);
end

% Create Markov Matrix Pb
Pb = tau*k(r);

% Initial particle distribution
Xc0 = c0(r);    % initial concentration distribution
X0 = round(Xc0.*V / (Xc0*V'/Np)); % initial particle distribution
Np = sum(X0);   % updated number of particles
parray = 1:Np;

% Arrange particle distribution by particle number
Xp0 = zeros(1,Np);
n = 1;
for i = 1:Nl
    npl = X0(i);
    for p = 1:npl
        Xp0(n) = i;
        n = n+1;
    end
end

% Initialise distribution array
Fdist = zeros(trials,Nts);

% Loop over trials
for trial = 1:trials

    % Initialise Information
    Xp = Xp0;
    active = ones(1,Np);
    release = zeros(1,Np);

    % Increment time
    for t = 2:Nts

        % Loop over particles
        for p = parray

            % Check for active particle
            if isequal(active(p),0)
                continue
            end

            % Extract lattice site number
            i = Xp(p);

            % Release event
            if isequal(i,Nl)
                if isequal(coating,'fully-permeable')
                    active(p) = 0;
                    release(p) = 1;
                    continue
                elseif isequal(coating,'semi-permeable')
                    rn = rand;
                    if rn < Pr
                        active(p) = 0;
                        release(p) = 1;
                        continue
                    end
                end
            end

            % Binding event
            rn = rand;
            if rn < Pb(i)
                active(p) = 0;
                continue
            end

            % Movement event
            rn = rand;
            if isequal(i,1)
                if rn < Ps(i,i)
                    j = i;
                else
                    j = i+1;
                end
            elseif isequal(i,Nl)
                if rn < Ps(i,i)
                    j = i;
                else
                    j = i-1;
                end
            else
                if rn < Ps(i,i-1)
                    j = i-1;
                elseif rn < Ps(i,i-1) + Ps(i,i)
                    j = i;
                else
                    j = i+1;
                end
            end
            Xp(p) = j;
        end

        % Record fraction array at this time
        fraction = sum(release) / Np;
        Fdist(trial,t) = fraction;

        % Check for active particles
        if isequal(sum(active),0)
            Fdist(trial,t:end) = sum(release) / Np;
            break
        end

    end

end