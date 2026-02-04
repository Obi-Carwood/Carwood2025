function [c,t,F] = numerical(R,D,k,P,c0,T,Nr,Nt,itol,d,coating)

% The function "numerical" performs the numerical solve for the
% governing mechanistic model and determines the fraction of drug release.
%
% Inputs:
% (R) - device radius
% (D) - diffusivity function
% (k) - reaction-rate function
% (P) - mass transfer coefficient
% (c0) - initial concentration function
% (T) - endtime
% (Nr) - number of spatial nodes
% (Nt) - number of temperal nodes
% (itol) - integration absolute tolerance
% (d) - device geometry dimension
% (coating) - permeability type
%
% Outputs:
% [c] - output concentration (solution)
% [t] - corresponding time array
% [F] - fraction of drug release profile array

% Finite difference between nodes
h = R / (Nr-1);

% Create radius vector
r = linspace(0,R,Nr);

% West and east radii vectors
west = [r(1) (r(1:Nr-1)+r(2:Nr))/2];
east = [(r(1:Nr-1)+r(2:Nr))/2 r(Nr)];

% Finite volumes
V = (east.^d - west.^d)/d;

% Arrays for simplicity
FunW = D(west).*west.^(d-1);
FunE = D(east).*east.^(d-1);

% Initialise A
A = zeros(Nr,Nr);

% Boundary node i = 1
A(1,1) = -FunE(1)/(V(1)*h) - k(r(1));
A(1,2) = FunE(1)/(V(1)*h);

% Interior nodes i = 2, ..., N-1
for i = 2:Nr-1
    A(i,i-1) = FunW(i)/(V(i)*h);
    A(i,i) = -(FunW(i) + FunE(i))/(V(i)*h) - k(r(i));
    A(i,i+1) = FunE(i)/(V(i)*h);
end

% Boundary node i = Nr
if isequal(coating,'fully-permeable')
    A(end,:) = [];
    A(:,end) = [];
elseif isequal(coating,'semi-permeable')
    A(Nr,Nr-1) = FunW(Nr)/(V(Nr)*h);
    A(Nr,Nr) = -FunW(Nr)/(V(Nr)*h) - k(r(i)) - P*R^(d-1)/V(Nr);
end

% Finite time difference
dt = T / (Nt-1);

% Time array
t = linspace(0,T,Nt);

% Model Size
[Nm,~] = size(A);

% Adjusted radius vector
if isequal(coating,'fully-permeable')
    rm = r(1:end-1);
else
    rm = r;
end

% Initialise Solution
c = zeros(Nt,Nm);
c(1,:) = c0(rm);

% Solve
At = eye(Nm) - dt*A;
SAt = sparse(At);
for n = 1:Nt-1
    c(n+1,:) = SAt \ c(n,:)';
end

% Implement boundary solutions
if isequal(coating, 'fully-permeable')
    c = [c [c0(R); zeros(Nt-1,1)]];
end

% Concentration gradient at r=R (6th order)
grad = -D(R) * (147*c(1:Nt,Nr) - 360*c(1:Nt,Nr-1) + 450*c(1:Nt,Nr-2) - ...
         400*c(1:Nt,Nr-3) + 225*c(1:Nt,Nr-4) - 72*c(1:Nt,Nr-5) + 10*c(1:Nt,Nr-6)) / (60*h);

% Constant term
C = R^(d-1)/integral(@(r) r.^(d-1) .* c0(r),0,R,"AbsTol",itol);

% Compute integral
F = zeros(size(t));
for ti = 2:length(t)
    F(ti) = C * trapz(t(1:ti), grad(1:ti));
end

end