function [Finf] = analytical(D,k,c0,P,Ns,d,coating,tolerances)

% The function "analytical" computes the total fraction of drug-release
% based on the analytical expressions.
%
% Inputs:
% (D) - diffusivity function
% (k) - reaction-rate function
% (c0) - initial concentration function
% (P) - mass transfer coefficient
% (Ns) - solution truncation term
% (d) - device geometry dimension
% (coating) - permeability type
% (tolerances) - [integral tol,bimethod maxiters,bimethod tol]
%
% Output:
% [Finf] - total fraction of drug release

% Extract tolerance information
atol = tolerances(1);
rtol = tolerances(2);
maxiters = tolerances(3);
btol = tolerances(4);

% Determine Eigenvalues
lams = eigvals(Ns,D,P,d,coating,maxiters,btol);

% Create Eigenfunction and Eigenfunction Derivative
Xn = @(r,n) eigfuncXn(r,lams(n),d);
Xdn = @(r,n) eigfuncXdn(r,lams(n),d);

% Derivative of parameter functions
r = sym('r');
Dd = diff(D(r),r); Dd = matlabFunction(Dd,'Vars',r);
kd = diff(k(r),r); kd = matlabFunction(kd,'Vars',r);

% Check for zero derivative
dzero = false;
x = sym('x');
if isequal(Dd(x),0) && isequal(kd(x),0)
    D = D(1);
    k = k(1);
    dzero = true;
end

% Initialise A and b matrices
A = zeros(Ns);
b = zeros(Ns,1);

% Create A and b by looping over rows and columns
warning off
for m = 1:Ns
    for n = 1:Ns
        if dzero
            A(m,n) = -(lams(n)^2*D+k)*integral(@(r) r.^(d-1)...
                .*Xn(r,n).*Xn(r,m),0,1,"AbsTol",atol,"RelTol",rtol);
        else
            A(m,n) = integral(@(r) r.^(d-1).*(Dd(r).*Xdn(r,n).*Xn(r,m)...
                -(lams(n)^2*D(r)+k(r)).*Xn(r,n).*Xn(r,m)),0,1,"AbsTol",atol,"RelTol",rtol);
        end
    end
    b(m) = -integral(@(r) r.^(d-1).*c0(r).*Xn(r,m),0,1,"AbsTol",atol,"RelTol",rtol);
end

% Solve for alpha terms
alpha = A\b;

% Integral term from initial concentration
c0Int = integral(@(r) r.^(d-1).*c0(r),0,1,"AbsTol",atol,"RelTol",rtol);

% Compute sum for Finf
S = 0;
for n = 1:Ns

    % Solve for given dimension input
    if isequal(d,1)
        A = 2*alpha(n)*sqrt(lams(n)) / (c0Int*sqrt(sin(2*lams(n)) + 2*lams(n)));
        if dzero
            S = S + A*k*integral(@(r) cos(lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        else
            S = S + A*integral(@(r) k(r).*cos(lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        end

    elseif isequal(d,2)
        A = sqrt(2)*alpha(n) / (c0Int*sqrt(besselj(0,lams(n))^2 + besselj(1,lams(n))^2));
        if dzero
            S = S + A*k*integral(@(r) r.* besselj(0,lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        else
            S = S + A*integral(@(r) r.* k(r) .* besselj(0,lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        end

    elseif isequal(d,3)
        A = 2*alpha(n)*sqrt(lams(n)) / (c0Int*sqrt(2*lams(n) - sin(2*lams(n))));
        if dzero
            S = S + A*k*integral(@(r) r.* sin(lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        else
            S = S + A*integral(@(r) r.* k(r) .* sin(lams(n)*r),0,1,"AbsTol",atol,"RelTol",rtol);
        end

    end
end

warning on

% Compute Finf
Finf = 1 - S;

end



%% Compute Eigenvalues Function
function [lams] = eigvals(Ns,D,P,d,coating,maxiters,tol)

% The function "eigvals" computes the first Ns eigenvalues for the model
% provided the dimension and coating.
%
% Inputs:
% (Ns) - number of eigenvalues to compute
% (D) - diffusivity function
% (P) - mass-transfer coefficient
% (d) - device geometry dimension
% (coating) - permeability type
% (maxiters) - maximum iterations
% (tol) - tolerance
%
% Output:
% [lams] - vector array of eigenvalues

if isequal(coating,'fully-permeable')
    if isequal(d,1) % d = 1 | Fully-Permeable
        lams = pi*((1:Ns) - 1/2);

    elseif isequal(d,2) % d = 2 | Fully-Permeable
        f = @(x) besselj(0,x);  % root function
        lams = zeros(1,Ns);     % initialise

        % Apply bisection method
        for n = 1:Ns
            a0 = (n-1)*pi;  % left bound
            a1 = n*pi;      % right bound
            lams(n) = bimethod(f,a0,a1,maxiters,tol);
        end

    elseif isequal(d,3) % d = 3 | Fully-Permeable
        lams = pi*(1:Ns);

    end

elseif isequal(coating,'semi-permeable')
    if isequal(d,1) % d = 1 | Semi-Permeable
        f = @(x) D(1)*x.*sin(x)-P*cos(x);   % root function
        lams = zeros(1,Ns);                 % initialise

        % Apply bisection method (first)
        a0 = 0;     % left bound
        a1 = pi/2;  % right bound
        lams(1) = bimethod(f,a0,a1,maxiters,tol);

        % Apply bisection method (rest)
        for n = 2:Ns
            a0 = (n-3/2)*pi;    % left bound
            a1 = (n-1/2)*pi;    % right bound
            lams(n) = bimethod(f,a0,a1,maxiters,tol);
        end

    elseif isequal(d,2) % d = 2 | Semi-Permeable
        f = @(x) D(1)*x.*besselj(1,x)-P*besselj(0,x);   % root function
        lams = zeros(1,Ns);                             % initialise

        % Apply bisection method
        for n = 1:Ns
            a0 = (n-1)*pi;  % left bound
            a1 = n*pi;      % right bound
            lams(n) = bimethod(f,a0,a1,maxiters,tol);
        end

    elseif isequal(d,3) % d = 3 | Semi-Permeable
        f = @(x) (D(1)-P)*sin(x)-D(1)*x.*cos(x);    % root function
        lams = zeros(1,Ns);                         % initialise

        % Apply bisection method
        for n = 1:Ns
            a0 = (n-1/2)*pi;    % left bound
            a1 = (n+1/2)*pi;    % right bound
            lams(n) = bimethod(f,a0,a1,maxiters,tol);
        end

    end
end
end

%% Bisection Method Function
function [root] = bimethod(f,a0,a1,maxiters,tol)

% The function "bimethod" computes the root for the given
% function (f) by implementing the bisection method for bounds
% over (a0,a1).
%
% Inputs:
% (f) - function to find roots
% (a0,a1) - bounds to find root within
% (maxiters) - maximum iterations
% (tol) - absolute tolerance between bounds
%
% Output:
% [root] - approximate root to function (f)

i = 0;  % initialise iterations
while i < maxiters && abs(a1-a0) > tol

    root = (a0+a1)/2;   % find midpoint

    % update bounds
    if sign(f(root)) == sign(f(a0))
        a0 = root;
    else
        a1 = root;
    end

    % iterate
    i = i+1;

end
end

%% Eigenfunction Function
function [Xn] = eigfuncXn(r,lambda,d)

% The function "eigfuncXn" gives the eigenfunction for the
% analytical method developed from the sturm-louiville
% problem.
%
% Inputs:
% (r) - independent variable
% (lambda) - corresponding eigenvalue
% (d) - device geometry dimension
%
% Outputs:
% [Xn] - eigenfunction

% Compute function from given dimension and lambda input
if isequal(d,1)
    X0 = 2*sqrt(lambda) / sqrt(sin(2*lambda) + 2*lambda);
    Xn = X0*cos(lambda*r);

elseif isequal(d,2)
    X0 = sqrt(2) / sqrt(besselj(0,lambda)^2 + besselj(1,lambda)^2);
    Xn = X0*besselj(0,lambda*r);

elseif isequal(d,3)
    X0 = 2*sqrt(lambda) / sqrt(2*lambda - sin(2*lambda));
    Xn = X0*sin(lambda*r)./r;

end
end

%% Eigenfunction Derivative Function
function [Xdn] = eigfuncXdn(r,lambda,d)

% The function (eigfuncXdn) gives the derivative of the
% eigenfunction for the analytical method developed from the
% sturm-louiville problem.
%
% Inputs:
% (r) - independent variable
% (lambda) - corresponding eigenvalue
% (d) - device geometry dimension
%
% Output:
% [Xdn] - eigenfunction derivative

% Compute function from given dimension and lambda input
if isequal(d,1)
    X0 = 2*sqrt(lambda) / sqrt(sin(2*lambda) + 2*lambda);
    Xdn = -lambda*X0*sin(lambda*r);

elseif isequal(d,2)
    X0 = sqrt(2) / sqrt(besselj(0,lambda)^2 + besselj(1,lambda)^2);
    Xdn = -lambda*X0*besselj(1, lambda*r);

elseif isequal(d,3)
    X0 = 2*sqrt(lambda) / sqrt(2*lambda - sin(2*lambda));
    Xdn = X0 * (lambda*cos(lambda*r)./r - sin(lambda*r)./r.^2);

end
end