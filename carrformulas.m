function [Finf] = carrformulas(R,Rc,D,k,P,d,coating,system)

% The function "carrformulas" outputs the formulas for the total
% fraction of drug release based on Elliot Carr's 2024 paper:
% "Total fraction of drug released from diffusion-controlled delivery 
% systems with binding reactions".
%
% Inputs:
% (R) - device radius
% (Rc) - core-shell radius
% (D) - diffusivity
% (k) - reaction-rate
% (P) - mass transfer coefficient
% (d) - device geometry dimension
% (coating) - permeability type
% (system) - device system (monolithic or coreshell)
%
% Outputs:
% [Finf] - total fraction of drug released

% Non-dim variables
lambda = sqrt(k*R^2/D);
Pcal = P*R/D;
Rcalc = Rc/R;
Rcals = 1 - Rcalc;

if isequal(system,'monolithic')

    if isequal(coating,'fully-permeable')

        if isequal(d,1)
            Finf = sinh(lambda)/(lambda*cosh(lambda));

        elseif isequal(d,2)
            Finf = 2*besseli(1,lambda)/(lambda*besseli(0,lambda));

        elseif isequal(d,3)
            Finf = 3*(lambda*cosh(lambda)-sinh(lambda))/(lambda^2*sinh(lambda));

        end

    elseif isequal(coating,'semi-permeable')

        if isequal(d,1)
            Finf = Pcal*sinh(lambda)/(lambda*(lambda*sinh(lambda)+Pcal*cosh(lambda)));

        elseif isequal(d,2)
            Finf = 2*Pcal*besseli(1,lambda)/(lambda*(lambda*besseli(1,lambda)+Pcal*besseli(0,lambda)));

        elseif isequal(d,3)
            Finf = 3*Pcal*(lambda*cosh(lambda)-sinh(lambda))/(lambda^2*(lambda*cosh(lambda)+(Pcal-1)*sinh(lambda)));

        end

    end

elseif isequal(system,'coreshell')

    if isequal(coating,'fully-permeable')

        if isequal(d,1)
            Finf = 1/(cosh(lambda*Rcals));

        elseif isequal(d,2)
            Finf = 1/(Rcalc*lambda*(besselk(0,lambda)*besseli(1,Rcalc*lambda)+besseli(0,lambda)*besselk(1,Rcalc*lambda)));

        elseif isequal(d,3)
            Finf = lambda/(sinh(lambda*Rcals)+lambda*Rcalc*cosh(lambda*Rcals));

        end

    elseif isequal(coating,'semi-permeable')

        if isequal(d,1)
            Finf = Pcal/(Pcal*cosh(lambda*Rcals)+lambda*sinh(lambda*Rcals));

        elseif isequal(d,2)
            Finf = Pcal/(Rcalc*lambda*((Pcal*besselk(0,lambda)-lambda*besselk(1,lambda))*besseli(1,Rcalc*lambda)+(Pcal*besseli(0,lambda)+lambda*besseli(1,lambda))*besselk(1,Rcalc*lambda)));

        elseif isequal(d,3)
            Finf = Pcal*lambda/((Pcal-1+Rcalc*lambda^2)*sinh(lambda*Rcals)+lambda*(Rcals+Pcal*Rcalc)*cosh(lambda*Rcals));

        end

    end

end