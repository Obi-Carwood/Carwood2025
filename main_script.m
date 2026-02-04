%% MPhil - Obi Carwood %%

clear, clc, close all

% ======================================================================= %
%% Figures and Results Settings %%
% ======================================================================= %

% == Enable Plots and Results == %
enableN = true;      % enable numerical results and plots
enableS = true;      % enable stochastic results and plots
enableA = true;      % enable analytical results and plots
enableG = false;     % enable general profile plots
enableP = false;     % enable parameter functions plots

% == Saving Figures == %
savefigRP = false;   % save release profile figure
pathRP = cd;         % release profile save path
savefigPF = false;   % save parameter functions figure
pathPF = cd;         % parameter functions save path
savefigGP = false;   % save general profile figure
pathGP = cd;

% == Compare to Carr (2024) == %
carrmono = false;    % compare to (Carr,2024) monolithic formulas
carrcore = false;    % compare to (Carr,2024) coreshell formulas

% == Dimension Set == %
dimSet = [1 2 3];    % dimensions to plot over (1-slab 2-cylinder 3-sphere)

% ======================================================================= %
% == Plotting Specifications == %
% ======================================================================= %

% Combined
fsz = 20;   % fontsize

% Numerical
Nlw = 3;                    % linewidth
Ncolourd1 = [0 130 200];    % colour d=1
Ncolourd2 = [0 122 67];     % colour d=2
Ncolourd3 = [203 42 42];    % colour d=3
Ncolour = [Ncolourd1;Ncolourd2;Ncolourd3]/255;

% Stochastic
Slw1 = 1;               % linewidth (upper/lower)
Slw2 = 0.5;             % linewidth (vertical lines)
redfac = 50;            % reduction factor (plotting)
redfac2 = 350;          % further reduction factor (plotting vertical lines)
op = 1;                 % opacity
bg = [255 255 255];     % background colour
Scolourd1 = op*Ncolourd1 + (1-op)*bg;   % colour d=1
Scolourd2 = op*Ncolourd2 + (1-op)*bg;   % colour d=2
Scolourd3 = op*Ncolourd3 + (1-op)*bg;   % colour d=3
Scolour = [Scolourd1;Scolourd2;Scolourd3]/255;

% Analytical
Alw = 3;                % linewidth
Acolourd1 = [0 0 0];    % colour d=1
Acolourd2 = [0 0 0];    % colour d=2
Acolourd3 = [0 0 0];    % colour d=3
Acolour = [Acolourd1;Acolourd2;Acolourd3]/255;

% Parameter Functions
lfsz = 25;                  % label fontsize
Plw = 3;                    % linewidth
Pcolour1 = [76 47 202];     % colour 1
Pcolour2 = [23 195 178];    % colour 2
Pcolour3 = [56 168 94];     % colour 3
Pcolour4 = [194 172 30];    % colour 4
Pcolour = [Pcolour1; Pcolour2; Pcolour3; Pcolour4]/255;

% General Profiles
capfsz = 18;                % caption fontsize
labfsz = 18;                % label fontsize
Glw = 4;                    % linewidth
Gcolour1 = [0 0 0];         % colour 1
Gcolour2 = [255 0 0];       % colour 2
Gcolour = [Gcolour1; Gcolour2]/255;

% ======================================================================= %
%% Input Parameters %%
% ======================================================================= %

% == Permeability == %
% coating = 'fully-permeable';
coating = 'semi-permeable';

% == Model Parameters == %
R = 1e-4;       % device radius
P = 5e-8;       % mass transfer coefficient
T = 1e5;        % end-time

% == Numerical Parameters == %
Nr = 5001;      % number of spatial nodes
Nt = 1e4;       % number of temporal nodes
itolN = 1e-12;  % integration tolerance

% == Analytical Constraints == %
Ns = 400;       % solution truncation term
atolA = 1e-12;  % integration absolute tolerance
rtolA = 1e-12;  % integration relative tolerance
maxiters = 1000;% root finding maximum iterations
btol = 1e-12;   % root finding tolerance

% == Stochastic Constraints == %
Nts = 1.5e6;    % number of time-steps
Np = 200;       % number of particles
Nl = 51;        % number of lattice sites
uP = 0.9;       % upper stochastic percentile
lP = 0.1;       % lower stochastic percentile
trials = 100;   % number of stochastic trials

% ======================================================================= %
% == Parameter Functions == %
% ======================================================================= %

% == Initial Concentraton == %
C0 = 0.4;
c0 = @(r) C0*ones(size(r));

% Alpha value
alphaSet = [1e-4 20 80 1e4];
alpha = alphaSet(1);

% Diffusivity Bounds
Dmin = 1e-13;
Dmax = 1e-11;

% Reaction-rate bounds (kmin will change to 0 if comparing to (Carr,2024) coreshell)
kmin = 8e-5;
kmax = 1e-4;

% Custom Parameter Functions
enableC = false;                % enable custom functions
D = @(r) 1e-13*ones(size(r));   % diffusivity
k = @(r) 1e-4*ones(size(r));    % reaction-rate

% Extra
Rc = R/2;   % core radius (coreshell system)

% ======================================================================= %
%% Results %%
% ======================================================================= %

% Pre-allocate for storing data
NUMERICAL = zeros(2,Nt,3);
ANALYTICAL = zeros(1,3);
STOCHASTIC = zeros(3,ceil(Nts/redfac),3);

% Loop over dimension set
firstloop = true;   % loop count
for d = dimSet

    % == Create Parameter Functions == %

    % Smooth Step-Wise Parameter Functions
    if ~enableC

        % Average Values
        Davg = 1/2^d*Dmax + (1-1/2^d)*Dmin;
        kavg = 1/2^d*kmin + (1-1/2^d)*kmax;

        % Determine Sigma
        D = @(sigma,r) Dmax + (Dmin-Dmax)*(1/2 + 1/pi*atan(alpha*(r-sigma)/R));
        sigma = fzero(@(sigma) d/R^d * integral(@(r) r.^(d-1) .* D(sigma,r),0,R) - Davg,R/2);

        % Parameter Functions
        D = @(r) Dmax + (Dmin-Dmax) * (1/2 + 1/pi * atan(alpha*(r - sigma)/R));
        k = @(r) kmin + (kmax-kmin) * (1/2 + 1/pi * atan(alpha*(r - sigma)/R));

        % Small Alpha or Monolithic
        if alpha < 1e-2 || carrmono
            D = @(r) Davg*ones(size(r));
            k = @(r) kavg*ones(size(r));
            c0 = @(r) C0*ones(size(r));
        end
    end
    
    % Coreshell
    if carrcore
        kmin = 0;
        D = @(r) Dmax + (Dmin-Dmax) * (1/2 + 1/pi * atan(alpha*(r - Rc)/R));
        k = @(r) kmin + (kmax-kmin) * (1/2 + 1/pi * atan(alpha*(r - Rc)/R));
        c0 = @(r) C0*(1-heaviside(r-Rc));
    end

    % == Non-Dimensionalise == %

    % Dimensionalised Parameters
    rn = R;                             % radius
    Dn = max(D(linspace(0,R,Nr)));      % diffusivity
    tn = rn^2/Dn;                       % time
    cn = max(c0(linspace(0,R,Nr)));     % concentration
    kn = Dn/rn^2;                       % reaction-rate
    c0n = cn*tn*Dn/rn^2;                % initial concentration
    Pn = Dn/rn;                         % mass transfer coefficient

    % Non-Dimensionalised Parameters
    Dh = @(rh) D(rh*rn) / Dn;       % diffusivity
    kh = @(rh) k(rh*rn) / kn;       % reaction-rate
    c0hN = @(rh) c0(rh*rn) / cn;    % initial concentration Numerical
    c0hA = @(rh) c0(rh*rn) / c0n;   % initial concentration Analytical
    Ph = P / Pn;                    % mass transfer coefficient
    Rh = R / rn;                    % radius
    Th = T / tn;                    % time

    % == Results == %

    % Stochastic
    if enableS
        Fdist = stochastic(R,D,k,c0,P,T,Nts,Nl,Np,trials,d,coating);
    end

    % Numerical
    if enableN
        [ch,th,FN] = numerical(Rh,Dh,kh,Ph,c0hN,Th,Nr,Nt,itolN,d,coating);

        % Dimensionalise Numerical Results
        c = ch*cn;
        t = th*tn;
    end

    % Analytical
    if enableA
        [FinfA] = analytical(Dh,kh,c0hA,Ph,Ns,d,coating,[atolA rtolA maxiters btol]);
    end

    % Carr's Formula
    if enableC
        Dconst = D(1);
        kconst = k(1);
    else
        Dconst = Davg;
        kconst = kavg;
    end
    if carrmono
        system = 'monolithic';
        FinfC = carrformulas(R,Rc,Dconst,kconst,P,d,coating,system);
    elseif carrcore
        system = 'coreshell';
        Dshell = Dmin;
        kshell = kmax;
        FinfC = carrformulas(R,Rc,Dshell,kshell,P,d,coating,system);
    end

    % == Display Information == %
    if firstloop
        disp(['alpha = ', num2str(alpha), ', ', coating])
        disp('----------------------------')
        firstloop = false;
    end
    disp(['d = ', num2str(d), ':'])
    if enableA
        disp(['Finf = ', num2str(FinfA)])
    end
    if carrmono || carrcore
        disp(['Finf (Carr) = ', num2str(FinfC), ' (',system,')'])
        disp(['Abs Error (with Carr) = ', num2str(abs(FinfA-FinfC))])
    elseif enableN
        disp(['F(T) = ', num2str(FN(end))])
    end
    if enableN && enableA
        disp(['Abs Error = ', num2str(abs(FinfA-FN(end)))])
    end
    disp('----------------------------')

    % == Store Information == %

    % Stochastic
    if enableS
        U = quantile(Fdist,uP,1);
        L = quantile(Fdist,lP,1);
        tS = 0:length(Fdist(1,:))-1;
        tau = T/length(tS);
        tS = tau*tS;
        STOCHASTIC(1,:,d) = U(1:redfac:end);
        STOCHASTIC(2,:,d) = L(1:redfac:end);
        STOCHASTIC(3,:,d) = tS(1:redfac:end);
    end

    % Numerical
    if enableN
        NUMERICAL(1,:,d) = FN;
        NUMERICAL(2,:,d) = t;
    end

    % Analytical
    if enableA
        ANALYTICAL(d) = FinfA;
    end

end

% ======================================================================= %
%% Figures %%
% ======================================================================= %

figcount = 0;

% Release Profiles
if enableN || enableA || enableS
    figcount = figcount+1;
    figure(figcount)
    hold on
    box on

    for d = sort(dimSet,'descend')

        % Stochastic
        if enableS
            U = STOCHASTIC(1,:,d);
            L = STOCHASTIC(2,:,d);
            t = STOCHASTIC(3,:,d);
            plot(t,U,LineWidth=Slw1,Color=Scolour(d,:))
            plot(t,L,LineWidth=Slw1,Color=Scolour(d,:))
            U = U(1:redfac2:end);
            L = L(1:redfac2:end);
            t = t(1:redfac2:end);
            plot([t; t], [U; L], LineWidth=Slw2, Color=Scolour(d,:))
        end

    end

    for d = sort(dimSet,'descend')

        % Numerical
        if enableN
            F = NUMERICAL(1,:,d);
            t = NUMERICAL(2,:,d);
            plot(t,F,'LineWidth',Nlw,'Color',Ncolour(d,:))
        end

        % Analytical
        if enableA
            tA = [0 T];
            Finf = ANALYTICAL(d);
            plot(tA,Finf*ones(size(tA)),'--','LineWidth',Alw,'Color',Acolour(d,:))
            xshift = T/11;
            text(-xshift, Finf, '{\boldmath$F_{\infty}$}', 'interpreter','latex',...
                'FontSize',18,'Color',Ncolour(d,:),'FontWeight','bold')
        end

    end

    % Axes
    xlim([0 T]), set(gca,'XTick',[0 T])
    set(gca,'XTickLabel',{'0','$T$'},'TickLabelInterpreter','latex','fontsize',fsz)
    ylim([0 1]), set(gca,'YTick',[0 1])
    set(gca,'YTickLabel',{'0','1'},'TickLabelInterpreter','latex','fontsize',fsz)
    xl = xlabel('Time [$t$]',Interpreter='latex',FontSize=labfsz);
    xl.Position(2) = -0.05;
    yl = ylabel('Fraction of Drug Released [$F(t)$]',Interpreter='latex',FontSize=labfsz);
    yl.Position(1) = yl.Position(1)-T/30;

    % Title (inside plot)
    if isequal(coating,'fully-permeable')
        perm = 'Fully-Permeable';
    elseif isequal(coating,'semi-permeable')
        perm = 'Semi-Permeable';
    end
    if ~enableC
        text(T/5,0.05,['$\alpha =', num2str(alpha), '$ (',perm,')'],Interpreter="latex",FontSize=fsz)
    else
        text(T/2,0.05,['(',perm,')'],Interpreter="latex",FontSize=fsz)
    end


    % Save Figure
    if savefigRP
        filename = [coating,num2str(alpha)];
        print(gcf,fullfile(pathRP,filename),'-depsc2')
    end
end

% Parameter Profiles
if enableP

    for d = sort(dimSet,'descend')

        % Open Diffusivity Figure
        figcount = figcount+1;
        df = figure(figcount);
        clf(df)
        hold on

        % Open Reaction-Rate Figure
        figcount = figcount+1;
        rf = figure(figcount);
        clf(rf)
        hold on

        inx = 0;
        for alpha = sort(alphaSet,'descend')

            % Average Values
            Davg = 1/2^d*Dmax + (1-1/2^d)*Dmin;
            kavg = 1/2^d*kmin + (1-1/2^d)*kmax;

            % Determine Sigma
            D = @(sigma,r) Dmax + (Dmin-Dmax)*(1/2 + 1/pi*atan(alpha*(r-sigma)/R));
            sigma = fzero(@(sigma) d/R^d * integral(@(r) r.^(d-1) .* D(sigma,r),0,R) - Davg,R/2);

            % Parameter Functions
            D = @(r) Dmax + (Dmin-Dmax) * (1/2 + 1/pi * atan(alpha*(r - sigma)/R));
            k = @(r) kmin + (kmax-kmin) * (1/2 + 1/pi * atan(alpha*(r - sigma)/R));

            inx = inx+1;

            % == Diffusivity Plot == %

            figure(df)
            box on

            % Plot Functions
            fplot(D,[0 R],Color=Pcolour(inx,:),LineWidth=Plw)

            % Axes
            xlim([0 R]), set(gca,'XTick',[0 R/2 R])
            set(gca,'XTickLabel',{'0','$R/2$','$R$'},'TickLabelInterpreter','latex','fontsize',fsz)
            rnge = abs(Dmax - Dmin); mgin = 0.05;
            ylim([Dmin-mgin*rnge Dmax+mgin*rnge]), set(gca,'YTick',[Dmin Davg Dmax])
            set(gca,'YTickLabel',{'$D_{\rm{min}}$','$D_{\rm{avg}}$','$D_{\rm{max}}$'},...
                'TickLabelInterpreter','latex','fontsize',fsz)
            xlabel('$r$',Interpreter='latex',FontSize=lfsz)
            ylabel('$D(r)$',Interpreter='latex',FontSize=lfsz)

            % == Reaction-Rate Plot == %

            figure(rf)
            box on

            % Plot Functions
            fplot(k,[0 R],Color=Pcolour(inx,:),LineWidth=Plw)

            % Axes
            xlim([0 R]), set(gca,'XTick',[0 R/2 R])
            set(gca,'XTickLabel',{'0','$R/2$','$R$'},'TickLabelInterpreter','latex','fontsize',fsz)
            rnge = abs(kmax - kmin); mgin = 0.05;
            ylim([kmin-mgin*rnge kmax+mgin*rnge]), set(gca,'YTick',[kmin kavg kmax])
            set(gca,'YTickLabel',{'$k_{\rm{min}}$','$k_{\rm{avg}}$','$k_{\rm{max}}$'},...
                'TickLabelInterpreter','latex','fontsize',fsz)
            xlabel('$r$',Interpreter='latex',FontSize=lfsz)
            ylabel('$k(r)$',Interpreter='latex',FontSize=lfsz)
        end

        % Save Figure
        if savefigPF

            % Diffusivity
            figure(df)
            filename = ['diffusivity',num2str(d)];
            print(gcf,fullfile(pathPF,filename),'-depsc2')

            % Reaction-Rate
            figure(rf)
            filename = ['reaction-rate',num2str(d)];
            print(gcf,fullfile(pathPF,filename),'-depsc2')
        end
    end
end

% General Profiles
if enableG

    % with binding
    [~,th,F] = numerical(Rh,Dh,kh,Ph,c0hN,Th,Nr,Nt,itolN,d,coating);
    t = th*tn;

    % plot (with binding)
    figcount = figcount+1;
    wb = figure(figcount);
    clf(wb)
    hold on
    box on
    plot(t,F,LineWidth=Glw,Color=Gcolour(1,:))
    plot(t,F(end)*ones(size(t)),'--',LineWidth=Glw,Color=Gcolour(2,:))
    ylim([0 1])
    % axes
    xlim([0 T]), set(gca,'XTick',[0 T])
    set(gca,'XTickLabel',{'0','$T$'},'TickLabelInterpreter','latex','fontsize',labfsz)
    ylim([0 1]), set(gca,'YTick',[0 F(end) 1])
    set(gca,'YTickLabel',{'0','$F_{\infty}$','1'},'TickLabelInterpreter','latex','fontsize',labfsz)
    xl = xlabel('Time [$t$]',Interpreter='latex',FontSize=labfsz);
    xl.Position(2) = -0.05;
    ylabel('Fraction of Drug Released [$F(t)$]',Interpreter='latex',FontSize=labfsz)
    text(T/2,0.4,'with binding',Interpreter='latex',FontSize=capfsz,FontWeight='bold')

    % without binding
    [~,th,F] = numerical(Rh,Dh,@(r) 0,Ph,c0hN,Th,Nr,Nt,itolN,d,coating);
    t = th*tn;

    % plot (without binding)
    wob = figure(2);
    clf(wob)
    hold on
    box on
    plot(t,F,LineWidth=Glw,Color=Gcolour(1,:))
    ylim([0 1])
    % axes
    xlim([0 T]), set(gca,'XTick',[0 T])
    set(gca,'XTickLabel',{'0','$T$'},'TickLabelInterpreter','latex','fontsize',labfsz)
    ylim([0 1]), set(gca,'YTick',[0 1])
    set(gca,'YTickLabel',{'0','1'},'TickLabelInterpreter','latex','fontsize',labfsz)
    xl = xlabel('Time [$t$]',Interpreter='latex',FontSize=labfsz);
    xl.Position(2) = -0.05;
    ylabel('Fraction of Drug Released [$F(t)$]',Interpreter='latex',FontSize=labfsz)
    text(T/2,0.4,'without binding',Interpreter='latex',FontSize=capfsz,FontWeight='bold')

    % Save Figure
    if savefigGP

        % with binding
        figure(wb)
        filename = 'General_withbinding';
        print(gcf,fullfile(pathGP,filename),'-depsc2')

        % without binding
        figure(wob)
        filename = 'General_withoutbinding';
        print(gcf,fullfile(pathGP,filename),'-depsc2')
    end

end