function pdePHK
%Adapted from PDEAuxRevKan, which solve a 3-component gradient mechanism 
%from 6/16/20 to 7/14/20. Changing variable names, 11/5/20.
% Altered Sep. 2021 to use Ka, rather than Kd, in Hill term.
% P is an upstream precursor of H. Could be a TF in general terms.
% H is an HD ZIP III. Since these are known in conifers generally. Not REV
% specifically. 
% K is Kanadi, or the basal inhibitor of the central HD ZIP III. 
%
%   The equation is to hold on an interval 0 <= x <= 3 for times t >= 0.
%   Boundary
%   conditions are specified so that it will be the solution of the initial-
%   boundary value problem.  The initial values are 0.
%   
%
%   The problem is coded in subfunctions SDDPDE, SDDIC, and SDDBC.
%
%   See also PDEPE, FUNCTION_HANDLE.

%   Lawrence F. Shampine and Jacek Kierzenka
%   Copyright 1984-2014 The MathWorks, Inc.

% m = 0 is for 'truly' 1D. Consider m = 1 for representation of
% 1D slice from disc. m = 2 would be for spherically symmetrical
% release from a point.

m = 0;
gridx = 300;
endt = 200;
gridt = 3;

x = linspace(0,3,gridx);
t = linspace(0,endt,gridt);
%disp(size(x));

C.DP = 5e-1;
C.C0P = 1e2;
C.decP = 1;
% old parameter, not used as of 9/20/21, C.prodP = 1;

% C.se is source extent on x = 0,3
C.se = 0.1;

C.DH = 1e-4;
%C.C0H = 0;
C.decH = 1e0;
C.prodH = 6e2;
C.Hill = 3;
C.KaHill = 3.2e1;
C.inhH = 2.4e0;


C.DK = 1e-4;
C.C0K = 100;
C.prodK = 1e2;
C.decK = 1;
C.inhK = 1e0;

C.ICH = (C.prodH/C.decH)*(C.C0P^C.Hill/(C.KaHill^C.Hill+C.C0P^C.Hill))*(1/(1+C.inhH*C.C0K));


eqn = @(x,t,u,dudx) PHKpde(x,t,u,dudx,C);
ic = @(x) PHKic(x,C);
bc = @(xl,ul,xr,ur,t) PHKbc(xl,ul,xr,ur,t,C);

sol = pdepe(m,eqn,ic,bc,x,t);
% Extract the first solution component as u.  This is not necessary
% for a single equation, but makes a point about the form of the output.
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);


% A surface plot is often a good way to study a solution.
%figure;
%surf(x,t,u1);
%title(['Num. soln:',num2str(gridx),' mesh points; ',num2str(gridt),' time steps']);
%xlabel('Distance x');
%ylabel('Time t');
%zlabel('PREC');

%figure;
%surf(x,t,u2);
%title(['Num. soln:',num2str(gridx),' mesh points; ',num2str(gridt),' time steps']);
%xlabel('Distance x');
%ylabel('Time t');
%zlabel('HDZ');

%figure;
%surf(x,t,u3);
%title(['Num. soln:',num2str(gridx),' mesh points; ',num2str(gridt),' time steps']);
%xlabel('Distance x');
%ylabel('Time t');
%zlabel('KAN');

% A solution profile can also be illuminating.
figure;
p=plot(x,u2(end,:),'-',x,u3(end,:),'-',x,u1(end,:),'-');
p(1).LineWidth = 4;
p(2).LineWidth = 4;
p(3).LineWidth = 1;
ylim([0; 100]);
title(['Solutions at t = ', num2str(endt)]);
legend('HDZ','KAN','PREC');
xlabel('Distance x');
ylabel(['Conc.']);



% --------------------------------------------------------------------------

function [c,f,s] = PHKpde(x,t,u,DuDx,C)

c = [1; 1; 1];
f = [C.DP; C.DH; C.DK] .* DuDx;
if x <= C.se
    s = [C.C0P - C.decP*u(1); (C.prodH*(u(1).^C.Hill)/(C.KaHill^C.Hill+u(1).^C.Hill)*(1/(1+C.inhH*u(3)))) - C.decH*u(2); C.prodK/(1+C.inhK*u(2)) - C.decK*u(3)];
else
    s = [- C.decP*u(1); (C.prodH*(u(1).^C.Hill)/(C.KaHill^C.Hill+u(1).^C.Hill)*(1/(1+C.inhH*u(3)))) - C.decH*u(2); C.prodK/(1+C.inhK*u(2)) - C.decK*u(3)];
end

   

% --------------------------------------------------------------------------

function u0 = PHKic(x,C)
u0 = [C.C0P; C.ICH; C.C0K];

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = PHKbc(xl,ul,xr,ur,t,C)

% LHS Neumann, exp slope
% pl = C.C0*sqrt(C.k/C.D);

% LHS DRC: C0
%pl = ul - C.C0;

%LHS DRC Aux; Neumann, no-flux Rev, Kan(i.e. setting LH bdry at plant tip, 1D slice of 
% radially symmetric pattern):
pl = [ul(1) - C.C0P; 0; 0];
ql = [0; 1; 1];

% RHS Neumann, no-flux for all:
%pr = [C.C0rev*sqrt(C.decrev/C.Drev)*exp(-sqrt(C.decrev/C.Drev)); 0];
pr = [0; 0; 0];

% RHS Neumann, exp slope (i.e. a soft RH bdry, the bdry is a position on a
% 'stalk' extending to the right
%pr = C.C0*sqrt(C.k/C.D)*exp(-sqrt(C.k/C.D));

%RHS Neumann no-flux for all:
qr = [1; 1; 1];


