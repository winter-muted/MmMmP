function [] = finElemMeth()
N = 8;      % number of elements

L = 5;  % length of beam
a = 0;  % start of beam
b = L;  % end of beam

n = N + 1;      % number of nodes
h = (b - a)/N;  % element size
w = 0.025;      % applied force/length in lb/ft
T = 1;    % Tension on string
h

f = -w*h/2;     % f term

o1 = ones(n,1);     % blank vector
M = spdiags([-o1,2*o1,-o1],-1:1,n,n);
M(1,2) = 0; M(2,1) = 0; M(n,n-1) = 0; M(n-1,n) = 0;
M(1,1) = 1; M(n,n) = 1;
K = T/h*M;  % coefficient matrix


fQ = o1*(2*f);
fQ(1) = 0;
fQ(n) = 0;

K
fQ

u = K\fQ;   % solution

min(u)

% x(:,1) = a:h:b;
%
% X = 0:.01:L;
% Y = -w/(2*T)*X.*(L-X);  % analytical solution
% y = -w/(2*T)*x.*(L-x);  % analytical solution for residual
% subplot(2,1,1)
% plot(x,u,'-+',X,Y)
% xlabel('x (ft)')
% ylabel('y (ft)')
% legend('FEM','analytical')
% subplot(2,1,2)
% plot(x,abs(u-y))
% xlabel('x (ft)')
% ylabel('residual (ft)')


end
