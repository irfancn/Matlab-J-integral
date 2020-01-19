% J integral for the static 3PB PMMA

clear

load 'dic'
notch = [247, 455];             % notch of the specimen
frames = 10;                    % required instances
warning('off', 'all')

% Material props & initial params
E = 3e9;
nu = 0.35;
mu = E / (2 *(1 + nu) );
k = (3 - nu)/(1 + nu);  % plane stress
% k = 3 - 4*nu;           % plane strain
C = E/(1-nu*nu) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% spacing between the displacement data points.
space = data_dic_save.dispinfo.spacing + 1;
rat = space * data_dic_save.dispinfo.pixtounits;

% domain of interest
n = 40;
notch = round( notch/space );
rect = [notch(1)+1-1*n notch(1)+1*n; notch(2)+1-n notch(2)+n];

[x, y] = meshgrid( rat*linspace( -n, n, 2*n ), ...
    rat*linspace( -n, n, 2*n ) );

% Points under consideration for J integral.
rout = 0.8*n * rat;         % outer radius of domain
strip = 10 * rat;           % strip width
ang = 90 * pi/180;          % wedge angle 
tip = [rout, rat*n];
instr = 1.2e-3;
xm = x - min(x(:)) - tip(1);
ym = y - min(y(:)) - tip(2);
[tm, rm] = cart2pol(xm, ym);
ele = find( rm > rout-strip & rm < rout & tm > -pi+ang/2 & tm < pi-ang/2);
P = zeros(2*n);
P(ele) = 1;

% Defining Q, for the strip
q1 = -1/strip * cos(tm);    % smooth continuous functions
q2 = -1/strip * sin(tm);

Kj = zeros( frames, 2 );
for i = 1 : frames
% extracting the strain terms
e11 = data_dic_save.strains(i).plot_exx_cur_formatted;
e22 = data_dic_save.strains(i).plot_eyy_cur_formatted;
e12 = data_dic_save.strains(i).plot_exy_cur_formatted;

e11 = e11( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
e22 = e22( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
e12 = e12( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );

% evaluating the displacment derivatives
u1 = data_dic_save.displacements(i).plot_v_cur_formatted;
u2 = data_dic_save.displacements(i).plot_v_cur_formatted;
u1m1 = u1 + flipud( u1 );
u1m2 = u1 - flipud( u1 );
u2m1 = u2 - flipud( u2 );
u2m2 = u2 + flipud( u2 );

u2m1 = u2m1( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
u2m1 = u2m1 - mean(u2m1(:));
u2m1( u2m1 > 5*std(u2m1(:)) | u2m1 < -5*std(u2m1(:)) ) = 0;
u21m1 = 0.5 * ( u2m1(2:2*n-1, 1:2*n-2) + u2m1(2:2*n-1, 3:2*n) - ...
    2*u2m1(2:2*n-1, 2:2*n-1) )./diff( x(2:end-1, 1:2*n-1), 1,2);
u21tm1 = zeros(2*n); u21tm1(2:end-1, 2:end-1) = u21m1;
u21m1 = u21tm1 - mean(u21tm1(:));
u21m1 = smoothn( u21m1, 'robust' );

u2m2 = u2m2( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
u2m2 = u2m2 - mean(u2m2(:));
u2m2( u2m2 > 5*std(u2m2(:)) | u2m2 < -5*std(u2m2(:)) ) = 0;
u21m2 = 0.5 * ( u2m2(2:2*n-1, 1:2*n-2) + u2m2(2:2*n-1, 3:2*n) - ...
    2*u2m2(2:2*n-1, 2:2*n-1) )./diff( x(2:end-1, 1:2*n-1), 1,2);
u21tm2 = zeros(2*n); u21tm2(2:end-1, 2:end-1) = u21m2;
u21m2 = u21tm2 - mean(u21tm2(:));
u21m2 = smoothn( u21m2, 'robust' );

u1m1 = u1m1( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
u1m1 = u1m1 - mean(u1m1(:));
u1m1( u1m1 > 5*std(u1m1(:)) | u1m1 < -5*std(u1m1(:)) ) = 0;
u12m1 = 0.5 * ( u1m1(2:2*n-1, 1:2*n-2) + u1m1(2:2*n-1, 3:2*n) - ...
    2*u1m1(2:2*n-1, 2:2*n-1) )./diff( y(2:end, 2:2*n-1), 1,1);
u12tm1 = zeros(2*n); u12tm1(2:end-1, 2:end-1) = u12m1;
u12m1 = u12tm1 - mean(u12tm1(:));
u12m1 = smoothn( u12m1, 'robust' );

u1m2 = u1m2( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
u1m2 = u1m2 - mean(u1m2(:));
u1m2( u1m2 > 5*std(u1m2(:)) | u1m2 < -5*std(u1m2(:)) ) = 0;
u12m2 = 0.5 * ( u1m2(2:2*n-1, 1:2*n-2) + u1m2(2:2*n-1, 3:2*n) - ...
    2*u1m2(2:2*n-1, 2:2*n-1) )./diff( y(2:end, 2:2*n-1), 1,1);
u12tm2 = zeros(2*n); u12tm2(2:end-1, 2:end-1) = u12m2;
u12m2 = u12tm2 - mean(u12tm2(:));
u12m2 = smoothn( u12m2, 'robust' );

u2 = data_dic_save.displacements(i).plot_v_cur_formatted;
u2 = u2( rect(2,1):rect(2,2), rect(1,1):rect(1,2) );
u21 = 0.5*(u2(:, 3:end) + u2(:, 1:end-2) - 2*u2(:, 2:end-1))/ (x(2,2) - x(1,1));
u21( u21 > 5*std(u21(:)) | u21 < -5*std(u21(:)) ) = 0;
u21t = zeros(2*n); u21t(:, 2:end-1) = u21;  
u21 = smoothn( u21t, 'robust' );

e1m1 = 0.5*( e11 + flipud( e11 ) );
e2m1 = 0.5*( e22 - flipud( e22 ) );
e3m1 = e12;%u12m1 + u21m1;

e1m2 = 0.5*( e11 - flipud( e11 ) );
e2m2 = 0.5*( e22 + flipud( e22 ) );
e3m2 = e12;%u12m2 + u21m2;

% mode I and mode II contributio to J1 and J2
J1 = 0; J2 = 0;
for j = 1 : 2*n
for k = 1 : 2*n
if( P(j,k) )
    S1 = C * [e1m1(j, k), e2m1(j, k), e3m1(j, k)]';
    S2 = C * [e1m2(j, k), e2m2(j, k), e3m2(j, k)]';
    W1 = 0.5 * sum( S1.*[e1m1(j, k), e2m1(j, k), e3m1(j, k)]' );
    W2 = 0.5 * sum( S2.*[e1m2(j, k), e2m2(j, k), e3m2(j, k)]' );
    Jin1 = rat^2 * ( q1(j,k)*( S1(1)*e1m1(j,k) + S1(3)*u21m1(j,k) - 2*W1 )...
        + q2(j,k)*( S1(3)*e1m1(j,k) + S1(2)*u21m1(j,k) ) );
    J1 = J1 + Jin1;
    Jin2 = rat^2 * ( q1(j,k)*( S2(1)*e1m2(j,k) + S2(3)*u21m2(j,k) - 2*W2 )...
        + q2(j,k)*( S2(3)*e1m2(j,k) + S2(2)*u21m2(j,k) ) );
    J2 = J2 + Jin2;
end
end
end

% evaluating the SIF from the J-integral
Kj(i, :) = real( sqrt( [J1 J2] * E ) );
end

% plotting the K1 and K2 evaluated from the J integral
figure, plot( Kj )
xlabel('Frames'), ylabel('SIF (MPa s^{1/2})')

% plotting the contour chosen for the J integral
figure, contourf( e22, 100 ), hold on, contour( P ), axis equal
title('Domain of interst for the J-integral evaluation')