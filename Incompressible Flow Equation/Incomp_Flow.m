clear all
close all
clc

% Maximum number of iterations
max_iter = 1000;
tolerance = 10^-15;
lambda = 0.3; %under relaxation factor
iter = 0;
angle = pi/3;
ro = 0.8;
mu = 5*10^-5;
cp = 1.03;
k = 0.036;

T_inf = 300;
h_inf = 15;

% Define the vertices of the quadrilateral
vertices = [0, 0; 1, 0; 1+cos(angle), sin(angle); cos(angle), sin(angle)];

% Discretization parameters
ny = 6; % Number of points along y
nx = 6; % Number of points along x

% Initialize matrices for grid points
X = zeros(ny, nx);
Y = zeros(ny, nx);

% Interpolation parameters
xi = linspace(0, 1, ny);
eta = linspace(0, 1, nx);

% Iterate over the grid points
for i = 1:ny
    for j = 1:nx
        % Interpolation along xi and eta axes
        xi_val = xi(i);
        eta_val = eta(j);

        % Transfinite interpolation formula for each coordinate
        X(i, j) = (1 - xi_val) * ((1 - eta_val) * vertices(1, 1) + eta_val * vertices(2, 1)) + ...
            xi_val * ((1 - eta_val) * vertices(4, 1) + eta_val * vertices(3, 1));
        Y(i, j) = (1 - eta_val) * ((1 - xi_val) * vertices(1, 2) + xi_val * vertices(4, 2)) + ...
            eta_val * ((1 - xi_val) * vertices(2, 2) + xi_val * vertices(3, 2));
    end
end

% Plot the quadrilateral grid
figure;
plot(X, Y, 'b-', X', Y', 'b-');
xlabel('X');
ylabel('Y');
title('Transfinite Interpolation for Quadrilateral Grid');
axis equal;
grid on;

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Centroids
Xc = zeros(ny-1,nx-1);
Yc = zeros(ny-1,nx-1);

for i=1:ny-1
    for j=1:nx-1
        Xc(i,j) = mean ([X(i,j), X(i,j+1),X(i+1,j),X(i+1,j+1)]);
        Yc(i,j) = mean ([Y(i,j), Y(i,j+1),Y(i+1,j),Y(i+1,j+1)]);
    end
end


%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Face Centroids
Xfh = zeros(ny,nx-1);
Yfh = zeros(ny,nx-1);
Xfv = zeros(ny-1,nx);
Yfv = zeros(ny-1,nx);

for i=1:ny-1
    for j=1:nx
        Xfv(i,j) = mean ([X(i,j),X(i+1,j)]);
        Yfv(i,j) = mean ([Y(i,j),Y(i+1,j)]);

    end
end

for i=1:ny
    for j=1:nx-1
        Xfh(i,j) = mean ([X(i,j),X(i,j+1)]);
        Yfh(i,j) = mean ([Y(i,j),Y(i,j+1)]);
    end
end


%scatter(Xfh,Yfh,'filled','go');
hold on
scatter(Xfv,Yfv,'filled','bo');

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Centroids and centroids at Boundary Faces

Xfc = zeros(ny+1,nx+1);
Yfc = zeros(ny+1,nx+1);

Xfc(1,1) = vertices(1,1);
Yfc(1,1) = vertices(1,2);
Xfc(1,nx+1) = vertices(2,1);
Yfc(1,nx+1) = vertices(2,2);
Xfc(ny+1,nx+1) = vertices(3,1);
Yfc(ny+1,nx+1) = vertices(3,2);
Xfc(ny+1,1) = vertices(4,1);
Yfc(ny+1,1) = vertices(4,2);


for i=1:ny-1
    for j=1:nx-1

        %AB
        Xfc(1,j+1) = Xfh(1,j);
        Yfc(1,j+1) = Yfh(1,j);

        %CD
        Xfc(ny+1,j+1) = Xfh(ny,j);
        Yfc(ny+1,j+1) = Yfh(ny,j);

        %AD
        Xfc(i+1,1) = Xfv(i,1);
        Yfc(i+1,1) = Yfv(i,1);

        %BC
        Xfc(i+1,nx+1) = Xfv(i,nx);
        Yfc(i+1,nx+1) = Yfv(i,nx);

        %Internal Centroids
        Xfc(i+1,j+1) = Xc(i,j);
        Yfc(i+1,j+1) = Yc(i,j);
    end
end

%hold on
%scatter(Xfc,Yfc,'filled','ro');



%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Centroids distance

de = zeros(ny+1,nx);
dn = zeros(ny,nx+1);


for i=1:ny+1
    for j=nx
        de(i,j) = sqrt((Xfc(i,j+1) - Xfc(i,j))^2 + (Yfc(i,j+1) - Yfc(i,j))^2);
    end
end

for i=1:ny
    for j=1:nx+1
        dn(i,j) = sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2);
    end
end


%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Areas (Volumes)

A = zeros(ny-1,nx-1);

p1 = zeros(ny-1,nx-1);
p2 = zeros(ny-1,nx-1);

for i=1:ny-1
    for j=1:nx-1

        p1(i,j) = abs((X(i+1,j)-X(i,j))*(Y(i,j+1)-Y(i,j))-(X(i,j+1)-X(i,j))*(Y(i+1,j)-Y(i,j)));
        p2(i,j) = abs((X(i+1,j+1)-X(i,j+1))*(Y(i+1,j+1)-Y(i+1,j))-(X(i+1,j+1)-X(i+1,j))*(Y(i+1,j+1)-Y(i,j+1)));
        A(i,j) = 0.5*(p1(i,j)+p2(i,j));
    end
end



%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Faces Areas

Afe = zeros(ny-1,nx);
Afn = zeros(ny,nx-1);


for i=1:ny-1
    for j=1:nx

        Afe(i,j) = sqrt((X(i+1,j) - X(i,j))^2 + (Y(i+1,j) - Y(i,j))^2);

    end
end

for i=1:ny
    for j=1:nx-1

        Afn(i,j) = sqrt((X(i,j+1) - X(i,j))^2 + (Y(i,j+1) - Y(i,j))^2);

    end
end

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Normal Vectors

Sxe = zeros(ny-1,nx); % East Normal Vector coordinates
sxe = zeros(ny-1,nx);
Sye = zeros(ny-1,nx);
sye = zeros(ny-1,nx);

Sxn = zeros(ny,nx-1); % North Normal Vector coordinates
sxn = zeros(ny,nx-1);
Syn = zeros(ny,nx-1);
syn = zeros(ny,nx-1);


for i=1:ny-1
    for j=1:nx
        Sxe(i,j) = Y(i+1,j)-Y(i,j);
        Sye(i,j) = -X(i+1,j)+X(i,j);
        sxe(i,j) = Sxe(i,j)/sqrt((Y(i+1,j)-Y(i,j))^2 + (-X(i+1,j)+X(i,j))^2);
        sye(i,j) = Sye(i,j)/sqrt((Y(i+1,j)-Y(i,j))^2 + (-X(i+1,j)+X(i,j))^2);

    end
end


for i=1:ny
    for j=1:nx-1
        Sxn(i,j) = -Y(i,j+1)+Y(i,j);
        Syn(i,j) = X(i,j+1)-X(i,j);
        sxn(i,j) = Sxn(i,j)/sqrt((-Y(i,j+1)+Y(i,j))^2 + (X(i,j+1)-X(i,j))^2);
        syn(i,j) = Syn(i,j)/sqrt((-Y(i,j+1)+Y(i,j))^2 + (X(i,j+1)-X(i,j))^2);

    end
end

quiver(Xfv,Yfv,Sxe,Sye,'k');
%quiver(Xfh,Yfh,Sxn,Syn,'filled','k');

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% East E - T Vectors

Exe = zeros(ny-1,nx);
Eye = zeros(ny-1,nx);
Txe = zeros(ny-1,nx);
Tye = zeros(ny-1,nx);
exe = zeros(ny-1,nx);
eye = zeros(ny-1,nx);
mod_ee = zeros(ny-1,nx);
mod_Ee = zeros(ny-1,nx);

for i=1:ny-1
    for j=1:nx


        mod_ee(i,j) = sqrt((Xfc(i+1,j+1)-Xfc(i+1,j))^2+(Yfc(i+1,j+1)-Yfc(i+1,j))^2);
        exe(i,j) = (Xfc(i+1,j+1)-Xfc(i+1,j))/mod_ee(i,j);
        eye(i,j) = (Yfc(i+1,j+1)-Yfc(i+1,j))/mod_ee(i,j);
        Exe(i,j) = dot([Sxe(i,j) Sye(i,j)],[Sxe(i,j) Sye(i,j)])*exe(i,j)/dot([exe(i,j) eye(i,j)],[Sxe(i,j) Sye(i,j)]);
        Eye(i,j) = dot([Sxe(i,j) Sye(i,j)],[Sxe(i,j) Sye(i,j)])*eye(i,j)/dot([exe(i,j) eye(i,j)],[Sxe(i,j) Sye(i,j)]);
        mod_Ee(i,j) = sqrt(Exe(i,j)^2 + Eye(i,j)^2);
        Txe(i,j) = Sxe(i,j) - Exe(i,j);
        Tye(i,j) = Sye(i,j) - Eye(i,j);


    end
end

hold on
quiver(Xfv,Yfv,Exe,Eye,'b');
quiver(Xfv,Yfv,Txe,Tye,'b');


%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% North E - T Vectors

Exn = zeros(ny,nx-1);
Eyn = zeros(ny,nx-1);
Txn = zeros(ny,nx-1);
Tyn = zeros(ny,nx-1);
exn = zeros(ny,nx-1);
eyn = zeros(ny,nx-1);
mod_en = zeros(ny,nx-1);
mod_En = zeros(ny,nx-1);


for i=1:ny
    for j=1:nx-1

        mod_en(i,j) = sqrt((Xfc(i+1,j+1)-Xfc(i,j+1))^2+(Yfc(i+1,j+1)-Yfc(i,j+1))^2);
        exn(i,j) = (Xfc(i+1,j+1)-Xfc(i,j+1))/mod_en(i,j);
        eyn(i,j) = (Yfc(i+1,j+1)-Yfc(i,j+1))/mod_en(i,j);
        Exn(i,j) = dot([Sxn(i,j) Syn(i,j)],[Sxn(i,j) Syn(i,j)])*exn(i,j)/dot([exn(i,j) eyn(i,j)],[Sxn(i,j) Syn(i,j)]);
        Eyn(i,j) = dot([Sxn(i,j) Syn(i,j)],[Sxn(i,j) Syn(i,j)])*eyn(i,j)/dot([exn(i,j) eyn(i,j)],[Sxn(i,j) Syn(i,j)]);
        mod_En(i,j) = sqrt(Exn(i,j)^2 + Eyn(i,j)^2);
        Txn(i,j) = Sxn(i,j) - Exn(i,j);
        Tyn(i,j) = Syn(i,j) - Eyn(i,j);

    end
end

%hold on
%quiver(Xfv,Yfv,Exn,Eyn,'filled','g');
%quiver(Xfv,Yfv,Txn,Tyn,'filled','g');

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% gDiff

g_Diffe = zeros(ny-1,nx);
g_Diffn = zeros(ny,nx-1);

for i=1:ny-1
    for j=1:nx

        g_Diffe(i,j) = sqrt(Exe(i,j)^2 + Eye(i,j)^2)/sqrt((Xfc(i,j+1) - Xfc(i,j))^2 + (Yfc(i,j+1) - Yfc(i,j))^2);

    end
end

for i=1:ny
    for j=1:nx-1

        g_Diffn(i,j) = sqrt(Exn(i,j)^2 + Eyn(i,j)^2)/sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2);

    end
end

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Interpolation factors

ge = zeros(ny-1,nx); % for east faces
gn = zeros(ny,nx-1); % for north faces


for i=1:ny-1

    ge(i,1) = 1;
    ge(i,nx) = 1;

    for j=2:nx-1
        ge(i,j) = (A(i,j-1))/(A(i,j-1) + A(i,j));
    end
end


for j=1:nx-1

    gn(1,j) = 1;
    gn(ny,j) = 1;

    for i=2:ny-1
        gn(i,j) = (A(i-1,j))/(A(i-1,j) + A(i,j));
    end
end

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


u = zeros(ny+1,nx+1);      % velocity in x-direction
v = zeros(ny+1,nx+1);      % velocity in y-direction

uold = zeros(ny+1,nx+1);

vold = zeros(ny+1,nx+1);

un = zeros(ny,nx-1);       % velocity in x-direction
vn = zeros(ny,nx-1);       % velocity in y-direction

ue = zeros(ny-1,nx);       % velocity in x-direction
ve = zeros(ny-1,nx);       % velocity in y-direction

mdotn = zeros(ny,nx-1);    % mass flow rate for north faces

mdote = zeros(ny-1,nx);    % mass flow rate for east faces

grad_uf = ones(ny-1,nx);

grad_ue = zeros(ny-1,nx);

grad_un = zeros(ny,nx-1);

grad_u = zeros(ny+1,nx+1);

grad_v = zeros(ny+1,nx+1);

grad_vf = ones(ny-1,nx);

grad_ve = zeros(ny-1,nx);

grad_vn = zeros(ny,nx-1);

grad_p = ones(ny+1,nx+1);

grad_pe = zeros(ny-1,nx);

grad_pn = zeros(ny,nx-1);

Sc = zeros(ny+1,nx+1);

Dcu = zeros(ny+1,nx+1);

Dcv = zeros(ny+1,nx+1);

Dfu = zeros(ny-1,nx);

Dfv = zeros(ny,nx-1);

p = zeros(ny+1,nx+1);

p_prime = zeros(ny+1,nx+1);

p_pold = zeros(ny+1,nx+1);

%==================================================================================================================================================================================================================================================
% u - velocity calculation

%while error>tolerance && iter<max_iter
for iter=1:max_iter

    aE = zeros(ny+1,nx+1);
    aW = zeros(ny+1,nx+1);
    aN = zeros(ny+1,nx+1);
    aS = zeros(ny+1,nx+1);
    aC = zeros(ny+1,nx+1);
    bC = zeros(ny+1,nx+1);


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % u-velocity at the faces

    % East


    for i=1:ny-1

        ue(i,1) = u(i+1,1);
        ue(i,nx) = u(i+1,nx+1);

        for j=2:nx-1
            ue(i,j) = (1 - ge(i,j))*u(i+1,j) + ge(i,j)*u(i+1,j+1);
        end
    end

    % North


    for j=1:nx-1

        un(1,j) = u(1,j+1);
        un(ny,j) = u(ny+1,j+1);

        for i=2:ny-1
            un(i,j) = (1 - gn(i,j))*u(i,j+1) + gn(i,j)*u(i+1,j+1);
        end
    end


    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of u-velocity at the centroids


    for i=2:ny
        for j=2:nx

            grad_u(i,j) = (-ue(i-1,j-1)*Sxe(i-1,j-1) + ue(i-1,j)*Sxe(i-1,j) - un(i-1,j-1)*Sxn(i-1,j-1) + un(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);

        end
    end

    % West Boundary
    grad_u(:,1) = grad_u(:,2);

    % South Boundary
    grad_u(1,:) = grad_u(2,:);

    % East Boundary
    grad_u(:,nx+1) = grad_u(:,nx);

    % North Boundary
    grad_u(ny+1,:) = grad_u(ny,:);

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of u-velocity at the faces

    for i=1:ny-1

        grad_uf(i,1) = grad_u(i+1,1);
        grad_uf(i,nx) = grad_u(i+1,nx+1);


        for j=2:nx-1

            grad_uf(i,j) = (1 - ge(i,j))*grad_u(i+1,j) + ge(i,j)*grad_u(i+1,j+1);


        end
    end

    grad_ue(:,:) = grad_uf(:,:);

    %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Loop for inner elements:

    for i=3:ny-1
        for j=3:nx-1

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------========================================================================================================================================================================
    % Loop for WEST elements without vertices

    for j=2
        for i=3:ny-1

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) +  aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aW(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);
        end
    end


    % Loop over SOUTH elements without vertices: u=0 (no slip condition)
    % wall

    for i=2
        for j=3:nx-1

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aS(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    % Loop over NORTH elements without vertices

    for i=ny
        for j=3:nx-1

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aN(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    % Loop over EAST elements without vertices: u=0 (no slip condition)

    for j=nx
        for i=3:ny-1

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j) + aE(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aE(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j)  + aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 5*aW(i,j) - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aW(i,j) = 0;
            aS(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny
            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aS(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aE(i,j) = 0;
            aN(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aS(i,j) = 0;
            aE(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - NORTH CORNER

    for j=2
        for i=ny

            Non_Ortho = mu*(grad_ue(i-1,j)*Txe(i-1,j)) + mu*(grad_ue(i-1,j-1)*Txe(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j)  + aS(i,j) + aW(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + mu*(grad_ue(i-1,j)*Sxe(i-1,j) + grad_ue(i-1,j-1)*Sxe(i-1,j-1));
            aW(i,j) = 0;
            aN(i,j) = 0;
            Dcu(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % SOlVER: Gauss Seidel


    for i=2:ny
        for j=2:nx

            uold(i,j) = u(i,j);

            u(i,j) = u(i,j) + lambda*((bC(i,j) - aE(i,j)*u(i,j+1) - aN(i,j)*u(i+1,j) - aW(i,j)*u(i,j-1) - aS(i,j)*u(i-1,j))/aC(i,j) - u(i,j));

            error = abs(u(i,j) - uold(i,j));


        end
    end




    %==================================================================================================================================================================================================================================================
    % v - velocity calculation

    aE = zeros(ny+1,nx+1);
    aW = zeros(ny+1,nx+1);
    aN = zeros(ny+1,nx+1);
    aS = zeros(ny+1,nx+1);
    aC = zeros(ny+1,nx+1);
    bC = zeros(ny+1,nx+1);


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % v-velocity at the faces

    % East


    for i=1:ny-1

        ve(i,1) = v(i+1,1);
        ve(i,nx) = v(i+1,nx+1);

        for j=2:nx-1
            ve(i,j) = (1 - ge(i,j))*v(i+1,j) + ge(i,j)*v(i+1,j+1);
        end
    end

    % North


    for j=1:nx-1

        vn(1,j) = v(1,j+1);
        vn(ny,j) = v(ny+1,j+1);

        for i=2:ny-1
            vn(i,j) = (1 - gn(i,j))*v(i,j+1) + gn(i,j)*v(i+1,j+1);
        end
    end


    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of v-velocity at the centroids


    for i=2:ny
        for j=2:nx

            grad_v(i,j) = (-ve(i-1,j-1)*Sxe(i-1,j-1) + ve(i-1,j)*Sxe(i-1,j) - vn(i-1,j-1)*Sxn(i-1,j-1) + vn(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);

        end
    end

    % West Boundary
    grad_v(:,1) = grad_v(:,2);

    % South Boundary
    grad_v(1,:) = grad_v(2,:);

    % East Boundary
    grad_v(:,nx+1) = grad_v(:,nx);

    % North Boundary
    grad_v(ny+1,:) = grad_v(ny,:);

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of v-velocity at the faces

    for i=1:ny-1

        grad_vf(i,1) = grad_v(i+1,1);
        grad_vf(i,nx) = grad_v(i+1,nx+1);


        for j=2:nx-1

            grad_vf(i,j) = (1 - gn(i,j))*grad_v(i,j+1) + gn(i,j)*grad_v(i+1,j+1);


        end
    end

grad_ve(:,:) = grad_vf(:,:);
    %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Loop for inner elements:

    for i=3:ny-1
        for j=3:nx-1

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);


        end
    end

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------========================================================================================================================================================================
    % Loop for WEST elements without vertices

    for j=2
        for i=3:ny-1

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) +  aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aW(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);
        end
    end


    % Loop over SOUTH elements without vertices: u=0 (no slip condition)
    % wall

    for i=2
        for j=3:nx-1

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));            aE(i,j) = - mu * g_Diffe(i-1,j);
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aS(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    % Loop over NORTH elements without vertices

    for i=ny
        for j=3:nx-1

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aN(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    % Loop over EAST elements without vertices: u=0 (no slip condition)

    for j=nx
        for i=3:ny-1

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j) + aE(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aE(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j)  + aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 0*aW(i,j) - grad_p(i,j)*A(i-1,j-1) + product;
            aW(i,j) = 0;
            aS(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny
            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aS(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aE(i,j) = 0;
            aN(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aS(i,j) = 0;
            aE(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - NORTH CORNER

    for j=2
        for i=ny

            Non_Ortho = mu*(grad_ve(i-1,j)*Tye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Tye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Tyn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Tyn(i-1,j-1));
            product = mu*(grad_ve(i-1,j)*Sye(i-1,j)) + mu*(grad_ve(i-1,j-1)*Sye(i-1,j-1)) + mu*(grad_vn(i,j-1)*Syn(i,j-1)) + mu*(grad_vn(i-1,j-1)*Syn(i-1,j-1));
            aE(i,j) = - mu * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - mu * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - mu * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - mu * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j)  + aS(i,j) + aW(i,j) + aN(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - grad_p(i,j)*A(i-1,j-1) + product;
            aW(i,j) = 0;
            aN(i,j) = 0;
            Dcv(i,j) = A(i-1,j-1)/aC(i,j);

        end
    end
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % SOlVER: Gauss Seidel


    for i=2:ny
        for j=2:nx

            vold(i,j) = v(i,j);

            v(i,j) = v(i,j) + lambda*((bC(i,j) - aE(i,j)*v(i,j+1) - aN(i,j)*v(i+1,j) - aW(i,j)*v(i,j-1) - aS(i,j)*v(i-1,j))/aC(i,j) - v(i,j));

            error = abs(v(i,j) - vold(i,j));


        end
    end



    %==================================================================================================================================================================================================================================================
    % Rhie-Chow Interpolation, mass flow rates update



    % Dfu


    for i=1:ny-1

        Dfu(i,1) = Dcu(i+1,1);
        Dfu(i,nx) = Dcu(i+1,nx+1);

        for j=2:nx-1
            Dfu(i,j) = (1 - ge(i,j))*Dcu(i+1,j) + ge(i,j)*Dcu(i+1,j+1);
        end
    end

    % Dfn


    for j=1:nx-1

        Dfv(1,j) = Dcv(1,j+1);
        Dfv(ny,j) = Dcv(ny+1,j+1);

        for i=2:ny-1
            Dfv(i,j) = (1 - gn(i,j))*Dcv(i,j+1) + gn(i,j)*Dcv(i+1,j+1);
        end
    end

    for i=1:ny-1
        for j=1:nx
            ue(i,j) = ((1 - ge(i,j))*u(i+1,j) + ge(i,j)*u(i+1,j+1)) - Dfu(i,j)*(grad_pe(i,j) - ((1 - ge(i,j))*grad_p(i+1,j) + ge(i,j)*grad_p(i+1,j+1)));
            mdote(i,j) = ue(i,j)*ro*sqrt(Sxe(i,j)^2+Sye(i,j)^2);
        end
    end

    for j=1:nx-1
        for i=1:ny
            un(i,j) = un(i,j) - Dfv(i,j)*(grad_pn(i,j) - ((1 - gn(i,j))*grad_p(i,j+1) + gn(i,j)*grad_p(i+1,j+1)));
            mdotn(i,j) = un(i,j)*ro*sqrt(Sxn(i,j)^2+Syn(i,j)^2);
        end
    end



    %==================================================================================================================================================================================================================================================
    % pressure correction

    aE = zeros(ny+1,nx+1);
    aW = zeros(ny+1,nx+1);
    aN = zeros(ny+1,nx+1);
    aS = zeros(ny+1,nx+1);
    aC = zeros(ny+1,nx+1);
    bC = zeros(ny+1,nx+1);

    %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Loop for inner elements:
    
    Dfe(:,:) = Dfv(:,:);

    for i=3:ny-1
        for j=3:nx-1

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------========================================================================================================================================================================
    % Loop for WEST elements without vertices

    for j=2
        for i=3:ny-1

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));

        end
    end


    % Loop over SOUTH elements without vertices: u=0 (no slip condition)
    % wall

    for i=2
        for j=3:nx-1

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aW(i,j) + aN(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end

    % Loop over NORTH elements without vertices

    for i=ny
        for j=3:nx-1

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end

    % Loop over EAST elements without vertices: u=0 (no slip condition)

    for j=nx
        for i=3:ny-1

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aN(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));

        end
    end


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aW(i,j) + aN(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - NORTH CORNER

    for j=2
        for i=ny

            aE(i,j) = - ro*cp * ((Xfc(i,j+1)-Xfc(i,j))*Dfe(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2) +  (Yfc(i,j+1) - Yfc(i,j))*Dfv(i-1,j)*sqrt(Sxe(i-1,j)^2+Sye(i-1,j)^2));
            aW(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i,j-1))*Dfe(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i,j-1))*Dfv(i-1,j-1)*sqrt(Sxe(i-1,j-1)^2+Sye(i-1,j-1)^2));
            aN(i,j) = - ro*cp * ((Xfc(i+1,j)-Xfc(i,j))*Dfe(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2) +  (Yfc(i+1,j) - Yfc(i,j))*Dfv(i,j-1)*sqrt(Sxn(i,j-1)^2+Syn(i,j-1)^2));
            aS(i,j) = - ro*cp * ((Xfc(i,j)-Xfc(i-1,j))*Dfe(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2) +  (Yfc(i,j) - Yfc(i-1,j))*Dfv(i-1,j-1)*sqrt(Sxn(i-1,j-1)^2+Syn(i-1,j-1)^2));
            aC(i,j) = - (aE(i,j) + aS(i,j));
            bC(i,j) = - (mdote(i-1,j) - mdote(i-1,j-1) + mdotn(i,j-1) - mdotn(i-1,j-1));


        end
    end
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % SOlVER: Gauss Seidel


    for i=2:ny
        for j=2:nx

            p_pold(i,j) = p_prime(i,j);

            p_prime(i,j) = p_prime(i,j) + lambda*((bC(i,j) - aE(i,j)*p_prime(i,j+1) - aN(i,j)*p_prime(i+1,j) - aW(i,j)*p_prime(i,j-1) - aS(i,j)*p_prime(i-1,j))/aC(i,j) - p_prime(i,j));

            error = abs(p_prime(i,j) - p_pold(i,j));


        end
    end

    for i=1:ny-1
        for j=1:nx
            ue(i,j) = ue(i,j) - Dfu(i,j)*grad_ppe(i,j);
            mdote(i,j) = mdote(i,j) + -ro*cp*Dfu(i,j)*sqrt(Sxe(i,j)^2+Sye(i,j)^2)*grad_ppe(i,j);
        end
    end

    for i=1:ny
        for j=1:nx-1
            un(i,j) = un(i,j) - Dfu(i,j)*grad_ppn(i,j);
            mdotn(i,j) = mdotn(i,j) + -ro*cp*Dfu(i,j)*sqrt(Sxn(i,j)^2+Syn(i,j)^2)*grad_ppn(i,j);
        end
    end

    for i=1:ny-1
        for j=1:nx
            ve(i,j) = ve(i,j) - Dfv(i,j)*grad_ppe(i,j);

        end
    end

    for i=1:ny
        for j=1:nx-1
            vn(i,j) = vn(i,j) - Dfv(i,j)*grad_ppn(i,j);

        end
    end

    for i=2:ny
        for j=2:nx
            p(i,j) = p(i,j) + p_prime(i,j);
        end
    end



    %iter = iter+1;
end

%==================================================================================================================================================================================================================================================
% Energy Equation

aE = zeros(ny+1,nx+1);
aW = zeros(ny+1,nx+1);
aN = zeros(ny+1,nx+1);
aS = zeros(ny+1,nx+1);
aC = zeros(ny+1,nx+1);
bC = zeros(ny+1,nx+1);
grad_Txe = zeros(ny-1,nx);
grad_Tye = zeros(ny-1,nx);
grad_Txn = zeros(ny,nx-1);
grad_Tyn = zeros(ny,nx-1);
Tfe = zeros (ny-1,nx);
Tfn = zeros(ny,nx-1);
Sc = zeros (size(Xfc));
grad_Tx = zeros(ny+1,nx+1);
grad_Ty = zeros(ny+1,nx+1);
fluxCb = zeros(ny-1,nx);
fluxVb = zeros(ny-1,nx);
T = zeros(ny+1,nx+1);
Told = zeros(ny+1,nx+1);


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% HERE THE LOOP SHOULD START

while error>tolerance && iter<max_iter


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Temperature at the faces

    % East


    for i=1:ny-1

        Tfe(i,1) = T(i+1,1);
        Tfe(i,nx) = T(i+1,nx+1);

        for j=2:nx-1
            Tfe(i,j) = (1 - ge(i,j))*T(i+1,j) + ge(i,j)*T(i+1,j+1);
        end
    end

    % North


    for j=1:nx-1

        Tfn(1,j) = T(1,j+1);
        Tfn(ny,j) = T(ny+1,j+1);

        for i=2:ny-1
            Tfn(i,j) = (1 - gn(i,j))*T(i,j+1) + gn(i,j)*T(i+1,j+1);
        end
    end


    %-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Source term at centroids

    for i=1:ny+1
        for j=1:nx+1

            Sc(i,j) = 0;

        end
    end

    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of temperature at the centroids


    for i=2:ny
        for j=2:nx
            grad_Tx(i,j) = (-Tfe(i-1,j-1)*Sxe(i-1,j-1) + Tfe(i-1,j)*Sxe(i-1,j) - Tfn(i-1,j-1)*Sxn(i-1,j-1) + Tfn(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);
            grad_Ty(i,j) = (-Tfe(i-1,j-1)*Sye(i-1,j-1) + Tfe(i-1,j)*Sye(i-1,j) - Tfn(i-1,j-1)*Syn(i-1,j-1) + Tfn(i,j-1)*Syn(i,j-1))/A(i-1,j-1);
            %grad_Tx(i,j) = 1;
            %grad_Ty(i,j) = 1;

        end
    end

    % West Boundary
    grad_Tx(:,1) = grad_Tx(:,2);
    grad_Ty(:,1) = grad_Ty(:,2);

    % South Boundary
    grad_Tx(1,:) = grad_Tx(2,:);
    grad_Ty(1,:) = grad_Ty(2,:);

    % East Boundary
    grad_Tx(:,nx+1) = grad_Tx(:,nx);
    grad_Ty(:,nx+1) = grad_Ty(:,nx);

    % North Boundary
    grad_Tx(ny+1,:) = 0;
    grad_Ty(ny+1,:) = 0;


    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Gradient of temperature at the faces

    for i=1:ny-1

        grad_Txe(i,1) = grad_Tx(i+1,1);
        grad_Tye(i,1) = grad_Ty(i+1,1);
        grad_Txe(i,nx) = grad_Tx(i+1,nx+1);
        grad_Tye(i,nx) = grad_Ty(i+1,nx+1);


        for j=2:nx-1

            grad_Txe(i,j) = (1 - ge(i,j))*grad_Tx(i+1,j) + ge(i,j)*grad_Tx(i+1,j+1);
            grad_Tye(i,j) = (1 - ge(i,j))*grad_Ty(i+1,j) + ge(i,j)*grad_Ty(i+1,j+1);

        end
    end

    for j=1:nx-1

        grad_Txn(1,j) = grad_Tx(1,j+1);
        grad_Tyn(1,j) = grad_Ty(1,j+1);
        grad_Txn(ny,j) = grad_Tx(ny+1,j+1);
        grad_Tyn(ny,j) = grad_Ty(ny+1,j+1);

        for i=2:ny-1

            grad_Txn(i,j) = (1-gn(i,j))*grad_Tx(i,j+1) + gn(i,j)*grad_Tx(i+1,j+1);
            grad_Tyn(i,j) = (1-gn(i,j))*grad_Ty(i,j+1) + gn(i,j)*grad_Ty(i+1,j+1);

        end
    end

    %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Loop for inner elements:

    for i=3:ny-1
        for j=3:nx-1

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho;


        end
    end

    % Loop for WEST elements without vertices: 2 different conditions

    for i=1:ny-1
        for j=1:nx
            distance = sqrt((Xfc(i+1,j+1)-Xfc(i+1,j))^2 + (Yfc(i+1,j+1)-Yfc(i+1,j))^2);
            fluxCb (i,j) = (h_inf * Afe(i,j) * gamma_fe(i,j) * mod_Ee(i,j)/distance)/(h_inf * Afe(i,j) + gamma_fe(i,j) * mod_Ee(i,j)/distance);
            fluxVb (i,j) = - fluxCb(i,j) * T_inf - (h_inf * Afe(i,j) * gamma_fe(i,j) * (grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)))/(h_inf * Afe(i,j) + gamma_fe(i,j) * mod_Ee(i,j)/distance);
        end
    end

    for j=2
        for i=3:ny-1

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aS(i,j)) + fluxCb(i-1,j);
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - fluxVb(i-1,j)*T_inf;
        end
    end


    % Loop over SOUTH elements without vertices: DIRICHLET T = 400 K

    for i=2
        for j=3:nx-1

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 400*aS(i,j);
            aS(i,j) = 0;

        end
    end

    % Loop over NORTH elements without vertices: Dirichlet T = 300 K

    for i=ny
        for j=3:nx-1

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 300*aN(i,j);
            aN(i,j) = 0;


        end
    end

    % Loop over EAST elements without vertices: Zero flux condition

    for j=nx
        for i=3:ny-1

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho;
            aE(i,j) = 0;
        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 500*aW(i,j) - 400*aS(i,j);
            aW(i,j) = 0;
            aS(i,j) = 0;

        end
    end


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 300*aN(i,j);
            aE(i,j) = 0;
            aN(i,j) = 0;


        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 400*aW(i,j);
            aE(i,j) = 0;
            aS(i,j) = 0;
        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - NORTH CORNER

    for j=2
        for i=ny

            Non_Ortho = k*(grad_Txe(i-1,j)*Txe(i-1,j)+grad_Tye(i-1,j)*Tye(i-1,j)) + k*(grad_Txe(i-1,j-1)*Txe(i-1,j-1) + grad_Tye(i-1,j-1)*Tye(i-1,j-1)) + k*(grad_Txn(i,j-1)*Txn(i,j-1)+grad_Tyn(i,j-1)*Tyn(i,j-1)) + k*(grad_Txn(i-1,j-1)*Txn(i-1,j-1)+grad_Tyn(i-1,j-1)*Tyn(i-1,j-1));
            aE(i,j) = - k * g_Diffe(i-1,j) - max(-mdote(i-1,j),0);
            aW(i,j) = - k * g_Diffe(i-1,j-1) - max(mdote(i-1,j-1),0);
            aN(i,j) = - k * g_Diffn(i,j-1) - max(-mdotn(i,j-1),0);
            aS(i,j) = - k * g_Diffn(i-1,j-1) - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aS(i,j)) + fluxCb(i-1,j);
            bC(i,j) = Sc(i,j) * A(i-1,j-1) + Non_Ortho - 300*aN(i,j) - fluxVb(i-1,j)*T_inf;
            aW(i,j) = 0;
            aN(i,j) = 0;
        end
    end
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % SOlVER: Gauss Seidel


    for i=2:ny
        for j=2:nx

            Told(i,j) = T(i,j);

            T(i,j) = T(i,j) + lambda*((bC(i,j) - aE(i,j)*T(i,j+1) - aN(i,j)*T(i+1,j) - aW(i,j)*T(i,j-1) - aS(i,j)*T(i-1,j))/aC(i,j) - T(i,j));

            error = abs(T(i,j) - Told(i,j));


        end
    end


    iter = iter+1;
end








