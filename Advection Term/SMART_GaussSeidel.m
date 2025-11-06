clear all
close all
clc

% Maximum number of iterations
max_iter = 1000;
tolerance = 10^-15;
lambda = 0.5; %under relaxation factor
iter = 0;


density = 1;

% Define the vertices of the quadrilateral
vertices = [0, 0; 1, 0; 1, 1; 0, 1];

% Discretization parameters
ny = 20; % Number of points along y
nx = 18; % Number of points along x

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

%--------------------------------------------------------------------------------------------------------------------------
% Centroids
Xc = zeros(ny-1,nx-1);
Yc = zeros(ny-1,nx-1);

for i=1:ny-1
    for j=1:nx-1
        Xc(i,j) = mean ([X(i,j), X(i,j+1),X(i+1,j),X(i+1,j+1)]);
        Yc(i,j) = mean ([Y(i,j), Y(i,j+1),Y(i+1,j),Y(i+1,j+1)]);
    end
end


%------------------------------------------------------------------------------------------------------------
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
%scatter(Xfv,Yfv,'filled','bo');

%------------------------------------------------------------------------------------------------------------
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

hold on
scatter(Xfc,Yfc,'filled','ro');



%------------------------------------------------------------------------------------------------------------
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



%------------------------------------------------------------------------------------------------------------------
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

%------------------------------------------------------------------------------------------------------------------
% Normal Vectors

Sxe = zeros(ny-1,nx);       % East Normal Vector coordinates
sxe = zeros(ny-1,nx);
Sye = zeros(ny-1,nx);
sye = zeros(ny-1,nx);

Sxn = zeros(ny,nx-1);       % North Normal Vector coordinates
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

%quiver(Xfv,Yfv,Sxe,Sye,'k');
%quiver(Xfh,Yfh,Sxn,Syn,'filled','k');

%--------------------------------------------------------------------------------------------------------------------
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
%quiver(Xfv,Yfv,Exe,Eye,'b');
%quiver(Xfv,Yfv,Txe,Tye,'b');


%-------------------------------------------------------------------------------------------------------------------
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



%------------------------------------------------------------------------------------------------------------------
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

%------------------------------------------------------------------------------------------------------------------
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

%------------------------------------------------------------------------------------------------------------------
% Velocity at Centroids

u = zeros(ny+1,nx+1);       % velocity in x-direction
v = zeros(ny+1,nx+1);       % velocity in y-direction

% NB: could have been done easier but to make to code more general

for i=1:ny+1
    for j=1:nx+1
        u(i,j) = sqrt(2)/2;
        v(i,j) = sqrt(2)/2;
    end
end


%------------------------------------------------------------------------------------------------------------------
% Velocity at horizontal (north) faces

un = zeros(ny,nx-1);       % velocity in x-direction
vn = zeros(ny,nx-1);       % velocity in y-direction

% NB: could have been done easier but to make to code more general

for i=1:ny
    for j=1:nx-1
        un(i,j) = (1 - gn(i,j))*u(i,j+1) + gn(i,j)*u(i+1,j+1);
        vn(i,j) = (1 - gn(i,j))*v(i,j+1) + gn(i,j)*v(i+1,j+1);
    end
end


%------------------------------------------------------------------------------------------------------------------
% Velocity at vertical (east) faces

ue = zeros(ny-1,nx);       % velocity in x-direction
ve = zeros(ny-1,nx);       % velocity in y-direction

% NB: could have been done easier but to make to code more general

for i=1:ny-1
    for j=1:nx
        ue(i,j) = (1 - ge(i,j))*u(i+1,j) + ge(i,j)*u(i+1,j+1);
        ve(i,j) = (1 - ge(i,j))*v(i+1,j) + ge(i,j)*v(i+1,j+1);

    end
end

%------------------------------------------------------------------------------------------------------------------
% Mass flow rate at horizontal (north) faces

mdotn = zeros(ny,nx-1);       % mass flow rate for north faces


for i=1:ny
    for j=1:nx-1
        mdotn(i,j) = density * (un(i,j)*Sxn(i,j) + vn(i,j)*Syn(i,j));
    end
end


%------------------------------------------------------------------------------------------------------------------
% Mass flow rate at vertical (east) faces

mdote = zeros(ny-1,nx);       % mass flow rate for east faces

for i=1:ny-1
    for j=1:nx
        mdote(i,j) = density *(ue(i,j)*Sxe(i,j) + ve(i,j)*Sye(i,j));
    end
end

%-------------------------------------------------------------------------------------------------------------------
% Initial guess for  phi
phife = zeros (ny-1,nx);
phifn = zeros(ny,nx-1);

aE = zeros(ny+1,nx+1);
aW = zeros(ny+1,nx+1);
aN = zeros(ny+1,nx+1);
aS = zeros(ny+1,nx+1);
aC = zeros(ny+1,nx+1);
bC = zeros(ny+1,nx+1);

% T = 400 * ones (ny+1,nx+1,2); % at the centroids and the boundary faces centroids
% T(1,:) = 320;
error = 900;
phi = zeros(ny+1,nx+1);
phi(:,1) = 1;
phi(1,:) = 0;
phiold = zeros(ny+1,nx+1);
smart = zeros(ny+1,nx+1);
phiCe = zeros(ny-1,nx);
phiCetild = zeros(ny-1,nx);
phifetild = zeros(ny-1,nx);
phifes = zeros(ny-1,nx);
phiDe = zeros(ny-1,nx);
phiUe = zeros(ny-1,nx);
phiCntild = zeros(ny,nx-1);
phifntild = zeros(ny,nx-1);
phifns = zeros(ny,nx-1);
phiCn= zeros(ny,nx-1);
phiDn= zeros(ny,nx-1);
phiUn= zeros(ny,nx-1);

%-------------------------------------------------------------------------------------------------------------------

% HERE THE LOOP SHOULD START

while error>tolerance && iter<max_iter


    %----------------------------------------------------------------------------------------------------------------------------
    % Phi at the faces

    % East



    for i=1:ny-1

        phife(i,1) = phi(i+1,1);
        phife(i,nx) = phi(i+1,nx+1);

        for j=2:nx-1
            phife(i,j) = (1 - ge(i,j))*phi(i+1,j) + ge(i,j)*phi(i+1,j+1);
        end
    end

    % North



    for j=1:nx-1

        phifn(1,j) = phi(1,j+1,1);
        phifn(ny,j) = phi(ny+1,j+1,1);

        for i=2:ny-1
            phifn(i,j) = (1 - gn(i,j))*phi(i,j+1) + gn(i,j)*phi(i+1,j+1);
        end
    end


    %-----------------------------------------------------------------------------------------------------------------------------
    % Source term at centroids

    %Sc = zeros (size(Xfc));


    %------------------------------------------------------------------------------------------------------------------------
    % Gradient of temperature at the centroids

    % grad_phix = zeros(ny+1,nx+1);
    % grad_phiy = zeros(ny+1,nx+1);

    % for i=2:ny
    %     for j=2:nx
    %         grad_phix(i,j) = (-phife(i-1,j-1)*Sxe(i-1,j-1) + phife(i-1,j)*Sxe(i-1,j) - phifn(i-1,j-1)*Sxn(i-1,j-1) + phifn(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);
    %         grad_phiy(i,j) = (-phife(i-1,j-1)*Sye(i-1,j-1) + phife(i-1,j)*Sye(i-1,j) - phifn(i-1,j-1)*Syn(i-1,j-1) + phifn(i,j-1)*Syn(i,j-1))/A(i-1,j-1);
    %
    %     end
    % end
    %
    % % West Boundary
    % grad_phix(:,1) = grad_phix(:,2);
    % grad_phiy(:,1) = grad_phiy(:,2);
    %
    % % South Boundary
    % grad_phix(1,:) = grad_phix(2,:);
    % grad_phiy(1,:) = grad_phiy(2,:);
    %
    % % East Boundary
    % grad_phix(:,nx+1) = 0;
    % grad_phiy(:,nx+1) = 0;
    %
    % % North Boundary
    % grad_phix(ny+1,:) = 0;
    % grad_phiy(ny+1,:) = 0;
    %
    %
    % %------------------------------------------------------------------------------------------------------------------------
    % % Gradient of phi at the faces
    %
    % grad_phixe = zeros(ny-1,nx);
    % grad_phiye = zeros(ny-1,nx);
    % grad_phixn = zeros(ny,nx-1);
    % grad_phiyn = zeros(ny,nx-1);
    %
    %
    %
    % for i=1:ny-1
    %
    %     grad_phixe(i,1) = grad_phix(i+1,1);
    %     grad_phiye(i,1) = grad_phiy(i+1,1);
    %     grad_phixe(i,nx) = grad_phix(i+1,nx+1);
    %     grad_phiye(i,nx) = grad_phiy(i+1,nx+1);
    %
    %
    %     for j=2:nx-1
    %
    %         grad_phixe(i,j) = (1 - ge(i,j))*grad_phix(i+1,j) + ge(i,j)*grad_phix(i+1,j+1);
    %         grad_phiye(i,j) = (1 - ge(i,j))*grad_phiy(i+1,j) + ge(i,j)*grad_phiy(i+1,j+1);
    %
    %     end
    % end
    %
    % for j=1:nx-1
    %
    %     grad_phixn(1,j) = grad_phix(1,j+1);
    %     grad_phiyn(1,j) = grad_phiy(1,j+1);
    %     grad_phixn(ny,j) = grad_phix(ny+1,j+1);
    %     grad_phiyn(ny,j) = grad_phiy(ny+1,j+1);
    %
    %     for i=2:ny-1
    %
    %         grad_phixn(i,j) = (1-gn(i,j))*grad_phix(i,j+1) + gn(i,j)*grad_phix(i+1,j+1);
    %         grad_phiyn(i,j) = (1-gn(i,j))*grad_phiy(i,j+1) + gn(i,j)*grad_phiy(i+1,j+1);
    %
    %     end
    % end

    %-------------------------------------------------------------------------------------------------------------
    % coefficients

for i=2:ny
    for j=2:nx

    smart(i,j) = (mdote(i-1,j)*(phifes(i-1,j) - phife(i-1,j)) - mdote(i-1,j-1)*(phifes(i-1,j-1) - phife(i-1,j-1)) + mdotn(i,j-1)*(phifns(i,j-1) - phifn(i,j-1)) - mdotn(i-1,j-1)*(phifns(i-1,j-1) - phifn(i-1,j-1)));

    end
end

    % Loop for inner elements:

    for i=3:ny-1
        for j=3:nx-1

            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = 0 - smart(i,j);

        end
    end

    % Loop for WEST elements without vertices: DIRICHLET T = 400 K

    for j=2
        for i=3:ny-1
            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = 0 - aW(i,j) - smart(i,j);
            aW(i,j) = 0;


        end
    end


    % Loop over SOUTH elements without vertices: DIRICHLET T = 320 K

    for i=2
        for j=3:nx-1
            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = 0 - 0*aS(i,j) - smart(i,j);
            aS(i,j) = 0;

        end
    end

    % Loop over NORTH elements without vertices: No FLUX CONDITION

    for i=ny
        for j=3:nx-1
            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = 0 - smart(i,j);
            aN(i,j) = 0;
            phi(i+1,j) = phi(i,j);

        end
    end

    % Loop over EAST elements without vertices:

    for j=nx
        for i=3:ny-1

            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = 0 - smart(i,j);
            aE(i,j) = 0;
            phi(i,j+1) = phi(i,j);


        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2
            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) =  -(aE(i,j) + aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = 0 - aW(i,j) - 0*aS(i,j) - smart(i,j);
            aW(i,j) = 0;
            aS(i,j) = 0;

        end
    end


    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny

            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aS(i,j));
            bC(i,j) = 0 - smart(i,j);
            aE(i,j) = 0;
            aN(i,j) = 0;
            phi(i+1,j) = phi(i,j);
            phi(i,j+1) = phi(i,j);


        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = 0 - 0*aS(i,j) - smart(i,j);
            aE(i,j) = 0;
            aS(i,j) = 0;
            phi(i,j+1) = phi(i,j);


        end
    end

    % WEST - NORTH CORNER

    for j=2
        for i=ny
            aE(i,j) =  - max(-mdote(i-1,j),0);
            aW(i,j) =  - max(mdote(i-1,j-1),0);
            aN(i,j) =  - max(-mdotn(i,j-1),0);
            aS(i,j) =  - max(mdotn(i-1,j-1),0);
            aC(i,j) = -(aE(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = 0 - aW(i,j) - smart(i,j);
            aW(i,j) = 0;
            aN(i,j) = 0;
            phi(i+1,j) = phi(i,j);

        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % SMART scheme implementation

    % Loop for internal faces: east

    for i=1:ny-1
        for j=2:nx-1
            if (mdote(i,j)>0)
                phiCe(i,j) = phi(i+1,j);
                phiDe(i,j) = phi(i+1,j+1);
                phiUe(i,j) = phi(i+1,j-1);
            else
                phiCe(i,j) = phi(i+1,j+1);
                phiDe(i,j) = phi(i+1,j);
                phiUe(i,j) = phi(i+1,j+2);
            end

            phiCetild(i,j) = (phiCe(i,j) - phiUe(i,j))/(phiDe(i,j) - phiUe(i,j));

            if phiCetild(i,j)<0 || phiCetild(i,j)>1
                phifetild(i,j) = phiCetild(i,j);
            elseif phiCetild(i,j)> 5/6 && phiCetild(i,j)<1
                phifetild(i,j) = 1;
            elseif  phiCetild(i,j)>= 0 && phiCetild(i,j)<= 5/6
                phifetild(i,j) = 3*phiCetild(i,j)/4 +3/8;      % QUICK scheme
            end


            phifes(i,j) = phifetild(i,j)*(phiDe(i,j) - phiUe(i,j)) + phiUe(i,j);
        end
    end

    % West Boundary

    for i=1:ny-1
        for j=1
            if (mdote(i,j)>0)
                phiCe(i,j) = phi(i+1,j);
                phiDe(i,j) = phi(i+1,j+1);
                phiUe(i,j) = phi(i+1,j);
            else
                phiCe(i,j) = phi(i+1,j);
                phiDe(i,j) = phi(i+1,j);
                phiUe(i,j) = phi(i+1,j+1);
            end

            phiCetild(i,j) = (phiCe(i,j) - phiUe(i,j))/(phiDe(i,j) - phiUe(i,j));

            if phiCetild(i,j)<0 || phiCetild(i,j)>1
                phifetild(i,j) = phiCetild(i,j);
            elseif phiCetild(i,j)> 5/6 && phiCetild(i,j)<1
                phifetild(i,j) = 1;
            elseif  phiCetild(i,j)>= 0 && phiCetild(i,j)<= 5/6
                phifetild(i,j) = 3*phiCetild(i,j)/4 +3/8;      % QUICK scheme
            end

            phifes(i,j) = phifetild(i,j)*(phiDe(i,j) - phiUe(i,j)) + phiUe(i,j);
        end
    end
    %
    % East Boundary

    for i=1:ny-1
        for j=nx
            if (mdote(i,j)>0)
                phiCe(i,j) = phi(i+1,j+1);
                phiDe(i,j) = phi(i+1,j+1);
                phiUe(i,j) = phi(i+1,j);
            else
                phiCe(i,j) = phi(i+1,j+1);
                phiDe(i,j) = phi(i+1,j);
                phiUe(i,j) = phi(i+1,j+1);
            end
            phiCetild(i,j) = (phiCe(i,j) - phiUe(i,j))/(phiDe(i,j) - phiUe(i,j));

            if phiCetild(i,j)<0 || phiCetild(i,j)>1
                phifetild(i,j) = phiCetild(i,j);
            elseif phiCetild(i,j)> 5/6 && phiCetild(i,j)<1
                phifetild(i,j) = 1;
            elseif  phiCetild(i,j)>= 0 && phiCetild(i,j)<= 5/6
                phifetild(i,j) = 3*phiCetild(i,j)/4 +3/8;      % QUICK scheme
            end

            phifes(i,j) = phifetild(i,j)*(phiDe(i,j) - phiUe(i,j)) + phiUe(i,j);
        end
    end


    % South Boundary

    for i=1
        for j=1:nx-1
            if (mdotn(i,j)>0)
                phiCn(i,j) = phi(i,j+1);
                phiDn(i,j) = phi(i+1,j+1);
                phiUn(i,j) = phi(i,j+1);
            else
                phiCn(i,j) = phi(i,j+1);
                phiUn(i,j) = phi(i+1,j+1);
                phiDn(i,j) = phi(i,j+1);
            end

            phiCntild(i,j) = (phiCn(i,j) - phiUn(i,j))/(phiDn(i,j) - phiUn(i,j));

            if phiCntild(i,j)<0 || phiCntild(i,j)>1
                phifntild(i,j) = phiCntild(i,j);
            elseif phiCntild(i,j)> 5/6 && phiCntild(i,j)<1
                phifntild(i,j) = 1;
            elseif  phiCntild(i,j)>= 0 && phiCntild(i,j)<= 5/6
                phifntild(i,j) = 3*phiCntild(i,j)/4 +3/8;      % QUICK scheme
            end

            phifns(i,j) = phifntild(i,j)*(phiDn(i,j) - phiUn(i,j)) + phiUn(i,j);
        end
    end

    % North Boundary

    for i=1
        for j=1:nx-1
            if (mdotn(i,j)>0)
                phiCn(i,j) = phi(i+1,j+1);
                phiDn(i,j) = phi(i+1,j+1);
                phiUn(i,j) = phi(i,j+1);
            else
                phiCn(i,j) = phi(i,j+1);
                phiDn(i,j) = phi(i+1,j+1);
                phiUn(i,j) = phi(i,j+1);
            end

            phiCntild(i,j) = (phiCn(i,j) - phiUn(i,j))/(phiDn(i,j) - phiUn(i,j));

            if phiCntild(i,j)<0 || phiCntild(i,j)>1
                phifntild(i,j) = phiCntild(i,j);
            elseif phiCntild(i,j)> 5/6 && phiCntild(i,j)<1
                phifntild(i,j) = 1;
            elseif  phiCntild(i,j)>= 0 && phiCntild(i,j)<= 5/6
                phifntild(i,j) = 3*phiCntild(i,j)/4 +3/8;      % QUICK scheme
            end

            phifns(i,j) = phifntild(i,j)*(phiDn(i,j) - phiUn(i,j)) + phiUn(i,j);
        end
    end
    % Loop for internal faces: north

    for i=2:ny-1
        for j=1:nx-1
            if(mdotn(i,j)>0)
                phiCn(i,j) = phi(i,j+1);
                phiDn(i,j) = phi(i+1,j+1);
                phiUn(i,j) = phi(i-1,j+1);
            else
                phiCn(i,j) = phi(i+1,j+1);
                phiDn(i,j) = phi(i,j+1);
                phiUn(i,j) = phi(i+2,j+1);
            end

            phiCntild(i,j) = (phiCn(i,j) - phiUn(i,j))/(phiDn(i,j) - phiUn(i,j));

            if phiCntild(i,j)<0 || phiCntild(i,j)>1
                phifntild(i,j) = phiCntild(i,j);
            elseif phiCntild(i,j)> 5/6 && phiCntild(i,j)<1
                phifntild(i,j) = 1;
            elseif  phiCntild(i,j)>= 0 && phiCntild(i,j)<= 5/6
                phifntild(i,j) = 3*phiCntild(i,j)/4 +3/8;      % QUICK scheme
            end

            phifns(i,j) = phifntild(i,j)*(phiDn(i,j) - phiUn(i,j)) + phiUn(i,j);
        end
    end

    %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Solver

    for i=2:ny
        for j=2:nx

            phiold(i,j) = phi(i,j);

            phi(i,j) = phi(i,j) + lambda*((bC(i,j) - aE(i,j)*phi(i,j+1) - aN(i,j)*phi(i+1,j) - aW(i,j)*phi(i,j-1) - aS(i,j)*phi(i-1,j))/aC(i,j) - phi(i,j));

            error = abs(phi(i,j) - phiold(i,j));

        end
    end


    iter = iter+1;
end


phi

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plotting

if mod(nx,2) == 0     % checking if nx is odd or even
    middle_column = phi(2:end-1 , (nx)/2 +1);
else
    middle_column = phi(2:end-1 , (nx+1)/2 +1);
end

row_number = 2:ny;

hold off
figure;

plot(row_number, middle_column);
xlabel('Element number along Y-axis');
ylabel('φ');

