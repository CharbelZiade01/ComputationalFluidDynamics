clear all
close all
clc

% Maximum number of iterations
max_iter = 1000000000;
tolerance = 10^-10;
lambda = 1; %under relaxation factor

T_inf = 300;
h_inf = 15;

% Define the vertices of the quadrilateral
vertices = [0, 0; 2, 0; 3, 5; -1, 3];

% Discretization parameters
ny = 6; % Number of points along y
nx = 5; % Number of points along x

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

%v1{ny-1,nx-1} = [0 0];
%v2{ny-1,nx-1} = [0 0];
%v3{ny-1,nx-1} = [0 0];
%v4{ny-1,nx-1} = [0 0];
p1 = zeros(ny-1,nx-1);
p2 = zeros(ny-1,nx-1);

for i=1:ny-1
    for j=1:nx-1
        %the area is the sum of two cross products over 2
        %v1 {i,j} = [X(i+1,j)-X(i,j) , Y(i+1,j)-Y(i,j)];
        %v2 {i,j} = [X(i,j+1)-X(i,j) , Y(i,j+1)-Y(i,j)];
        %v3 {i,j} = [X(i+1,j+1)-X(i,j+1) , Y(i+1,j+1)-Y(i,j+1)];
        %v4 {i,j} = [X(i+1,j+1)-X(i+1,j) , Y(i+1,j+1)-Y(i+1,j)];

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

Sxe = zeros(ny-1,nx); % East Normal Vector coordinates
sxe = zeros(ny-1,nx);
Sye = zeros(ny-1,nx);
sye = zeros(ny-1,nx);

Sxn = zeros(ny,nx-1); % North Normal Vector coordinates
sxn = zeros(ny,nx-1);
Syn = zeros(ny,nx-1);
syn = zeros(ny,nx-1);

% variables used to do a plot of the vectors
%Ue = zeros(ny-1,nx);
%Ve = zeros(ny-1,nx);
%Un = zeros(ny,nx-1);
%Vn = zeros(ny,nx-1);

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


%-------------------------------------------------------------------------------------------------------------------
% Initial guess for the temperature


T = 400 * ones (ny+1,nx+1,2); % at the centroids and the boundary faces centroids
T(1,:) = 320;
error = zeros(ny-1,nx-1);
iter = 0;

%-------------------------------------------------------------------------------------------------------------------

% HERE THE LOOP SHOULD START

for k = 1:max_iter
    %-------------------------------------------------------------------------------------------------------------------
    % Finding Gammas at centroids

    gamma_c = zeros (size(Xfc));

    for i=1:ny+1
        for j=1:nx+1
            gamma_c(i,j) = T(i,j,1)* (Xfc(i,j)^2 + exp(0.1*Yfc(i,j)))/400;
        end
    end

    %-------------------------------------------------------------------------------------------------------------------
    % Finding Gammas at face centroids using interpolation

    %--------------------------
    % East

    gamma_fe = zeros (ny-1,nx);

    for i=1:ny-1

        gamma_fe(i,1) = gamma_c(i,1);
        gamma_fe(i,nx) = gamma_c(i,nx+1);

        for j=2:nx-1
            gamma_fe(i,j) = (1 - ge(i,j))*gamma_c(i,j) + ge(i,j)*gamma_c(i,j+1);
        end
    end

    %-------------------------
    % North

    gamma_fn = zeros(ny,nx-1);

    for j=1:nx-1

        gamma_fn(1,j) = gamma_c(1,j);
        gamma_fn(ny,j) = gamma_c(ny+1,j);

        for i=2:ny-1
            gamma_fn(i,j) = (1 - gn(i,j))*gamma_c(i,j) + gn(i,j)*gamma_c(i+1,j);
        end
    end

    %----------------------------------------------------------------------------------------------------------------------------
    % Temperature at the faces

    % East

    Tfe = zeros (ny-1,nx);

    for i=1:ny-1

        Tfe(i,1) = T(i,1);
        Tfe(i,nx) = T(i,nx+1);

        for j=2:nx-1
            Tfe(i,j) = (1 - ge(i,j))*T(i,j,1) + ge(i,j)*T(i,j+1,1);
        end
    end

    % North

    Tfn = zeros(ny,nx-1);

    for j=1:nx-1

        Tfn(1,j) = T(1,j,1);
        Tfn(ny,j) = T(ny+1,j,1);

        for i=2:ny-1
            Tfn(i,j) = (1 - gn(i,j))*T(i,j,1) + gn(i,j)*T(i+1,j,1);
        end
    end


    %-----------------------------------------------------------------------------------------------------------------------------
    % Source term at centroids

    Sc = zeros (size(Xfc));

    for i=1:ny+1
        for j=1:nx+1
            Sc(i,j) = (T(i,j,1)* (2*Xfc(i,j) - 0.2*Yfc(i,j)))/400;
        end
    end

    %-------------------------------------------------------------------------------------------------------------------
    % Finding Source term at face centroids using interpolation

    % East

    Sfe = zeros (ny-1,nx);

    for i=1:ny-1

        Sfe(i,1) = Sc(i,1);
        Sfe(i,nx) = Sc(i,nx+1);

        for j=2:nx-1
            Sfe(i,j) = (1 - ge(i,j))*Sc(i,j) + ge(i,j)*Sc(i,j+1);
        end
    end

    % North

    Sfn = zeros(ny,nx-1);

    for j=1:nx-1

        Sfn(1,j) = Sc(1,j);
        Sfn(ny,j) = Sc(ny+1,j);

        for i=2:ny-1
            Sfn(i,j) = (1 - gn(i,j))*Sc(i,j) + gn(i,j)*Sc(i+1,j);
        end
    end

    %------------------------------------------------------------------------------------------------------------------------
    % Gradient of temperature at the centroids

    grad_Tx = zeros(ny+1,nx+1);
    grad_Ty = zeros(ny+1,nx+1);

    for i=2:ny
        for j=2:nx
            grad_Tx(i,j) = (-Tfe(i-1,j-1)*Sxe(i-1,j-1) + Tfe(i-1,j)*Sxe(i-1,j) - Tfn(i-1,j-1)*Sxn(i-1,j-1) + Tfn(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);
            grad_Ty(i,j) = (-Tfe(i-1,j-1)*Sye(i-1,j-1) + Tfe(i-1,j)*Sye(i-1,j) - Tfn(i-1,j-1)*Syn(i-1,j-1) + Tfn(i,j-1)*Syn(i,j-1))/A(i-1,j-1);

        end
    end


    %------------------------------------------------------------------------------------------------------------------------
    % Gradient of temperature at the faces

    grad_Txe = zeros(ny-1,nx);
    grad_Tye = zeros(ny-1,nx);
    grad_Txn = zeros(ny,nx-1);
    grad_Tyn = zeros(ny,nx-1);


    for i=1:ny-1
        for j=2:nx-1

            grad_Txe(i,j) = (1 - ge(i,j))*grad_Tx(i+1,j) + ge(i,j)*grad_Tx(i+1,j+1);
            grad_Tye(i,j) = (1 - ge(i,j))*grad_Ty(i+1,j) + ge(i,j)*grad_Ty(i+1,j+1);

        end
    end

    for i=2:ny-1
        for j=1:nx-1

            grad_Txn(i,j) = (1-gn(i,j))*grad_Tx(i,j+1) + gn(i,j)*grad_Tx(i+1,j+1);
            grad_Tyn(i,j) = (1-gn(i,j))*grad_Ty(i,j+1) + gn(i,j)*grad_Ty(i+1,j+1);

        end
    end

    %-------------------------------------------------------------------------------------------------------------
    % coefficients

    aE = zeros(ny-1,nx-1);
    aW = zeros(ny-1,nx-1);
    aN = zeros(ny-1,nx-1);
    aS = zeros(ny-1,nx-1);
    aC = zeros(ny-1,nx-1);
    bC = zeros(ny-1,nx-1);
    A_glob = zeros((ny-1)*(ny-1),(nx-1)*(nx-1));
    T_glob = zeros((ny-1)*(ny-1),1);
    B_glob = zeros((ny-1)*(ny-1),1);


    % Loop for inner elements:

    for i=2:ny-2
        for j=2:nx-2

            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = Sc(i+1,j+1) * A(i,j) + gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) + gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) + gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) + gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j));
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);

        end
    end

    % Loop for WEST elements without vertices:

    for j=1
        for i=2:ny-2
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aS(i,j) + aW(i,j)) + gamma_fe(i,j) * g_Diffe(i,j);
            bC(i,j) = Sc(i+1,j+1) * A(i,j) + gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) + gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) + gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) + gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j)) - aW(i,j)*400;
            aW(i,j) = 0;
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end


    % Loop over SOUTH elements without vertices:

    for i=1
        for j=2:nx-2
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j)+aS(i,j)) + gamma_fn(i,j) * g_Diffn(i,j);
            bC(i,j) = Sc(i+1,j+1) * A(i,j) - gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) + gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) + gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) + gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j)) - 320*aS(i,j) ;
            %aS(i,j) = 0;
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    % Loop over NORTH elements without vertices:

    for i=ny-1
        for j=2:nx-2
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j));
            bC(i,j) = Sc(i+1,j+1) * A(i,j) + gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) + gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) + gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) + gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j));
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    % Loop over EAST elements without vertices:

    Req = zeros(ny-1,nx);
    for j=nx-1
        for i=2:ny-2
            Req(i,j+1) = h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2)+h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1));
            Req(i,j+1) = Req(i,j+1)/(h_inf*Afe(i,j+1)+gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2));
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j)) + (Req(i,j+1)-h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1)))*T_inf;
            bC(i,j) = Sc(i+1,j+1) * A(i,j) - gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) - gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) - gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) - gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j)) + Req(i,j+1)*T_inf;
            aE(i,j) = 0;
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    aE(1,1) = - gamma_fe(1,2) * g_Diffe(1,2);
    aW(1,1) = - gamma_fe(1,1) * g_Diffe(1,1);
    aN(1,1) = - gamma_fn(2,1) * g_Diffn(2,1);
    aS(1,1) = - gamma_fn(1,1) * g_Diffn(1,1);
    aC(1,1) = -(aE(i,j)  + aN(i,j));
    bC(1,1) = Sc(2,2) * A(1,1) -(- gamma_fe(1,2)*(grad_Txe(1,2)*Txe(1,2)+grad_Tye(1,2)*Tye(1,2)) - gamma_fe(1,1)*(grad_Txe(1,1)*Txe(1,1) + grad_Tye(1,1)*Tye(1,1)) - gamma_fn(2,1)*(grad_Txn(2,1)*Txn(2,1)+grad_Tyn(2,1)*Tyn(2,1)) - gamma_fn(1,1)*(grad_Txn(1,1)*Txn(1,1)+grad_Tyn(1,1)*Tyn(1,1))) - 400*aW(1,1) - 320*aS(1,1);

   % for i=1
    %    for j=1

            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);

     %   end
    %end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx-1
        for i=ny-1
            Req(i,j+1) = h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2) + h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1));
            Req(i,j+1) = Req(i,j+1)/(h_inf*Afe(i,j+1)+gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2));
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1) - (h_inf*Afe(i,nx)*gamma_fe(i,nx)*mod_Ee(i,nx)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2))*T_inf/(h_inf*Afe(i,nx)+gamma_fe(i,nx)*mod_Ee(i,nx)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2));
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aW(i,j) + aS(i,j)) + (Req(i,j+1)-h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1)))*T_inf;
            bC(i,j) = Sc(i+1,j+1) * A(i,j) - gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) - gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) - gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) - gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j)) + Req(i,j+1)*T_inf;
          
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx-1
        for i=1
            Req(i,j+1) = h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2) + h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1));
            Req(i,j+1) = Req(i,j+1)/(h_inf*Afe(i,j+1)+gamma_fe(i,j+1)*mod_Ee(i,j+1)/sqrt((Xfc(i+1,j+2)-Xfc(i+1,j+1))^2 + (Yfc(i+1,j+2)-Yfc(i+1,j+1))^2));
            aE(1,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(1,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(1,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(1,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(1,j) = -(aW(i,j) + aN(i,j)) + (Req(i,j+1)-h_inf*Afe(i,j+1)*gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1) + grad_Tye(i,j+1)*Tye(i,j+1)))*T_inf;
            bC(1,j) = Sc(i+1,j+1) * A(i,j) - gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) - gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) - gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) - gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j) + grad_Tyn(i,j)*Tyn(i,j)) + Req(i,j+1)*T_inf - 320*aS(i,j);
            
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    % WEST - NORTH CORNER

    for i=ny-1
        for j=1
            aE(i,j) = - gamma_fe(i,j+1) * g_Diffe(i,j+1);
            aW(i,j) = - gamma_fe(i,j) * g_Diffe(i,j);
            aN(i,j) = - gamma_fn(i+1,j) * g_Diffn(i+1,j);
            aS(i,j) = - gamma_fn(i,j) * g_Diffn(i,j);
            aC(i,j) = -(aE(i,j) + aS(i,j));
            bC(i,j) = Sc(i+1,j+1) * A(i,j) - gamma_fe(i,j+1)*(grad_Txe(i,j+1)*Txe(i,j+1)+grad_Tye(i,j+1)*Tye(i,j+1)) - gamma_fe(i,j)*(grad_Txe(i,j)*Txe(i,j) + grad_Tye(i,j)*Tye(i,j)) - gamma_fn(i+1,j)*(grad_Txn(i+1,j)*Txn(i+1,j)+grad_Tyn(i+1,j)*Tyn(i+1,j)) - gamma_fn(i,j)*(grad_Txn(i,j)*Txn(i,j)+grad_Tyn(i,j)*Tyn(i,j)) - 400*aW(i,j);
            
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j) = aC(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j-1) = aW(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-1)*(nx-1)+j+1) = aE(i,j);
            %A_glob((i-1)*(nx-1)+j,(i-2)*(nx-1)+j) = aS(i,j);
            %A_glob((i-1)*(nx-1)+j,(i)*(nx-1)+j) = aN(i,j);
            %B_glob((i-1)*(nx-1)+j) = bC(i,j);
        end
    end

    for i=1:ny-1
        for j=1:nx-1
            
            T(i+1,j+1,2) =T(i+1,j+1,1) + lambda*((bC(i,j) - aE(i,j)*T(i+1,j+2,1) - aN(i,j)*T(i+2,j+1,1) - aW(i,j)*T(i+1,j,1) - aS(i,j)*T(i,j+1,1))/aC(i,j) - T(i+1,j+1,1));
           
            error(i,j) = abs(T(i+1,j+1,1) - T(i+1,j+1,2));
            if(error < tolerance)
                break
            end
            T(i+1,j+1,1) = T(i+1,j+1,2);
        end
    end

    if(error < tolerance)
        break
    end

end

T(2:ny,2:nx)
