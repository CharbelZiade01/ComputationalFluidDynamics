clear all
close all
clc

% Maximum number of iterations
max_iter = 10000;
tolerance = 10^-15;
lambda = 0.1; %under relaxation factor
iter = 0;


T_inf = 300;
h_inf = 20;

% Define the vertices of the quadrilateral
vertices = [0, 0; 0.6, 0; 0.6, 0.3; 0, 0.3];

% Discretization parameters
ny = 4; % Number of points along y
nx = 7; % Number of points along x

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
%hold on
%scatter(Xfv,Yfv,'filled','bo');

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

hold on
scatter(Xfc,Yfc,'filled','ro');



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



%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% gDiff

g_Diffe = zeros(ny-1,nx);
g_Diffn = zeros(ny,nx-1);

for i=1:ny-1
    for j=1:nx

        g_Diffe(i,j) = Afe(i,j)/sqrt((Xfc(i+1,j+1) - Xfc(i+1,j))^2 + (Yfc(i+1,j+1) - Yfc(i+1,j))^2);

    end
end

for i=1:ny
    for j=1:nx-1

        g_Diffn(i,j) = Afn(i,j)/sqrt((Xfc(i+1,j+1) - Xfc(i,j+1))^2 + (Yfc(i+1,j+1) - Yfc(i,j+1))^2);

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
gamma_c = zeros (size(Xfc));
gamma_fe = zeros (ny-1,nx);
gamma_fn = zeros(ny,nx-1);

% Finding Gammas at centroids

for i=1:ny+1
    for j=1:2
        gamma_c(i,j) = 10^-3;
    end
end

for i=1:ny+1
    for j=3:nx+1
        gamma_c(i,j) = 10^2;
    end
end
    % for i=1:ny+1
    %     for j=1:nx+1
    %         %gamma_c(i,j) = T(i,j)* (Xfc(i,j)^2 + exp(0.1*Yfc(i,j)))/400;
    %         gamma_c(i,j) = 1;
    %     end
    % end

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % Finding Gammas at face centroids using interpolation

    %--------------------------
    % East


    for i=1:ny-1

        gamma_fe(i,1) = gamma_c(i+1,1);
        gamma_fe(i,nx) = gamma_c(i+1,nx+1);

        for j=2:nx-1
            gamma_fe(i,j) = gamma_c(i+1,j)*gamma_c(i+1,j+1)/((1 - ge(i,j))*gamma_c(i+1,j) + ge(i,j)*gamma_c(i+1,j+1));
        end
    end

    %-------------------------
    % North


    for j=1:nx-1

        gamma_fn(1,j) = gamma_c(1,j+1);
        gamma_fn(ny,j) = gamma_c(ny+1,j+1);

        for i=2:ny-1
            gamma_fn(i,j) = gamma_c(i,j+1)*gamma_c(i+1,j+1)/((1 - gn(i,j))*gamma_c(i,j+1) + gn(i,j)*gamma_c(i+1,j+1));
        end
    end
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Initial guess for the temperature

Told = zeros(ny+1,nx+1);
T = 400 * ones (ny+1,nx+1); % at the centroids and the boundary faces centroids
T(1,:) = 320;
T(ny+1,:) = 320;
T(:,nx+1) = 300;
T(:,1) = 300;
error = 900;

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% coefficients

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


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% HERE THE LOOP SHOULD START
% 
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
            Sc(i,j) = (T(i,j)* (2*Xfc(i,j) - 0.2*Yfc(i,j)))/400;
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

            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = 0;

        end
    end

    % Loop for WEST elements without vertices: DIRICHLET T = 320 K

    for j=2
        for i=3:ny-1
            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j) +  aN(i,j) + aS(i,j) + aW(i,j));
            bC(i,j) = 0 -320*aW(i,j);
            aW(i,j) = 0;
        end
    end


    % Loop over SOUTH elements without vertices: Von Neumman q = 100

    for i=2
        for j=3:nx-1
            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j) + aN(i,j) + aW(i,j));
            bC(i,j) = 0 - 100*Afn(i-1,j-1);
            aS(i,j) = 0;

        end
    end

    % Loop over NORTH elements without vertices: mixed

    for i=ny
        for j=3:nx-1
            Req = h_inf*gamma_fn(i,j-1)*Afn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2));
            Req = Req/(h_inf + gamma_fn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2)));
            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = 0;
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j)) + Req;
            bC(i,j) = 0 + Req*T_inf;
            

        end
    end

    % Loop over EAST elements without vertices: Zero Flux

        
    for j=nx
        for i=3:ny-1

            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j));
            bC(i,j) = 0 ;
            aE(i,j) = 0;

        end
    end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - SOUTH CORNER

    for j=2
        for i=2
            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j)  + aN(i,j) + aW(i,j));
            bC(i,j) = -320*aW(i,j);
            aW(i,j) = 0;
            aS(i,j) = 0;

        end
    end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - NORTH CORNER

    for j=nx
        for i=ny
            Req = h_inf*gamma_fn(i,j-1)*Afn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2));
            Req = Req/(h_inf + gamma_fn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2)));
            aE(i,j) = 0;
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = 0;
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aW(i,j) + aS(i,j))  + Req;
            bC(i,j) = 0 + Req*T_inf;
            
            
        end
    end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % EAST - SOUTH CORNER

    for j=nx
        for i=2

            aE(i,j) = 0;
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aW(i,j) + aN(i,j)) ;
            bC(i,j) = 0 - 100*Afn(i-1,j-1);
            aS(i,j) = 0;
        

        end
    end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    % WEST - NORTH CORNER

    for j=2
        for i=ny
            Req = h_inf*gamma_fn(i,j-1)*Afn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2));
            Req = Req/(h_inf + gamma_fn(i,j-1)/(sqrt((Xfc(i+1,j) - Xfc(i,j))^2 + (Yfc(i+1,j) - Yfc(i,j))^2)));
            aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
            aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
            aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
            aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
            aC(i,j) = -(aE(i,j)  + aS(i,j)+ aW(i,j)) + Req ;
            bC(i,j) = 0 + Req*T_inf - 320*aW(i,j);
            aN(i,j) = 0;
            aW(i,j) = 0;

        end
    end
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% SOlVER: Gauss Seidel



    for i=2:ny
        for j=2:nx

            Told(i,j) = T(i,j);

            T(i,j) = Told(i,j) + lambda*((bC(i,j) - aE(i,j)*Told(i,j+1) - aN(i,j)*Told(i+1,j) - aW(i,j)*Told(i,j-1) - aS(i,j)*Told(i-1,j))/aC(i,j) - Told(i,j));

            error = abs(T(i,j) - Told(i,j));


        end
    end


    iter = iter+1;
end

T

