clear all
close all
clc

% Maximum number of iterations
max_iter = 1000000000000;
tolerance = 10^-20;
lambda = 1;           % under relaxation factor
iter = 0;
tiny = 10^-8;

k = 15;
cp = 500;
ro = 7800;

a = 0.1;
b = 0.2;
T_inf = 300;
h_inf = 15;

t = 60;         % Time
dt = 1;         % Time Step

% Define the vertices of the quadrilateral
vertices = [0, 0; a, 0; a, b; 0 , b];

% Discretization parameters
ny = 41; % Number of points along y
nx = 21; % Number of points along x

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
%quiver(Xfv,Yfv,Exe,Eye,'b');
%quiver(Xfv,Yfv,Txe,Tye,'b');


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
% Initial guess for the temperature

Told = zeros(ny+1,nx+1,t+1);        % To store the values from previous iteration (Gauss Seidel)
T = 1000*ones (ny+1,nx+1,t+1);      % at the centroids and the boundary faces centroids

                                    % ny and nx are the number of nodes

% Boundary Conditions

T(1,:,:) = tiny;                    % Dirichlet
T(ny+1,:,:) = tiny;                 % Dirichlet
T(:,1,:) = T(:,2,:);                % Von Neumann
T(:,end,:) = T(:,end-1,:);          % Von Neumann

error = 9999;

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
gamma_c = zeros (size(Xfc));
gamma_fe = zeros (ny-1,nx);
gamma_fn = zeros(ny,nx-1);
% Tfe = zeros (ny-1,nx);
% Tfn = zeros(ny,nx-1);
Sc = zeros (size(Xfc));
grad_Tx = zeros(ny+1,nx+1);
grad_Ty = zeros(ny+1,nx+1);



%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% HERE THE LOOP SHOULD START
     

 while error>tolerance && iter<max_iter

      for hour = 2:dt:t+1

            
            % Boundary Conditions
            T(:,1,:) = T(:,2,:);                % Von Neumann
            T(:,end,:) = T(:,end-1,:);          % Von Neumann


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        

        % Finding Gammas at centroids

        % for i=1:ny+1
        %     for j=1:2
        %         gamma_c(i,j) = 10^-3;
        %     end
        % end
        %
        % for i=1:ny+1
        %     for j=3:nx+1
        %         gamma_c(i,j) = 10^2;
        %     end
        % end
        for i=1:ny+1
            for j=1:nx+1
                %gamma_c(i,j) = Told(i,j)* (Xfc(i,j)^2 + exp(0.1*Yfc(i,j)))/400;
                gamma_c(i,j) = k;
            end
        end

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
        % Temperature at the faces

        % East


        % for i=1:ny-1
        %
        %     Tfe(i,1) = T(i+1,1);
        %     Tfe(i,nx) = T(i+1,nx+1);
        %
        %     for j=2:nx-1
        %         Tfe(i,j) = (1 - ge(i,j))*T(i+1,j) + ge(i,j)*T(i+1,j+1);
        %     end
        % end
        %
        % % North
        %
        %
        % for j=1:nx-1
        %
        %     Tfn(1,j) = T(1,j+1);
        %     Tfn(ny,j) = T(ny+1,j+1);
        %
        %     for i=2:ny-1
        %         Tfn(i,j) = (1 - gn(i,j))*T(i,j+1) + gn(i,j)*T(i+1,j+1);
        %     end
        % end


%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % Source term at centroids


        for i=1:ny+1
            for j=1:nx+1
                Sc(i,j) = 0*(Told(i,j)* (2*Xfc(i,j) - 0.2*Yfc(i,j)))/400;

            end
        end


%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % Gradient of temperature at the centroids


        % for i=2:ny
        %     for j=2:nx
        %         grad_Tx(i,j) = (-Tfe(i-1,j-1)*Sxe(i-1,j-1) + Tfe(i-1,j)*Sxe(i-1,j) - Tfn(i-1,j-1)*Sxn(i-1,j-1) + Tfn(i,j-1)*Sxn(i,j-1))/A(i-1,j-1);
        %         grad_Ty(i,j) = (-Tfe(i-1,j-1)*Sye(i-1,j-1) + Tfe(i-1,j)*Sye(i-1,j) - Tfn(i-1,j-1)*Syn(i-1,j-1) + Tfn(i,j-1)*Syn(i,j-1))/A(i-1,j-1);
        %         %grad_Tx(i,j) = 1;
        %         %grad_Ty(i,j) = 1;
        %
        %     end
        % end
        %
        % % West Boundary
        % grad_Tx(:,1) = grad_Tx(:,2);
        % grad_Ty(:,1) = grad_Ty(:,2);
        %
        % % South Boundary
        % grad_Tx(1,:) = grad_Tx(2,:);
        % grad_Ty(1,:) = grad_Ty(2,:);
        %
        % % East Boundary
        % grad_Tx(:,nx+1) = grad_Tx(:,nx);
        % grad_Ty(:,nx+1) = grad_Ty(:,nx);
        %
        % % North Boundary
        % grad_Tx(ny+1,:) = 0;
        % grad_Ty(ny+1,:) = 0;


        %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % Gradient of temperature at the faces

        % for i=1:ny-1
        %
        %     grad_Txe(i,1) = grad_Tx(i+1,1);
        %     grad_Tye(i,1) = grad_Ty(i+1,1);
        %     grad_Txe(i,nx) = grad_Tx(i+1,nx+1);
        %     grad_Tye(i,nx) = grad_Ty(i+1,nx+1);
        %
        %
        %     for j=2:nx-1
        %
        %         grad_Txe(i,j) = (1 - ge(i,j))*grad_Tx(i+1,j) + ge(i,j)*grad_Tx(i+1,j+1);
        %         grad_Tye(i,j) = (1 - ge(i,j))*grad_Ty(i+1,j) + ge(i,j)*grad_Ty(i+1,j+1);
        %
        %     end
        % end
        %
        % for j=1:nx-1
        %
        %     grad_Txn(1,j) = grad_Tx(1,j+1);
        %     grad_Tyn(1,j) = grad_Ty(1,j+1);
        %     grad_Txn(ny,j) = grad_Tx(ny+1,j+1);
        %     grad_Tyn(ny,j) = grad_Ty(ny+1,j+1);
        %
        %     for i=2:ny-1
        %
        %         grad_Txn(i,j) = (1-gn(i,j))*grad_Tx(i,j+1) + gn(i,j)*grad_Tx(i+1,j+1);
        %         grad_Tyn(i,j) = (1-gn(i,j))*grad_Ty(i,j+1) + gn(i,j)*grad_Ty(i+1,j+1);
        %
        %     end
        % end

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % Loop for inner elements:

        for i=3:ny-1
            for j=3:nx-1
                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt;

            end
        end

        % Loop for WEST elements without vertices: Zero Flux

        for j=2
            for i=3:ny-1

                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aE(i,j) +  aN(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt;
                aW(i,j) = 0;
            end
        end


        % Loop over SOUTH elements without vertices: Dirichlet T = 0K

        for i=2
            for j=3:nx-1

                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aE(i,j) + aN(i,j) + aW(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = 0*Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aS(i,j)*tiny;
                aS(i,j) = 0;

            end
        end

        % Loop over NORTH elements without vertices: Dirichlet T = 0K

        for i=ny
            for j=3:nx-1

                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aE(i,j) + aW(i,j) + aS(i,j) + aN(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aN(i,j)*tiny;
                aN(i,j) = 0;


            end
        end

        % Loop over EAST elements without vertices: Zero Flux


        for j=nx
            for i=3:ny-1


                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt;
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
                aC(i,j) = -(aE(i,j)  + aN(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aS(i,j)*tiny;
                aW(i,j) = 0;
                aS(i,j) = 0;

            end
        end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % EAST - NORTH CORNER

        for j=nx
            for i=ny


                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aW(i,j) + aS(i,j) + aN(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aN(i,j)*tiny;
                aE(i,j) = 0;
                aN(i,j) = 0;


            end
        end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % EAST - SOUTH CORNER

        for j=nx
            for i=2

                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aW(i,j) + aN(i,j) + aS(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aS(i,j)*tiny;
                aS(i,j) = 0;
                aE(i,j) = 0;


            end
        end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % WEST - NORTH CORNER

        for j=2
            for i=ny

                aE(i,j) = - gamma_fe(i-1,j) * g_Diffe(i-1,j);
                aW(i,j) = - gamma_fe(i-1,j-1) * g_Diffe(i-1,j-1);
                aN(i,j) = - gamma_fn(i,j-1) * g_Diffn(i,j-1);
                aS(i,j) = - gamma_fn(i-1,j-1) * g_Diffn(i-1,j-1);
                aC(i,j) = -(aE(i,j)  + aS(i,j)+ aN(i,j)) + (ro * cp * A(i-1,j-1))/dt;
                bC(i,j) = Sc(i,j) * A(i-1,j-1) + (ro * cp * A(i-1,j-1))*T(i,j,hour-1)/dt - aN(i,j)*tiny;
                aN(i,j) = 0;
                aW(i,j) = 0;

            end
        end
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % SOlVER: Gauss Seidel



        for i=2:ny
            for j=2:nx

                Told(i,j,hour) = T(i,j,hour);

                T(i,j,hour) = Told(i,j,hour) + lambda*((bC(i,j) - aE(i,j)*Told(i,j+1,hour) - aN(i,j)*Told(i+1,j,hour) - aW(i,j)*Told(i,j-1,hour) - aS(i,j)*Told(i-1,j,hour))/aC(i,j) - Told(i,j,hour));

                error = abs(T(i,j,hour) - Told(i,j,hour));


            end
        end


        iter = iter+1;
     end

    
end

T

