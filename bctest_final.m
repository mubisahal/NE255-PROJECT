%% 2D One-Group Neutron Transpot Equation Solver for Top and Left vacuum boundaries and Bottom and Right Reflective boundaries
% Solves using Discrete Ordinate Method 

clc, clear all

%% Input
%Geometry
    xR = 1;     %horizontal length
    yT = 1;     %vertical length

    Nx = 1000; Ny = 1000;   %Mesh/Grid sizing

%Cross Sections
    sigt = 1;       %total xs
    sigs0 = 0.5;    %scattering xs

%Source
    S = 0.*ones(Nx,Ny);       % Change 0 to set source in every grid at once, homogeneous source
   
    % Remove the "%" ssymbol to set source coming from one or more faces
    
    S(:,Ny)=10;       %setting source coming from right face       
    S(Nx,:)=10;        %setting source coming from bottom face      
    %S(:,1)=10;         %setting source coming from left face 
    %S(1,:)=10;        %setting source coming from top face 

    Q = zeros(Nx,Ny); % Q set to zero. Can be modified using above matrix commands
    
% Calculation Parameters
    maxiter = 1E1; %Maximum number of iterations
    tol = 1e-8;    % setting required tolerance to stop the code from doing all iterations

%% Material Properies
siga = sigt - sigs0;

%% Geometry and Angular Discretization
xL = 0;
yB = 0;

Nang = 6;  % SN6 quadrature level

dx = (xR - xL)/(Nx);
dy = (yT - yB)/(Ny);

x = xL+dx/2:dx:xR-dx/2;
y = yB+dy/2:dy:yT-dy/2;

[Xarr,Yarr] = meshgrid(x,y);  %2d grid creation

%% Angular Descritization (SN6)
mu= [ 0.9261808, 0.6815076, 0.6815076, 0.2666355, 0.2666355, 0.2666355,
    -0.9261808, -0.6815076, -0.6815076, -0.2666355, -0.2666355, -0.2666355,
 	0.9261808, 0.6815076, 0.6815076, 0.2666355, 0.2666355, 0.2666355,
	-0.9261808, -0.6815076, -0.6815076, -0.2666355, -0.2666355, -0.2666355 ];
		
eta= [ 0.2666355, 0.6815076, 0.2666355, 0.9261808, 0.6815076, 0.2666355,
    0.2666355, 0.6815076, 0.2666355, 0.9261808, 0.6815076, 0.2666355,
	-0.2666355, -0.6815076, -0.2666355, -0.9261808, -0.6815076, -0.2666355,
	-0.2666355, -0.6815076, -0.2666355, -0.9261808, -0.6815076, -0.2666355 ];
	

wi= [ 0.1761263, 0.1572071, 0.1572071, 0.1761263, 0.1572071, 0.1761263,
    0.1761263, 0.1572071, 0.1761263, 0.1761263, 0.1572071, 0.1761263,
    0.1761263, 0.1572071, 0.1761263, 0.1761263, 0.1572071, 0.1761263,
	0.1761263, 0.1572071, 0.1761263, 0.1572071, 0.1572071, 0.1761263 ];
                         

%% Generation of Angular Flux Array
angular_flux = zeros(Nx,Ny,Nang*(Nang+2)/2);

%% Initialization of Angular Flux Evaluated at Half-Points Arrays
angular_flux_half_x = zeros(Nx,Nx+1,Nang*(Nang+2)/2);
angular_flux_half_y = zeros(Ny+1,Ny,Nang*(Nang+2)/2);

%% Transport Sweep
for iter = 1:maxiter            % initiation of sweep loop for each iteration
    
    
    angular_flux_prev = angular_flux;
    
    scalar_flux = zeros(Nx,Ny);
    
    for k=1:Nang*(Nang+2)/2
        scalar_flux = scalar_flux + 0.25.*wi(k).*angular_flux(:,:,k);
    end
        
    
    Q = ( S + sigs0.*scalar_flux );       
                
       angular_flux_half_x = zeros(Nx,Nx+1,Nang*(Nang+2)/2);   %setting the values as zero for conveniently
       angular_flux_half_y = zeros(Ny+1,Ny,Nang*(Nang+2)/2);   %applying bc during each sweep direction
 
  
    for L = 1:Nang*(Nang+2)/2
        
        if ( mu(L) > 0 && eta(L) > 0 )      % sweep direction condition
             
              
            for j = 1:Ny
                
                for i = 1:Nx   % calculating angular flux at each node i,j
                    
                    angular_flux(j,i,L) = ( 2*mu(L)*angular_flux_half_x(j,i,L)/dx + ...
                        2*eta(L)*angular_flux_half_y(j,i,L)/dy + Q(i,j) )/...
                        ( 2*mu(L)/dx + 2*eta(L)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1,L) = 2*angular_flux(j,i,L) - angular_flux_half_x(j,i,L);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j+1,m,L) = 2*angular_flux(j,m,L) - angular_flux_half_y(j,m,L);
                    
                end
                
            end
            
        elseif ( mu(L) < 0 && eta(L) > 0 )  % if condition being repeated for next sweep
           
                loc = find ( mu(L) == -mu & eta(L) == eta );  %setting bc 
                                                              
                angular_flux_half_x(:,:,L) = angular_flux_half_x(:,:,loc);
          
                % sweep is from right to left and top to bottom, so x value
                % is changed. the direction is different with respect to
                % how the matrix is laid out. so the sweep is not
                % technically from the said direction.
            
            for j = 1:Ny
                
                for i = Nx:-1:1
            
                    angular_flux(j,i,L) = ( -2*mu(L)*angular_flux_half_x(j,i+1,L)/dx + ...
                        2*eta(L)*angular_flux_half_y(j,i,L)/dy + Q(i,j) )/...
                        ( -2*mu(L)/dx + 2*eta(L)/dy + sigt );
                    
                    angular_flux_half_x(j,i,L) = 2*angular_flux(j,i,L) - angular_flux_half_x(j,i+1,L);
                    
                end
                
                for m = Nx:-1:1
                    
                    angular_flux_half_y(j+1,m,L) = 2*angular_flux(j,m,L) - angular_flux_half_y(j,m,L);
                    
                end
                
            end
                    
        elseif ( mu(L) > 0 && eta(L) < 0 )          
                
                loc = find( eta(L) == -eta & mu(L) == mu );  %setting bc
                
                angular_flux_half_y(:,:,L) = angular_flux_half_y(:,:,loc);
                                   
            
            for j = Ny:-1:1
                
                for i = 1:Nx
                    
                    angular_flux(j,i,L) = ( 2*mu(L)*angular_flux_half_x(j,i,L)/dx - ...
                        2*eta(L)*angular_flux_half_y(j+1,i,L)/dy + Q(i,j) )/...
                        ( 2*mu(L)/dx - 2*eta(L)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1,L) = 2*angular_flux(j,i,L) - angular_flux_half_x(j,i,L);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m,L) = 2*angular_flux(j,m,L) - angular_flux_half_y(j+1,m,L);
                    
                end
                
            end
            
        elseif ( mu(L) < 0 && eta(L) < 0 )   
            
                loc = find ( mu(L) == -mu & eta(L) == -eta );  %setting bc
                
                angular_flux_half_x(:,:,L) = angular_flux_half_x(:,:,loc);
                angular_flux_half_y(:,:,L) = angular_flux_half_y(:,:,loc);
            
            for j = Ny:-1:1
                
                for i = Nx:-1:1
                    
                    angular_flux(j,i,L) = ( -2*mu(L)*angular_flux_half_x(j,i+1,L)/dx - ...
                        2*eta(L)*angular_flux_half_y(j+1,i,L)/dy + Q(i,j) )/...
                        ( -2*mu(L)/dx - 2*eta(L)/dy + sigt );
                    
                    angular_flux_half_x(j,i,L) = 2*angular_flux(j,i,L) - angular_flux_half_x(j,i+1,L);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m,L) = 2*angular_flux(j,m,L) - angular_flux_half_y(j+1,m,L);
                    
                end
                
            end
                    
            
        end
        
    end
    
    scalar_flux_prev = zeros(Nx,Ny);
    scalar_flux_new = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux_prev = scalar_flux_prev + 0.25.*wi(k).*angular_flux_prev(:,:,k);
        scalar_flux_new = scalar_flux_new + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    residual = norm(scalar_flux_new-scalar_flux_prev); 
    err(iter) = residual;
    
    fprintf('Residual: %e     Iteration: %i \n',residual,iter-1);
    
    if ( residual < tol ) %condition to stop the code if the residue is less than tolerance
        
        break
        
    end 
  
       
end

%% results plotting

int = 1;
    
    figure(1)
    mesh(x,y,scalar_flux);
    xlabel('x-direction');
    ylabel('y-direction');
    zlabel('flux');
    axis([xL xR yB yT min(min(scalar_flux)) max(max(scalar_flux))])
    