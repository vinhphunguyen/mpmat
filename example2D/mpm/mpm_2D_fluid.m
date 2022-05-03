
% MPM for 2D


close all
clear all
clc

addpath('A:\Numerical Study\MPM\mpm codes-w\MPM-Matlab\fem_util')
addpath('A:\Numerical Study\MPM\mpm codes-w\MPM-Matlab\fem-functions')



%% Problem definition
% Geometry
dim = 2;

% Material properties
E   = 1000;
nu  = 0.3;
G   = E/2/(1+nu);
rho = 1000;
g   = 9.8;

K_f     = 2e6;
lambda  = 0.001;


stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

Problem = 'DamBreak';

switch lower(Problem)
    case 'blocknbeam'
        
        %% Material point of a block and a beam
        
        nnx     = 17;
        nny     = 17;
        node1   = square_node_array([0.6 0.6],[1.0 0.6],[1.0 1.0],[0.6 1.0],nnx,nny);
        inc_x   = 1;
        inc_y   = nnx;
        node_pa = [1 2 nnx+2 nnx+1];
        elem1   = make_elem(node_pa,nnx-1,nny-1,inc_x,inc_y);
        
        nnx     = 65;
        nny     = 9;
        node2   = square_node_array([0 0.2],[1.6 0.2],[1.6 0.4],[0 0.4],nnx,nny);
        inc_x   = 1;
        inc_y   = nnx;
        node_pa = [1 2 nnx+2 nnx+1];
        elem2   = make_elem(node_pa,nnx-1,nny-1,inc_x,inc_y);
        
        
        % plot_mesh(node ,elem ,'Q4','k-',1.);
        % plot_mesh(node1,elem1,'Q4','k-',1.);
        
        
        NoP_1 = size(elem1,1);
        NoP_2 = size(elem2,1);
        NoP   = NoP_1 + NoP_2;
        
        Xp    = zeros(NoP,dim);    % Position
        Vp    = zeros(NoP,dim);    % Velocity
        Mp    = zeros(NoP,1);      % Mass
        Volp  = zeros(NoP,1);      % Volume
        Volp0 = zeros(NoP,1);      % Initial volume
        sigp  = zeros(NoP,3);      % Stress
        epsp  = zeros(NoP,3);      % Strain
        Momp  = zeros(NoP,1);      % Momentum
        rhop  = zeros(NoP,1);      % Density
        Fp    = zeros(NoP,4);      % Gradient deformation
        type  = ones(NoP,1); 
        
        v0    = 0.1;
        
        for ie = 1 : NoP_1
            coord       = node1(elem1(ie,:),:);
            a           = (coord(3,1)-coord(1,1))*(coord(3,2)-coord(1,2));
            Volp(ie)    = a;
            Mp(ie)      = a*rho;
            Xp(ie,:)    = mean(coord);
            Vp(ie,:)    = [0 -v0];
            
            rhop(ie)    = rho;
            Fp(ie,:)    = [1 0 0 1];
        end
        
        for ie = 1 : NoP_2
            coord           = node2(elem2(ie,:),:);
            a               = (coord(3,1)-coord(1,1))*(coord(3,2)-coord(1,2));
            Volp(ie+NoP_1)  = a;
            Mp(ie+NoP_1)    = a*rho;
            Xp(ie+NoP_1,:)  = mean(coord);
            Vp(ie+NoP_1,:)  = [0 0];
            
            rhop(ie+NoP_1)  = rho;
            Fp(ie+NoP_1,:)  = [1 0 0 1];
        end
        
        Volp0   = Volp;
        Xp(:,2) = Xp(:,2)+0.2;
        
    case 'dambreak'
  
        %% Material point of a water column
        
        nnx     = 17;
        nny     = 41;
        node    = square_node_array([0 0],[0.4 0],[0.4 1.0],[0 1.0],nnx,nny);
        inc_x   = 1;
        inc_y   = nnx;
        node_pa = [1 2 nnx+2 nnx+1];
        elem    = make_elem(node_pa,nnx-1,nny-1,inc_x,inc_y);
        
        % Height  = 1;
        % B		= 200 * rho * g * Height / 7;
        
        % plot_mesh(node ,elem ,'Q4','k-',1.);
        % plot_mesh(node1,elem1,'Q4','k-',1.);
        
        
        NoP   = size(elem,1);
        
        Xp    = zeros(NoP,dim);    % Position
        Vp    = zeros(NoP,dim);    % Velocity
        Mp    = zeros(NoP,1);      % Mass
        Volp  = zeros(NoP,1);      % Volume
        Volp0 = zeros(NoP,1);      % Initial volume
        sigp  = zeros(NoP,3);      % Stress
        epsp  = zeros(NoP,3);      % Strain
        Momp  = zeros(NoP,1);      % Momentum
        rhop  = zeros(NoP,1);      % Density
        Fp    = zeros(NoP,4);      % Gradient deformation
        IntE  = zeros(NoP,1);      % Internal energy
        type  = ones(NoP,1); 
                
        for ie = 1 : NoP
            coord       = node(elem(ie,:),:);
            a           = (coord(3,1)-coord(1,1))*(coord(3,2)-coord(1,2));
            Volp(ie)    = a;
            Xp(ie,:)    = mean(coord);
            Vp(ie,:)    = [0 0];
            
            rhop(ie)    = rho;
            % rhop(ie)    = rho * (1 + (9810 / B * (Height-Xp(ie,2))))^(1/7);
            IntE(ie)    = 1;
            Mp(ie)      = a*rhop(ie);
            Fp(ie,:)    = [1 0 0 1];
        end
        
        Volp0   = Volp;
        type(:) = 2;

end



%% Computational grid
L = 1.6;

nex = 16;      % number of elements along X direction
ney = 16;      % number of elements along Y direction

dx = L/nex;
dy = L/ney;

inc_u=1;
inc_v=nex+1;
node_pattern=[ 1 2 nex+3 nex+2 ];

nodes    = square_node_array([0 0],[L 0],[L L],[0 L],nex+1,ney+1);
elements = make_elem(node_pattern,nex,ney,inc_u,inc_v);

NoE = size(elements,1);
NoN = size(nodes,1);

L_Nodes = find(nodes(:,1)==0);
R_Nodes = find(nodes(:,1)==L);
T_Nodes = find(nodes(:,2)==L);
B_Nodes = find(nodes(:,2)==0);

% figure; hold on;
% plot(Xp(:,1),Xp(:,2),'b.','markersize',10);
% % for i = 1 : NoP
% %     text(Xp(i,1),Xp(i,2),num2str(i));
% % end
% plot_mesh(nodes ,elements ,'Q4','k-',1.);
% % for i = 1 : NoN
% %     text(nodes(i,1),nodes(i,2),num2str(i));
% % end

M       = zeros(NoN,1);
Mom     = zeros(NoN,dim);
fint    = zeros(NoN,dim);
fext    = zeros(NoN,dim);
force   = zeros(NoN,dim);
V       = zeros(NoN,dim);

POS = [150 150 800 800];

clc



%% Structure to store material points for each element
% Two data structures are used
% 1. particle -> element
% 2. element  -> particles
element_of_point = ones(NoP,1);
point_in_element = cell(NoE,1);

for ip = 1 : NoP
    x = Xp(ip,1);
    y = Xp(ip,2);
    e = floor(x/dx) + 1 + ney*floor(y/dy);
    element_of_point(ip) = e;
end

for ie = 1 : NoE
    id = find(element_of_point == ie);
    point_in_element{ie} = id;
end



%% Calculation cycle
etime  = 0.5;

dt  = 0.001;
TOL = 1e-8;     % mass tolerance

steps = floor(etime/dt);
pstep = 50;
sav   = floor(steps/pstep);
isav  = 0;

particle.time  = cell(sav,1);
particle.pos   = cell(sav,1);
particle.vel   = cell(sav,1);
particle.sig   = cell(sav,1);
particle.eps   = cell(sav,1);
mesh.node      = nodes;
mesh.element   = elements;



for tstep = 1 : steps
    
    M     = zeros(NoN,1);
    Mom   = zeros(NoN,dim);
    fint  = zeros(NoN,dim);
    fext  = zeros(NoN,dim);

    
    
%% Phase 1 : Material points to Nodes
    % Loop over computational cells
    
    for ie = 1 : NoE
        element = elements(ie,:);
        node    = nodes(element,:);
        pie     = point_in_element{ie};
        
        % Loop over particles in current computational cell
        for ip = 1 : length(pie)
            pid = pie(ip);

            % N and B calculation
            pt(1)       = (2*Xp(pid,1)-(node(1,1)+node(2,1)))/(node(2,1)-node(1,1));
            pt(2)       = (2*Xp(pid,2)-(node(2,2)+node(3,2)))/(node(3,2)-node(2,2));
           
            [N,dNdxi]   = lagrange_basis('Q4',pt);	% element shape functions
            J0          = node'*dNdxi;              % element Jacobian matrix
            invJ0       = inv(J0);
            dNdx        = dNdxi*invJ0;
            
            % Mass, Momentum and Internal/External forces
            stress      = sigp(pid,:);
            
            for i = 1 : length(element)
                id          = element(i);
                dNIdx       = dNdx(i,1);
                dNIdy       = dNdx(i,2);
                M(id)       = M(id)      + N(i)*Mp(pid);
                Mom(id,:)   = Mom(id,:)  + N(i)*Mp(pid)*Vp(pid,:);
                fint(id,1)  = fint(id,1) - Volp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy);
                fint(id,2)  = fint(id,2) - Volp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy);
                fext(id,1)  = 0;
                fext(id,2)  = fext(id,2) - N(i)*Mp(pid)*g;
            end  
            
        end
      
    end
    

  
%% Phase 2 - Lagrangian : Solve Momentum equation

    
    % Residual force
    force = fint + fext;
    
        % BC: n/a in this problem
    Mom ([L_Nodes R_Nodes],1) = 0;
    force([L_Nodes R_Nodes],1) = 0;
    Mom ([B_Nodes],2) = 0;
    force([B_Nodes],2) = 0;
    
    % Update of Momentum n Velocity
    %for i = 1 : NoN
        Mom(:) = Mom(:) + dt * force(:);
    %end    
    
    % BC: n/a in this problem
    
    
%% Phase 3 - Convective : Nodes to Material points   
    eKinp = 0;
    ePotp = 0;
    for ie = 1 : NoE
        element = elements(ie,:);
        node    = nodes(element,:);
        pie     = point_in_element{ie};
        
        % Loop over particles in current computational cell
        for ip = 1 : length(pie)
            pid = pie(ip);

            % N and B calculation
            pt(1)       = (2*Xp(pid,1)-(node(1,1)+node(2,1)))/(node(2,1)-node(1,1));
            pt(2)       = (2*Xp(pid,2)-(node(2,2)+node(3,2)))/(node(3,2)-node(2,2));
           
            [N,dNdxi]   = lagrange_basis('Q4',pt);      % element shape functions
            J0          = node'*dNdxi;                  % element Jacobian matrix
            invJ0       = inv(J0);
            dNdx        = dNdxi*invJ0;
            
            Lp = zeros(2,2);    % particle gradient velocity 
            
            for i = 1 : length(element)
                id = element(i);     
                vI = [0 0];
                if M(id) > TOL
                    Vp(pid,:) = Vp(pid,:) + dt * N(i)*force(id,:)/M(id);
                    Xp(pid,:) = Xp(pid,:) + dt * N(i)*Mom(id,:)  /M(id);
                    vI        = Mom(id,:)/M(id);  
                end                                                
                Lp = Lp + vI'*dNdx(i,:);         
            end
            
            % Gradient deformation and Volume
            F         = ([1 0;0 1] + Lp*dt) * reshape(Fp(pid,:),2,2);
            Fp(pid,:) = reshape(F,1,4);
            Volp(pid) = det(F)*Volp0(pid);
            
            % Strain
            deps        = dt * 0.5 * (Lp+Lp'); 
            epsp(pid,:) = epsp(pid,:) + [deps(1,1) deps(2,2) 2*deps(1,2)];
            
            
            % Stress
            if ( type(pid) == 1 )
                dsigma      = C * [deps(1,1); deps(2,2); 2*deps(1,2)] ;
                sigp(pid,:) = sigp(pid,:) + dsigma';
            else
%                 J = det(F);
%                 p = K_f*log(J)/J;
                
%                 depsv = ( deps(1,1) + deps(2,2) )/3;
%                 p     = -K_f*depsv;

%                 p = B * ((rhop(pid)/1000)^7) - 1;
                
                temp        = deps(1,1) + deps(2,2);
                rhop(pid)   = rhop(pid) / (1+temp);
                IntE(pid)   = IntE(pid) + sigp(pid,:)*[deps(1,1); deps(2,2); 2*deps(1,2)]/rhop(pid);
                p           = (1.4-1)*rhop(pid)*IntE(pid);
                sigp(pid,:) = -p*[1; 1; 0] + 2 * lambda * [deps(1,1); deps(2,2); 2*deps(1,2)] ...
                                           - 2/3 * lambda * temp*[1; 1; 0];
            end

            
            % Energy
            eKinp = eKinp + 0.5*Mp(pid)* ( Vp(pid,1)^2 + Vp(pid,2)^2 );
            ePotp = ePotp + 0.5*Volp(pid)*sigp(pid,:)*epsp(pid,:)';
                       
        end
      
    end
                    
  
    
%% Phase 4 : Redefine messh and advance to next time step
    % Update list of material points for each element
    for ip = 1 : NoP
        x = Xp(ip,1);
        y = Xp(ip,2);
        e = floor(x/dx) + 1 + ney*floor(y/dy);
        element_of_point(ip) = e;
    end
    
    for ie = 1 : NoE
        id = find(element_of_point == ie);
        point_in_element{ie} = id;
    end

    
    % Saving data
    if ( mod(tstep,pstep) == 0 )
        isav = isav + 1;
        particle.time{isav} = tstep*dt;
        particle.pos{isav}  = Xp;
        particle.vel{isav}  = Vp;
        particle.sig{isav}  = sigp;
        particle.eps{isav}  = epsp;
    end
    
    if ( mod(tstep,pstep) == 0 )
        figure('Position',POS); hold on;
        plot_mesh(nodes, elements, 'Q4', 'k-', 1.);
        plot(Xp(:,1),Xp(:,2),'b.','markersize',20);        
        text(0.2,1.4,[' Time : ',num2str(tstep*dt),' s'],'FontSize',20)
    end

    if max(max(abs(Xp))) > 2
        disp('Explode')
        save('mpmSolid','particle','mesh')
        stop
    end
    
end

save('mpmSolid','particle','mesh')

% figure; hold on;
% plot(Xp(:,1),Xp(:,2),'b.','markersize',10);
% plot_mesh(nodes ,elements ,'Q4','k-',1.);



















