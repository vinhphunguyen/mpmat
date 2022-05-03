function [aa] = hierarchicalGaussQuad(order,pts,nodel,mpts,pCoords,level,levelMax)
%
% Adaptive time integration (inspired by Finite Cell Method).
% Input:
%
% order: order of Gauss quadrature for a sub-cell
% pts  : coordinates of one cell (first level it is the element, next level
% it is one of sub-cells and so on)
% nodel: coordinates of one cell but it natural coordinates
% mpts : material points inside the element 
% pCoords: particle coordinates vector
% levelMax: maximum number of cell division.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 July 2015.

level = level + 1;

% first decompose the parent element [-1,1]x[-1,1]
% into 4 sub-cells

pt1 = pts(1,:);
pt2 = pts(2,:);
pt3 = pts(3,:);
pt4 = pts(4,:);
node    = square_node_array(pt1,pt2,pt3,pt4,3,3);
inc_u = 1;
inc_v = 3;
node_pattern = [ 1 2 3+2 3+1 ];
element      = make_elem(node_pattern,3-1,3-1,inc_u,inc_v);

pcoords = pCoords(mpts,:);

% hold on
% plot_mesh(node,element,'Q4','g.-',2);
% plot(Q(:,1),Q(:,2),'*');
% axis equal

% loop over sub-cells
Q = [];
W = [];
for e = 1:4
  sctr   = element(e,:); 
  coord  = node(sctr,:);
  coordl = nodel(sctr,:);
  in     = inpolygon(pcoords(:,1),pcoords(:,2),coord(:,1),coord(:,2));
  area   = (coordl(3,1) - coordl(1,1))*(coordl(3,2) - coordl(1,2));
  
  if ( any(in == 1) ) % sub-cell contains particles
    if level < levelMax
      [aa] = hierarchicalGaussQuad(order,coord,coordl,mpts,pCoords,level,levelMax);
    else
      [w,q] = quadrature(order,'GAUSS',2);      
      % transform quadrature points into the parent element      
      for n=1:length(w)
        N = lagrange_basis('Q4',q(n,:));
        Q = [Q;N'*coordl];
        W = [W;w(n)*area];
      end
    end 
  end  
end

aa.W = W;
aa.Q = Q;

% debug only

%plot(Q(:,1),Q(:,2),'*')




