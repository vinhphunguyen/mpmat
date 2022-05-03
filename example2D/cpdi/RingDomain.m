%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = RingDomain(Demand,Arg)
  BdBox = [-2 2 -2 2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox); 
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1   = dCircle(P,0.0,0.0,0.75);
  d2   = dCircle(P,0.0,0.0,1.25);
  Dist = dDiff(d2,d1);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
