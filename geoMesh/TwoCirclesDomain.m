%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = TwoCirclesDomain(Demand,Arg)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox); 
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1   = dCircle(P,0.2,0.2,0.2);
  d2   = dCircle(P,0.8,0.8,0.2);
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
