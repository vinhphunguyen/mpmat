function [stress,strain] = updateStressHooke(dstrain,C,sigma,eps)

stress = sigma + dstrain*C;
strain = eps   + dstrain;