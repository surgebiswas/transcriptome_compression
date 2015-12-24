function [Val, Ind1, Ind2]=max2d(Mat)
%-------------------------------------------------------------------------- 
% max2d function      2d maximum function. Return the maximum value 
%                   and index in a 2d matrix. 
% Input  : - Matrix. 
% Output : - Maximum value. 
%          - [I,J] index of maximum. 
% Tested : Matlab 5.3 
%     By : Eran O. Ofek                  October 2000 
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%-------------------------------------------------------------------------- 

[V1,I1] = max(Mat); 
[V2,I2] = max(V1); 

Val = V2; 
Ind1 = I1(I2);   % I 
Ind2 = I2;       % J 
end