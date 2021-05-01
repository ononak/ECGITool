
function [dpp] = PolyDerivative(pp)
%POLYDERIVATIVE SComputes the derivatie of the given polynomial
%
% Author: Onder Nazim Onak ononak@gmail.com
%
 ncoef = size(pp.coefs,1);
 ncol = size(pp.coefs,2) - 1;
 dcoefs = zeros(ncoef,ncol);
 
     for i=1:ncoef
         tmpcoef = polyder(pp.coefs(i,:)+1e-99);
         dcoefs(i,:) = tmpcoef;
     end
 
 dpp = mkpp(pp.breaks , dcoefs);

end