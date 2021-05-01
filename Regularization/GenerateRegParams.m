function [ reg_param ] = GenerateRegParams(ftm , npoints)
% GENERATEREGPARAMS Generates range og regularization parameters
% Input variables
%  ftm -> Forward transfer matrix.
%  npoint -> number of regularization parameter to be generated
%
% Output variables
% reg_param -> Regularization parameter list

% Reference: Per Christian Hansen regularization toolbox
%
% Author: Onder Nazim Onak ononak@gmail.com

[U sm V] = svd(ftm);
s =diag(sm);

smin_ratio = 16*eps;
  reg_param(npoints) = max([s(end),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  
  reg_param = reg_param*12000/max(reg_param);

end

