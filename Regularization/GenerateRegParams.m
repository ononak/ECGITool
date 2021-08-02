function [ reg_param ] = GenerateRegParams(ftm , npoints, R)
% GENERATEREGPARAMS Generates range og regularization parameters
% Input variables
%  ftm -> Forward transfer matrix.
%  npoint -> number of regularization parameter to be generated
%  R -> Regularization matrix (if it is not provided, R is assumed to be I matrix)
%
% Output variables
% reg_param -> Regularization parameter list

% Reference: Per Christian Hansen regularization toolbox
%
% Author: Onder Nazim Onak ononak@gmail.com

switch nargin
    case 2
        [~, sm, ~] = svd(ftm);
        s =diag(sm);
    case 3
        [~,~,~,C,S] = gsvd(ftm,R);
        s =[diag(C); diag(S)];        
    otherwise
        error('invalid number of input');
end


smin_ratio = 16*eps;
  reg_param(npoints,1) = max([s(end),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  
  %reg_param = reg_param*12000/max(reg_param);

end

