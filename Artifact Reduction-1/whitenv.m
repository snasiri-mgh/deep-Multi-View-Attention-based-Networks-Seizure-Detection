function [newVectors, whiteningMatrix, dewhiteningMatrix] = whitenv ...
    (vectors, E, D, s_verbose);

% Default value for 'verbose'
if nargin < 4, s_verbose = 'on'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

% ========================================================
% In some cases, rounding errors in Matlab cause negative
% eigenvalues (elements in the diagonal of D). Since it
% is difficult to know when this happens, it is difficult
% to correct it automatically. Therefore an error is 
% signalled and the correction is left to the user.
if any (diag (D) < 0),
  error (sprintf (['[ %d ] negative eigenvalues computed from the' ...
		   ' covariance matrix.\nThese are due to rounding' ...
		   ' errors in Matlab (the correct eigenvalues are\n' ...
		   'probably very small).\nTo correct the situation,' ...
		   ' please reduce the number of dimensions in the' ...
		   ' data\nby using the ''lastEig'' argument in' ...
		   ' function FASTICA, or ''Reduce dim.'' button\nin' ...
		   ' the graphical user interface.'], ...
		  sum (diag (D) < 0)));
end

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
whiteningMatrix = inv (sqrt (D)) * E';
dewhiteningMatrix = E * sqrt (D);

% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
if b_verbose, fprintf ('Whitening...\n'); end
newVectors =  whiteningMatrix * vectors;

% ========================================================
% Just some security...
if ~isreal(newVectors)
  error ('Whitened vectors have imaginary values.');
end

% Print some information to user
if b_verbose
  fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
    max (max (abs (cov (newVectors', 1) - eye (size (newVectors, 1))))));
end
