function [Aes, Wes, J]= fasticaamar(mixedsig, varargin)
disp('FAST')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some basic requirements of the data
if nargin == 0,
  error ('You must supply the mixed data as input argument.');
end

if length (size (mixedsig)) > 2,
  error ('Input data can not have more than two dimensions.');
end

if any (any (isnan (mixedsig))),
  error ('Input data contains NaN''s.');
end

if ~isa (mixedsig, 'double')
  fprintf ('Warning: converting input data into regular (double) precision.\n');
  mixedsig = double (mixedsig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = 3;
g                 = 'tanh';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 0;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch lower (varargin{i})
     case 'stabilization'
      stabilization = lower (varargin{i+1});
     case 'maxfinetune'
      maxFinetune = varargin{i+1};
     case 'samplesize'
      sampleSize = varargin{i+1};
     case 'verbose'
      verbose = lower (varargin{i+1});
      % silence this program also
      if strcmp (verbose, 'off'), b_verbose = 0; end
     case 'firsteig'
      firstEig = varargin{i+1};
     case 'lasteig'
      lastEig = varargin{i+1};
     case 'interactivepca'
      interactivePCA = lower (varargin{i+1});
     case 'approach'
      approach = lower (varargin{i+1});
     case 'numofic'
      numOfIC = varargin{i+1};
      % User has supplied new value for numOfIC.
      % We'll use this information later on...
      userNumOfIC = 1;
     case 'g'
      g = lower (varargin{i+1});
     case 'finetune'
      finetune = lower (varargin{i+1});
     case 'a1'
      a1 = varargin{i+1};
     case 'a2'
      a2 = varargin{i+1};
     case {'mu', 'myy'}
      myy = varargin{i+1};
     case 'epsilon'
      epsilon = varargin{i+1};
     case 'maxnumiterations'
      maxNumIterations = varargin{i+1};
     case 'initguess'
      % no use setting 'guess' if the 'initState' is not set
      initState = 'guess';
      guess = varargin{i+1};
     case 'displaymode'
      displayMode = lower (varargin{i+1});
     case 'displayinterval'
      displayInterval = varargin{i+1};
     case 'pcae'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      E = varargin{i+1};
     case 'pcad'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      D = varargin{i+1};
     case 'whitesig'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whitesig = varargin{i+1};
     case 'whitemat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = varargin{i+1};
     case 'dewhitemat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = varargin{i+1};
     case 'only'
      % if the user only wants to calculate PCA or...
      switch lower (varargin{i+1})
       case 'pca'
	only = 1;
       case 'white'
	only = 2;
       case 'all'
	only = 3;
      end
      
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print information about data
if b_verbose
  fprintf('Number of signals: %d\n', Dim);
  fprintf('Number of samples: %d\n', NumOfSampl);
end

% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
  if b_verbose,
    fprintf ('Whitened signal and corresponding matrises supplied.\n');
    fprintf ('PCA calculations not needed.\n');
  end;
else
  
  % OK, so first we need to calculate PCA
  % Check to see if we already have the PCA data
  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations supplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    % display notice if the user entered one, but not both, of E and D.
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;
    
    % Calculate PCA
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

% skip the rest if user only wanted PCA
if only > 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Whitening the data
  
  % Check to see if the whitening is needed...
  if jumpWhitening == 3,
    if b_verbose,
      fprintf ('Whitening not needed.\n');
    end;
  else
    
    % Whitening is needed
    % display notice if the user entered some of the whitening info, but not all.
    if (jumpWhitening > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump whitening:\n');
      fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
    end;
    
    % Calculate the whitening
    [Z, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
  end
end % if only > 1
[P,BT]=size(Z);
B = zeros(P);
for p = 1:P
    w = randn(P,1);
    %w =[0.1121;-1.4965;0.2293];
    w = w - B * B' * w;
    w = w / norm(w);
    wold = zeros(P,1);
    counter=0;
    for counter = 1 : maxNumIterations + 1,
        if counter == maxNumIterations + 1,
            fprintf('No convergence after %d steps\n', counter);
        break;
        end
    %wold = w; 
    hypTan=tanh(Z'*w);	
    w = (Z * hypTan - 1 * sum(1 - hypTan .^ 2)' * w) / BT;
   
    % Decorrelation:
    w = w - B*B'*w;
    w = w / norm(w);
   
    %test d'arret
        if norm(w - wold) < epsilon | norm(w + wold) < epsilon
            fprintf('Convergence after %d steps\n', counter); 
        break;
        end  
    wold = w; 
    end
    B(:,p) = w;
    Aes(:,p) = dewhiteningMatrix * w;
    Wes(p,:) = w' * whiteningMatrix;
    J(1,p)=counter;
end 


