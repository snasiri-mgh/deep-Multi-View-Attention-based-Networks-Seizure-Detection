function [Out1, Out2, Out3,maxit] = fastica(mixedsig, varargin)
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

[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose           = 'on';

firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      jumpPCA = jumpPCA + 1;
      E = varargin{i+1};
     case 'pcad'
      jumpPCA = jumpPCA + 1;
      D = varargin{i+1};
     case 'whitesig'
      jumpWhitening = jumpWhitening + 1;
      whitesig = varargin{i+1};
     case 'whitemat'
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = varargin{i+1};
     case 'dewhitemat'
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = varargin{i+1};
     case 'only'
      switch lower (varargin{i+1})
       case 'pca'
	only = 1;
       case 'white'
	only = 2;
       case 'all'
	only = 3;
      end
      
     otherwise
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_verbose
  fprintf('Number of signals: %d\n', Dim);
  fprintf('Number of samples: %d\n', NumOfSampl);
end



if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if jumpWhitening == 3
  if b_verbose,
    fprintf ('Whitened signal and corresponding matrises supplied.\n');
    fprintf ('PCA calculations not needed.\n');
  end;
else
  

  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations supplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;
    
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

if only > 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if jumpWhitening == 3,
    if b_verbose,
      fprintf ('Whitening not needed.\n');
    end;
  else
    

    if (jumpWhitening > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump whitening:\n');
      fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
    end;
    
    [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
  end
  
end % if only > 1

if only > 2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
  Dim = size(whitesig, 1);
  
  if numOfIC > Dim
    numOfIC = Dim;
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating only %d independent components\n', numOfIC);
      fprintf('(Can''t estimate more independent components than dimension of data)\n');
    end
  end
  
  [A, W,maxit] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
		  numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
		  maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
		  displayMode, displayInterval, verbose);
  
  if ~isempty(W)
    if b_verbose
      fprintf('Adding the mean back to the data.\n');
    end
    icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
    %icasig = W * mixedsig;
    if b_verbose & ...
	  (max(abs(W * mixedmean)) > 1e-9) & ...
	  (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
      fprintf('Note that the plots don''t have the mean added.\n');
    end
  else
    icasig = [];
  end

end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if only == 1    % only PCA
  Out1 = E;
  Out2 = D;
elseif only == 2  % only PCA & whitening
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else      % ICA
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = icasig;
    Out2 = A;
    Out3 = W;
  end
end
