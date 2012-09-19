function [A, varargout]= nonnegative_unmix( S, X, varargin )
% nonnegative_unmix
% is the non_negative solution in the least squares sense to the under- or
% overdetermined system of equations SA = X. further details can be found
% at MATLAB function reference of '\'
%
% usage:
%   A = nonnegative_unmix( S, X ...,'PropertyName',PropertyValue,...)
%
% with input and output arguments (Properties are optional):
% X                 (matrix)  n dimensional data matrix with the spectral
%                             information in the last dimension.
% S                 (matrix)  two dimensional reference matrix with the
%                             shape spectra x no_of_spectra
%
%                             The size of the last dimension of X and S
%                             must be the same.
%
% A                 (matrix)  output matrix.
%
% 'PropertyName',PropertyValue,
%                             sets properties to the specified property
%                             values mentioned below.
% Here is the set of implemented Properties and its allowed PropertyValues.
%
%  'verbose'        (boolean) Set to true if text output is to be
%                             generated. Default is false.
%  'GUI'            (boolean) Set to true if waitbar is to be used.
%                             Default is false.
%  'single'         (boolean) Set to true to use single precision instead
%                             of double (due to memory limits).
%                             Default is false.
%  'image'          (boolean) Set to true if the input matrix X is in shape
%                             of an image. Default is true.
%  'timer'          (boolean) Starts timer 'TIC'. Default is true.
%
% Examples:
% A = nonnegative_unmix( S, X )
% starts nonnegative_unmix with default values and returns the 'Concentration' matrix
%
% A = nonnegative_unmix( S, X, 'GUI',1, 'verbose',1, 'single',0, 'image',1, 'timer', 1);
% starts poissonNMF, where all 'Properties' are set
%
%
% 12.06.2007 (Version 1.0)
% Authors:
%   Andre Zeug (azeug@gwdg.de)
%
% This software is for non-commercial use only.

par.nonnegative_unmix_date='12.06.2007';
par.nonnegative_unmix_version='1.0';

%% Start - check and set parameters

if nargin<2 || ~isnumeric(S) || ~isnumeric(X),
  help nonnegative_unmix;
  A=[];
  return;
end

% parse input params
par.gui         = false;
par.verbose     = false;
par.singleprec  = false;
par.image       = true;
par.timer       = true;

for i=1:2:length(varargin),
  % Check that all argument types are strings
  if ~ischar(varargin{i}), error('Invalid input identifier names or input argument order.'); end
  % Lower/uppercase in identifier types doesn't matter:
  identifier=lower(varargin{i});     % identifier (lowercase)
  value=varargin{i+1};
  switch identifier
    case 'verbose'
      par.verbose = value;
    case 'gui'
      par.gui = value;
    case 'single'
      par.singleprec = value;
    case 'image'
      par.image = value;
    case 'timer'
      par.timer = value;
    otherwise
      %%% Unknown identifier
      error(['Invalid argument identifier ''' identifier '''!']);
  end
end
%% setup
% start timer for verbose and par.gui, TIC is global!!!
if par.timer, tic; end;              
t1=toc; 

% outprint
if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: starting calculation\n',t ); t1=t; end

% reshaping of X in case of image 
if par.image
  sz = size(X);
  X = reshape(X, [prod(sz(1:end-1)) sz(end)]);
end
X=X';
sz_X = size(X); % size of X
sz_S = size(S); % size of S

if sz_S(1)~=sz_X(1)
  S=S';
end

sz_S = size(S); % size of S

% outprint to check dimensions
if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: Matrix: %d x %d (%s), Reference: %d x %d (%s)  \n',...
    t, sz_X(1), sz_X(2), class(X), sz_S(1), sz_S(2), class(S) ); t1=t; end

% input matrix control
if sz_S(1)~=sz_X(1), error(['nonnegative_unmix: Matrix dimensions do not match!']); end;

% calculates the solution in the least squares sense to the under- or
% overdetermined system of equations SA = X.
% due to memory limitations the direct way might cause an error
%try
if length(X(:))<80000000
  % outprint to control time for calculation
  if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: Trying calculate estimate for A (%0.5f)\n',t ,t-t1); t1=t; end

  % distinguish single and double precession
  if par.singleprec
    if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: using SINGLE precision \n',t); t1=t; end
    if nargout > 1
      [A AError] = myBackSlash(single(S) ,  single(X));
    else
      A = single(S) \ single(X);
    end
  else
    if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: using DOUBLE precision \n',t); t1=t; end
    if nargout > 1
      [A AError] = myBackSlash(  double(S) , double(X));
    else
      A = double(S) \ double(X);
    end
  end
%catch
else
  % catch routime for error handling
  if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: Try failed, start subsequent method (%0.5f)\n',t ,t-t1); t1=t; end
  if par.gui, h = waitbar(0,'Please wait again...','Name','nonnegative_unmix: Wait even more!');end;

  % subsequential calculation
  if par.singleprec
    A = zeros(sz_S(2),sz_X(2),'single');
  else
    A = zeros(sz_S(2),sz_X(2),'double');
  end
  if nargout > 1
    if par.singleprec
      AError = zeros(sz_S(2),sz_X(2),'single');
    else
      AError = zeros(sz_S(2),sz_X(2),'double');
    end
  end
  for ii=0:1000000:sz_X(2) % subsequentially...
    if par.gui, waitbar(ii/sz_X(2),h,['Calculated estimate for A: ' num2str(100*ii/sz_X(2),'%0.2f') '%']);end;
    ii_end=min(sz_X(2),ii+1000000)-ii;
    % distinguish single and double precession
    if par.singleprec
      if nargout > 1
        [A(:,(1:ii_end)+ii) AError(:,(1:ii_end)+ii)] = myBackSlash( single(S) ,  single(X(:,(1:ii_end)+ii)));
      else
        A(:,(1:ii_end)+ii) = single(S) \ single(X(:,(1:ii_end)+ii));
      end
    else
      if nargout > 1
        [A(:,(1:ii_end)+ii) AError(:,(1:ii_end)+ii)] = myBackSlash(  double(S) , double(X(:,(1:ii_end)+ii)));
      else
        A(:,(1:ii_end)+ii) = double(S) \ double(X(:,(1:ii_end)+ii));
      end
    end
  end

  if par.gui, delete(h); drawnow; end;

end
%% ERROR handling 
% in certain cases the '\' operator returns NaN values, which can not be
% reproduced by calculating a second time.
% Might be removed in later versions.
nans=sum(isnan(A(:)));
if nans>0
  t=toc; fprintf('(%07.3fs) nonnegative_unmix:  FOUND %d NaNs ! ! !  (%0.5f)\n', t, nans, t-t1); t1=t;
  [~, yy]=find(isnan(A)==1);
  if par.singleprec,
    A(:,yy)=single(S)\single(X(:,yy));
  else
    A(:,yy)=double(S)\double(X(:,yy));
  end

end
nans=sum(isnan(A(:)));
if nans>0
  t=toc; fprintf('(%07.3fs) nonnegative_unmix:  STILL FOUND %d NaNs ! ! !  (%0.5f)\n', t, d, t-t1); t1=t;
end

%%
% outprint end of calculation
if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: Calculation of estimate for A finished (%0.5f)\n',t ,t-t1); t1=t; end

% spectra index
ind_S=1:sz_S(2);

% search for negative values and calculates estimate without negative
% component (recursively)
for i=1:sz_S(2)
  ind = find(A(i,:)<=0);
  if ~isempty(ind)
    A(i,ind) =0;
    subind_S=ind_S(ind_S ~= i );
    % starts nonnegative_unmix with the subset (without negative component)
    % with the conditions 'image' = false (not required) and 
    % 'timer' = false to avoid restart of the TIC function
    if nargout > 1
      AError(i,ind) =0;
      [A(subind_S,ind) AError(subind_S,ind)]=nonnegative_unmix(S(:,subind_S)', X(:,ind)',...
        'GUI',par.gui, 'verbose',par.verbose, 'single',par.singleprec, 'image',0, 'timer', 0);
    else
      A(subind_S,ind)=nonnegative_unmix(S(:,subind_S)', X(:,ind)',...
        'GUI',par.gui, 'verbose',par.verbose, 'single',par.singleprec, 'image',0, 'timer', 0);
    end

    
  end
end


% reshaping of A according to X 
if par.image
  A = reshape(A',[sz(1:end-1) sz_S(2)]);
  if nargout>1
    AError = reshape(AError',[sz(1:end-1) sz_S(2)]);
  end
end
%% transfer output
if nargout>1, 
  varargout(1) = {AError}; 
end


% final outprint
if par.verbose, t=toc; fprintf('(%07.3fs) nonnegative_unmix: Calculation of A finished (%0.5f)\n',t ,t-t1); t1=t; end

end

function [dData dError] = myBackSlash(A, B) 
  dData = A \ B; 
  y_predict = A * dData; 
  s = sqrt( sum((B - y_predict).^2)  / (numel(A(:)) -numel(dData(:)))); 
  XTX_inverse = (A'*A)^-1; 
  dError = s * sqrt( diag( XTX_inverse , 0) ); 
end

