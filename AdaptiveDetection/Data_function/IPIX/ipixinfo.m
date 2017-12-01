function [nc nrange nsweep ntxpol nadc cdfFileName]= ncdump(theNetCDFFile, theOutputFile)

% Derived from ncdump.
% Lists variables, units, descriptions, and scalar values of netcdf file.

if nargin < 1, help ncdump, theNetCDFFile = '*.cdf'; end
if nargin < 2, theOutputFile = 'stdout'; end   % stdout.

if isa(theNetCDFFile, 'ncitem')
    theNetCDFFile = name(parent(parent(theNetCDFFile)));
end

if any(theNetCDFFile == '*')
   theFilterSpec = '*.cdf';%theNetCDFFile;
   thePrompt = 'Select a NetCDF Input File:';
   [theFile, thePath] = uigetfile(theFilterSpec, thePrompt);
   if ~any(theFile), return, end
   theNetCDFFile = [thePath theFile];
end

if any(theOutputFile == '*')
%    theFilterSpec = theOutputFile;
%    thePrompt = 'Select a Text Output File:';
%    [theFile, thePath] = uiputfile(theFilterSpec, thePrompt);
%    if ~any(theFile), return, end
   theOutputFile = [theNetCDFFile '.txt'];
end

cdfFileName = theNetCDFFile;

nctypes = ['byte   '; 'char   '; 'short  '; ...
           'long   '; 'float  '; 'double '; ...;
           'unknown'; 'unknown'; 'unknown'];

nc = netcdf.open ( theNetCDFFile, 'NC_NOWRITE' );
%theNCid = ncid(nc);

if isempty(nc)
   disp([' ## Unable to open: ' theNetCDFFile])
   return
end

if strcmp(theOutputFile, 'stdout')
   fp = 1;
  elseif strcmp(theOutputFile, 'stderr')
   fp = 2;
  elseif isstr(theOutputFile)
   fp = fopen(theOutputFile, 'w');
  else
   fp = theOutputFile;
end

if fp < 0, close(nc), return, end

[ndims, nvars, ngatts, recdim] = netcdf.inq(nc);

[dims,ndims] = netcdf.inqDim(nc,0);

% Get name of global attribute
gattname = netcdf.inqAttName(nc,netcdf.getConstant('NC_GLOBAL'),0);
[gatts,ngattLen] = netcdf.inqAtt(nc,netcdf.getConstant('NC_GLOBAL'),gattname);

s = ['File: ' theNetCDFFile];
fprintf(fp, '%s\n', s);

s = '%% Variables:';
fprintf(fp, '\n%s\n', s);
s = '%% (none)';
if nvars < 1, fprintf(fp, '%s\n', s), end
for j = 1:nvars;
  [varname, theDatatype, dimids, natts] = netcdf.inqVar(nc,j-1);

  varname = strrep(varname, '''', '''''');
  ndims = length(dimids);
  
  theDatatype = ['nc' theDatatype];
  
  varid = netcdf.inqVarID(nc,varname);
 
  varData = netcdf.getVar(nc,varid);
  sz=size(varData);
  elements = prod(sz);   
  if elements==1
      theValue=num2str(varData);%(nc{j}(1)); 
  else 
    % unwrap angles
    if natts && strcmp(varname,'deg'),
      val=nc{j}(:);
      val=unwrap(val); if val(1)<0, val=val+360; end
    end
    % compact display of arrays and matrices 
    if prod(sz)==max(sz) & elements>=3,
      val=varData;
      [dummy,i]=max(sz);
      dimname = varname;
      theValue=['[' num2str(val(1)) ' ' num2str(val(2)) ' ..' dimname '.. ' num2str(val(end)) ']'];
    else 
      for i = 1:ndims
        dimname = varname;
        if i==1, s=['[' dimname];
        else s = [s ' x ' dimname]; end
      end
      s=[s ']'];
      theValue=s; 
    end
  end
  if natts
      attname=netcdf.inqAttName(nc,varid,0);%unit
      attUnit=netcdf.getAtt(nc,varid,attname);
      theUnit=attUnit; 
  else
      theUnit=''; 
  end
  fprintf(fp, '%24s = %10s %s\n',varname,theValue,strrep(theUnit,'\0',''));
end

s = '%% Dimensions:';
fprintf(fp, '\n%s\n', s);
if ndims < 1, disp('%% (none)'), end
for i = 1:ndims
    dimid = i-1;
    [dimname, dimlen] = netcdf.inqDim(nc,dimid);
    
    if strcmp(dimname,'nrange') nrange=dimlen, end
    if strcmp(dimname,'nsweep') nsweep=dimlen, end
    if strcmp(dimname,'ntxpol') ntxpol=dimlen, end
    if strcmp(dimname,'nadc')   nadc=dimlen, end
    if ndims==3 nadc=0, end
    %dimname = name(dims{i});
    dimname = strrep(dimname, '''', '''''');
    
    %dimlen = ncsize(dims{i});
    fprintf(fp, '%24s = %10s\n',dimname,int2str(dimlen));
end

s = '%% Global attributes:';
fprintf(fp, '\n%s\n', s);
s = '%% (none)';
if ngatts < 1,fprintf(fp, '%s\n', s); end
for i = 1:ngatts
  varid = -1;
  attnum = i-1;
%   attname =netcdf.inqattname(nc,netcdf.getConstant('NC_GLOBAL'),attnum);
    attname =netcdf.inqAttName(nc,netcdf.getConstant('NC_GLOBAL'),attnum);
  if any(attname ~= '_')
    while attname(1) == '_'
      attname = [attname(2:length(attname)) attname(1)];
    end
  end
  attname = strrep(attname, '''', '''''');
%   theDatatype = datatype(gatts{i});
%   attlen = ncsize(gatts{i});
%   attvalue = gatts{i}(:);
% Get information about the attribute.
  [theDatatype,attlen] = netcdf.inqAtt(nc,varid,attname);
  
  % Get value of attribute.
  attvalue = netcdf.getAtt(nc,varid,attname);  
  theDatatype = ['nc' theDatatype];
  s = attname;
  t = mat2str(attvalue);
  if length(t) > 0 & 0
    if t(1) == '[' & t(length(t)) == ']'
      t = [ '{' t(2:length(t)-1) '}'];
    end
  end
  if ~isstr(attvalue)
    if (0)
      f = [];
      k = 1:length(t)-1;
      if any(k), f = find(t(k) == t(k+1)); end
      if any(f), t(f) = []; end
       f = find(t == ' ');
      if any(f), t(f) = setstr(t(f) .* 0 + ','); end
      t = strrep(t, ',', ', ');
    end
  end
  fprintf(fp, '%24s = %10s\n',attname,strrep(attvalue,'\0',''));
end

if ischar(theOutputFile) & fp > 2, fclose(fp); end

% netcdf.Close(nc)

% if nargout > 0, theStatus = status; end

