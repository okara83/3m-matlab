function P=findpeaksb(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype,window,PeakShape,extra,AUTOZERO)

% Input arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity. Higher values will neglect smaller features.
% "smoothwidth" is the width of the smooth applied before peak detection;
% larger values ignore narrow peaks. If smoothwidth=0, no smoothing is
% performed.
% "seakgroup" is the number points around the top part of the peak that are
% taken for thge initial estimation. If Peakgroup=0 the local maximum is
% takes as the peak height and position.
% "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) 
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% "window" specifies the number of data points over which the peak is fit 
% to the model shape.   
% "PeakShape" specifies the model peak shape: "peakshape" = 1-20.
% (1=Gaussian, 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=sigmoid, 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian,
% 16=Fixed-position Gaussians; 17=Fixed-position Lorentzians;
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile.
%
% AUROZERO=0 does not subtract baseline from data segment. 
% AUROZERO=1 (default) subtracts interpolated linear baseline from data segment.
% AUROZERO=2, subtracts interpolated quadratic baseline from data segment.
% AUROZERO=3, corrects flat baseline without reference to the data segment.

%
% Example 1: Three Gaussians on a linear background
% x=1:.2:100;
% y=2-(x./50)+modelpeaks(x,3,1,[1 2 3],[20 50 80],[3 3 3])+.02*randn(size(x));
% disp('          Peak      Position      Height      Width        Area         % error')
% plot(x,y);NoBackgroundSubtraction=findpeaksb(x,y,.00005,.5,30,20,3,150,1,0,0)
% LinearBackgroundSubtraction=findpeaksb(x,y,.00005,.5,30,20,3,150,1,0,1)
%
% Example 2: Three Lorentzians on a curved background, quadratic background subtraction
% x=1:.2:100;
% y=modelpeaks(x,4,2,[10 2 2 3],[-10 20 50 80],[50 3 3 3])+.02*randn(size(x));
% disp('          Peak      Position      Height      Width        Area         % error')
% plot(x,y);NoBackgroundSubtraction=findpeaksb(x,y,.00005,.5,15,15,3,150,2,0,0)
% plot(x,y);QuadraticBackgroundSubtraction=findpeaksb(x,y,.00005,.5,15,15,3,150,2,0,2)

NumPeaks=1; % Number of models shapes to youe to fit each peak (usually 1)
NumTrials=1;
start=0;
fixedparameters=0; % Used for fixed parameters shapes only (11,16,12,and 17)
plots=0; % Supress plotting
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=y;
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>3,
                    startpoint=round(j-window/2);
                    if startpoint<1;startpoint=1;end
                    endpoint=round(j+window/2);
                    if endpoint>length(x);endpoint=length(x); end
                    XXX=x((startpoint):(endpoint));
                    YYY=y((startpoint):(endpoint));  
                    signal=[XXX',YYY'];
                    [FitResults,Error]=peakfit(signal,x(j),window,NumPeaks,PeakShape,extra,NumTrials,start,AUTOZERO,fixedparameters,plots);                
                    PeakX=FitResults(2);
                    PeakY=FitResults(3);
                    MeasuredWidth=FitResults(4);
                else
                    % if the peak is too narrow for least-squares technique
                    % to work well, just use the max value of y in the
                    % sub-group of points near peak.
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurements fails and
                % results in NaN, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip any peak that gives a NaN or is less that the
                    % AmpThreshold
                else % Otherwiase count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth Error];
                    peak=peak+1; % Move on to next peak
                end % if isnan(PeakX)...
            end % if y(j) > AmpTest...
        end  % if d(j)-d(j+1) > SlopeThreshold...
    end  % if sign(d(j)) > sign (d(j+1))...
end  % for j=....
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)

% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function d=deriv(a)
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. 
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end
% ----------------------------------------------------------------------
function [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots) 
global AA xxx PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO

format short g
format compact
warning off all

NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  Must be isignal(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % Must be isignal(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';

% Define values of any missing arguments
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 4
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 5
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 8
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 9
        AUTOZERO=autozero;
        FIXEDWIDTH=0;
        plots=1;
    case 10
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
        plots=1;
    case 11
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
    otherwise
end % switch nargin

% Default values for placeholder zeros
if NumTrials==0;NumTrials=1;end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
if FIXEDWIDTH==0, FIXEDWIDTH=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3,end
if AUTOZERO<0,AUTOZERO=0,end

% Remove linear baseline from data segment if AUTOZERO==1
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

% Assign ShapStrings
switch peakshape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Sigmoid';
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='BiLorentzian';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt profile';
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
%         plot(zxx,zyy);pause
        TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart);     
      case 19
          TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart);
      case 20
          % extra=extra+.2*randn().*extra;
          TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        otherwise
    end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)
%
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=sigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDWIDTH);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDWIDTH);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
       otherwise
  end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
heightsize=size(height');
AAsize=size(AA);
if heightsize(2)==AAsize(1),
   mmodel=height'*AA+baseline;
else
    mmodel=height*AA+baseline;
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
for m=1:NumPeaks,
    if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
    yi(m,:)=height(m)*AA(m,:); % (NEW) Place y values of individual model peaks into matrix yi
end
xi=xxx+xoffset; % (NEW) Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed lines
    if peakshape==16||peakshape==17
    else
        for marker=1:NumPeaks,
            markx=BestStart((2*marker)-1);
            subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
        end % for
    end % if peakshape
    plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
    hold off;
    axis([min(xx)+xoffset max(xx)+xoffset min(yy) max(yy)]);
    switch AUTOZERO,
        case 0
            title('peakfit 4.2   No baseline correction')
        case 1
            title('peakfit 4.2   Linear baseline subtraction')
        case 2
            title('peakfit 4.2   Quadratic subtraction baseline')
        case 3
            title('peakfit 4.2   Flat baseline correction')
    end
 
    switch peakshape
    case {4,20}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Shape Constant = ' num2str(extra) ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Time Constant = ' num2str(extra) ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     % Gaussian = ' num2str(extra) ] )
    case {14,15}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Asymmetry = ' num2str(extra) ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '% ' ] )
    end

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'b.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
end % if plots

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                else
                    FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
                end
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                else
                    FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% Display Fit Results on upper graph
if plots,
    subplot(2,1,1);
    startx=min(xx)+xoffset+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=(max(yy)-min(yy))./10;
    starty=max(yy)-dyy;
    FigureSize=get(gcf,'Position');
    if peakshape==9||peakshape==10||peakshape==19,  % Pulse and sigmoid shapes only
        text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
    else
        text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area'] );
    end
    % Display FitResults using sprintf
    for peaknumber=1:NumPeaks,
        for column=1:5,
            itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
    xposition=startx;
    yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
    if AUTOZERO==3,
       text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
    end
end % if plots

if NumArgOut==6,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(5,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDWIDTH);
        for peak=1:NumPeaks,
            BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:5);
        BootstrapSTD=BootstrapSTD(2:5);
        BootstrapIQR=BootstrapIQR(2:5);
        PercentRSD=PercentRSD(2:5);
        PercentIQR=PercentIQR(2:5);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==6,
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)


global PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO
format short g
format compact
warning off all
FIXEDWIDTH=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstar);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,zxx,zyy,-extra)),newstart);
      case 19
          TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart);
      case 20
        % extra=extra+.2*randn().*extra;
        TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch peakshape
    case 1
        A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 3
        A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 4
        A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 6
        A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
    case 7
        A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));  
    case 8
        A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';    
    case 9
        A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));  
    case 10
        A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 11
        A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
    case 12
        A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH); 
    case 13
        A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 14
        A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 15
        A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 16
        A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 17
        A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 18
        A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 19
        A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 20
        A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra); 
   end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)

for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12,  % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)

% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitlorentziancw(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix

% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
tlogistic assumes a logistic function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fitlognormal(lambda,t,y)
%	Fitting function for lognormal, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlognormal assumes a lognormal function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lognormal(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band signal.
% function P=findpeaksb(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype,window,PeakShape,extra,AUTOZERO)
% function  P=P=findpeaksb(x,y,SlopeThreshold,AmpThreshold,
% smoothwidth,peakgroup,smoothtype,window,PeakShape,extra, AUTOZERO)
% Function to locate and fit the positive peaks in a noisy x-y time series
% data set and fits each peak to a selected peak model with the peakfit.m
% function with selectable window width, peak shape, and optional
% background subtraction. Finds peaks by looking for downward
% zero-crossings in the first derivative that exceed SlopeThreshold.
% Returns the peak number, peak position, peak height, peak width, peak
% area, and the percent fitting error, for each peak.
%
% Input arguments "slopeThreshold", "ampThreshold" and "smoothwidth"
% control peak sensitivity. Higher values will neglect smaller features.
% "smoothwidth" is the width of the smooth applied before peak detection;
% larger values ignore narrow peaks. If smoothwidth=0, no smoothing is
% performed.
% "seakgroup" is the number points around the top part of the peak that are
% taken for thge initial estimation. If Peakgroup=0 the local maximum is
% takes as the peak height and position.
% "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) 
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% "window" specifies the number of data points over which the peak is fit 
% to the model shape.   
% "PeakShape" specifies the model peak shape: "peakshape" = 1-20.
% (1=Gaussian, 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=sigmoid, 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian,
% 16=Fixed-position Gaussians; 17=Fixed-position Lorentzians;
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile.
%
% AUROZERO=0 does not subtract baseline from data segment. 
% AUROZERO=1 (default) subtracts interpolated linear baseline from data segment.
% AUROZERO=2, subtracts interpolated quadratic baseline from data segment.
% AUROZERO=3, corrects flat baseline without reference to the data segment.
%
% Note: this funciton works best for isolated peaks that don't overlap; for
% overlapping peaks, findpeaksfit.m works better but is best suited to
% small numbers of peaks, whereas findpeaksb.m works with any number of
% peaks.
%
% Example 1: Three Gaussians on a linear background
% x=1:.2:100;
% y=2-(x./50)+modelpeaks(x,3,1,[1 2 3],[20 50 80],[3 3 3])+.02*randn(size(x));
% disp('          Peak      Position      Height      Width        Area         % error')
% plot(x,y);NoBackgroundSubtraction=findpeaksb(x,y,.00005,.5,30,20,3,150,1,0,0)
% LinearBackgroundSubtraction=findpeaksb(x,y,.00005,.5,30,20,3,150,1,0,1)
%
% Example 2: Three Lorentzians on a curved background, quadratic background subtraction
% x=1:.2:100;
% y=modelpeaks(x,4,2,[10 2 2 3],[-10 20 50 80],[50 3 3 3])+.02*randn(size(x));
% disp('          Peak      Position      Height      Width        Area         % error')
% plot(x,y);NoBackgroundSubtraction=findpeaksb(x,y,.00005,.5,15,15,3,150,2,0,0)
% plot(x,y);QuadraticBackgroundSubtraction=findpeaksb(x,y,.00005,.5,15,15,3,150,2,0,2)



NumPeaks=1; % Number of models shapes to youe to fit each peak (usually 1)
NumTrials=1;
start=0;
fixedparameters=0; % Used for fixed parameters shapes only (11,16,12,and 17)
plots=0; % Supress plotting
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=y;
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>3,
                    startpoint=round(j-window/2);
                    if startpoint<1;startpoint=1;end
                    endpoint=round(j+window/2);
                    if endpoint>length(x);endpoint=length(x); end
                    XXX=x((startpoint):(endpoint));
                    YYY=y((startpoint):(endpoint));  
                    signal=[XXX',YYY'];
                    [FitResults,Error]=peakfit(signal,x(j),window,NumPeaks,PeakShape,extra,NumTrials,start,AUTOZERO,fixedparameters,plots);                
                    PeakX=FitResults(2);
                    PeakY=FitResults(3);
                    MeasuredWidth=FitResults(4);
                else
                    % if the peak is too narrow for least-squares technique
                    % to work well, just use the max value of y in the
                    % sub-group of points near peak.
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurements fails and
                % results in NaN, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip any peak that gives a NaN or is less that the
                    % AmpThreshold
                else % Otherwiase count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth Error];
                    peak=peak+1; % Move on to next peak
                end % if isnan(PeakX)...
            end % if y(j) > AmpTest...
        end  % if d(j)-d(j+1) > SlopeThreshold...
    end  % if sign(d(j)) > sign (d(j+1))...
end  % for j=....
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values

% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function d=deriv(a)
% First derivative of vector using 2-point central difference.
%   O, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. 
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end
% ----------------------------------------------------------------------
function [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots)

global AA xxx PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO

format short g
format compact
warning off all

NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  Must be isignal(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % Must be isignal(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';

% Define values of any missing arguments
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 4
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 5
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 8
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 9
        AUTOZERO=autozero;
        FIXEDWIDTH=0;
        plots=1;
    case 10
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
        plots=1;
    case 11
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
    otherwise
end % switch nargin

% Default values for placeholder zeros
if NumTrials==0;NumTrials=1;end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
if FIXEDWIDTH==0, FIXEDWIDTH=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3,end
if AUTOZERO<0,AUTOZERO=0,end

% Remove linear baseline from data segment if AUTOZERO==1
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

% Assign ShapStrings
switch peakshape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Sigmoid';
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='BiLorentzian';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt profile';
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
%         plot(zxx,zyy);pause
        TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart);     
      case 19
          TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart);
      case 20
          % extra=extra+.2*randn().*extra;
          TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        otherwise
    end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)
%
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=sigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDWIDTH);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDWIDTH);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
       otherwise
  end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
heightsize=size(height');
AAsize=size(AA);
if heightsize(2)==AAsize(1),
   mmodel=height'*AA+baseline;
else
    mmodel=height*AA+baseline;
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
for m=1:NumPeaks,
    if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
    yi(m,:)=height(m)*AA(m,:); % (NEW) Place y values of individual model peaks into matrix yi
end
xi=xxx+xoffset; % (NEW) Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed lines
    if peakshape==16||peakshape==17
    else
        for marker=1:NumPeaks,
            markx=BestStart((2*marker)-1);
            subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
        end % for
    end % if peakshape
    plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
    hold off;
    axis([min(xx)+xoffset max(xx)+xoffset min(yy) max(yy)]);
    switch AUTOZERO,
        case 0
            title('peakfit 4.2   No baseline correction')
        case 1
            title('peakfit 4.2   Linear baseline subtraction')
        case 2
            title('peakfit 4.2   Quadratic subtraction baseline')
        case 3
            title('peakfit 4.2   Flat baseline correction')
    end
 
    switch peakshape
    case {4,20}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Shape Constant = ' num2str(extra) ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Time Constant = ' num2str(extra) ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     % Gaussian = ' num2str(extra) ] )
    case {14,15}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Asymmetry = ' num2str(extra) ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '% ' ] )
    end

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'b.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
end % if plots

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                else
                    FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
                end
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                else
                    FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% Display Fit Results on upper graph
if plots,
    subplot(2,1,1);
    startx=min(xx)+xoffset+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=(max(yy)-min(yy))./10;
    starty=max(yy)-dyy;
    FigureSize=get(gcf,'Position');
    if peakshape==9||peakshape==10||peakshape==19,  % Pulse and sigmoid shapes only
        text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
    else
        text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area'] );
    end
    % Display FitResults using sprintf
    for peaknumber=1:NumPeaks,
        for column=1:5,
            itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
    xposition=startx;
    yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
    if AUTOZERO==3,
       text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
    end
end % if plots

if NumArgOut==6,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(5,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDWIDTH);
        for peak=1:NumPeaks,
            BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:5);
        BootstrapSTD=BootstrapSTD(2:5);
        BootstrapIQR=BootstrapIQR(2:5);
        PercentRSD=PercentRSD(2:5);
        PercentIQR=PercentIQR(2:5);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==6,
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)

global PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO
format short g
format compact
warning off all
FIXEDWIDTH=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstar);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,zxx,zyy,-extra)),newstart);
      case 19
          TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart);
      case 20
        % extra=extra+.2*randn().*extra;
        TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch peakshape
    case 1
        A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 3
        A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 4
        A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 6
        A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
    case 7
        A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));  
    case 8
        A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';    
    case 9
        A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));  
    case 10
        A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 11
        A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
    case 12
        A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH); 
    case 13
        A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 14
        A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 15
        A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 16
        A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 17
        A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 18
        A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 19
        A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 20
        A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra); 
   end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)

for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12,  % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% 
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitlorentziancw(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
%  O, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fitlognormal(lambda,t,y)
%	Fitting function for lognormal, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlognormal assumes a lognormal function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lognormal(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band signal.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar

g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
% figure(2);plot(ey);figure(1);
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponental pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)

g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% sigmiods of the form Height./(1 + exp((t1 - t)/t2))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=sigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); % Modified sigmoid
% Bifurcated sigmoid
% lx=length(x);
% hx=val2ind(x,t1);
% g(1:hx)=1./(1 + exp((t1 - x(1:hx))./t2));
% g(hx+1:lx)=1./(1 + exp((t1 - x(hx+1:lx))./(1.3*t2)));
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
g=2*((m/100)*gaussian(x,pos,wid)+(1-(m/100))*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for Voigt profile function
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
y = abs(xx-pos)/gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function err = fitBiLorentzian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiLorentzian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = BiLorentzian(x,pos,wid,m)
% BiLorentzian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.

lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=lorentzian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=lorentzian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a)  returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a.

mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar

g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar

g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
% figure(2);plot(ey);figure(1);
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponental pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
 
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% sigmiods of the form Height./(1 + exp((t1 - t)/t2))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=sigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); % Modified sigmoid
% Bifurcated sigmoid
% lx=length(x);
% hx=val2ind(x,t1);
% g(1:hx)=1./(1 + exp((t1 - x(1:hx))./t2));
% g(hx+1:lx)=1./(1 + exp((t1 - x(hx+1:lx))./(1.3*t2)));
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.

global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%
g=2*((m/100)*gaussian(x,pos,wid)+(1-(m/100))*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for Voigt profile function
% 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
y = abs(xx-pos)/gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
% 
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function err = fitBiLorentzian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% .
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiLorentzian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = BiLorentzian(x,pos,wid,m)
% BiLorentzian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=lorentzian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=lorentzian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a)  returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a.

mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));



% A self-contained demonstration of the OKPeak function. In this example, 
% the signal contains a repeated pattern of two overlapping Gaussian peaks
% of width 12, with a 2:1 height ratio. These patterns occur a random
% intervals, and the noise level is about 10% of the average peak height.
% Using OKPeak's ensemble average function (Shift-E), the patterns can be
% averaged and the signal-to-noise ratio significantly improved. 
format short g
format compact
increment=1;
x=[1:increment:18000];

% For each simulated peak, compute the amplitude, position, and width
pos=[200:50:17800];  % Positions of the peaks (Change if desired)
amp=round(10.*randn(1,length(pos)));  % Amplitudes of the peaks  (Change if desired)
wid=12.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=1; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplitude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>9,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2)+.5.*exp(-((x-pos(k)-wid(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
% y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
% demodata=[x' y']; % Assembles x and y vectors into data matrix
% 
% figure(1);plot(x,y,'r')  % Graph the signal in red
% title('Detected peaks are numbered. Peak table is printed in Command Window')

% Initial values of variable parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.5*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=3;
SmoothWidth=15;  % SmoothWidth should be roughly equal the peak width (in points)
FitWidth=15; % FitWidth should be roughly equal to the peak widths (in points)

% % Lavel the x-axis with the parameter values
% xlabel(['SlopeThresh. = ' num2str(SlopeThreshold) '    AmpThresh. = ' num2str(AmpThreshold) '    SmoothWidth = ' num2str(SmoothWidth) '    FitWidth = ' num2str(FitWidth) ])

% Find the peaks
tic;
Measuredpeaks=OKpeak([x;y],0,AmpThreshold,0.00015625,SmoothWidth,FitWidth,300,60,0);
ElapsedTime=toc;
% Zoom in on average peak
AveragePeakIndex=val2ind(Measuredpeaks(:,3),15);
AveragePeakPosition=Measuredpeaks(AveragePeakIndex,2);
Measuredpeaks=OKpeak([x;y],0,AmpThreshold,0.00015625,SmoothWidth,FitWidth,AveragePeakPosition,73,0);
PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;

% Display results
disp('---------------------------------------------------------')
disp('Demonstration of OKpeak''s ensemble average function (Shift-E)')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
disp(['Number of peaks detected = ' num2str(length(Measuredpeaks)) ] )
disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
disp(' ')
disp('Press Shift-E to compute the ensemble average of all peaks, then look at ' )
disp('the plot of the ensemble average in Figure 2 and compare its signal-to-noise ' )
disp('ratio to that of the original signal. The theoretical signal-to-noise improvement ')
disp(['is the square root of the number of peaks averaged = ' num2str(sqrt(length(Measuredpeaks))) ] )
disp(' ')
disp('To curve-fit the ensemble-averaged signal, type this: ')
disp('[parameters,error]=peakfit([1:length(EnsembleAverage);EnsembleAverage],40,60,2)')
disp(' ')
disp('To fit one of the peaks in the original signal, type this: ')
disp('[FitResults,MeanFitError]=peakfit([x;y],AveragePeakPosition,70,2,1,0,10)')
disp(' ')
disp('Compare those two fitting results. The actual repeated signal pattern')
disp('consists of two Gaussians, 12 points apart, width 12, with a 2:1 height ratio.')
disp('The results for the averaged signal pattern should be much closer.')


% Comparison of the accuracy of peak measurement between findpeaks and
% peakfit.
increment=1;
x=1:increment:500;

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400];   % Positions of the peaks (Change if desired)
wid=[40 30 20 10];   % Widths of the peaks (Change if desired)
Noise=0.5; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positio ns
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
plot(x,y)
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=WidthPoints; 
FitWidth=WidthPoints;

% call findpeaks and peakfit 
P=findpeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
% For Lorentzian peaks, use the folloiwng line instead
% P=findpeaksL(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
[FitResults,MeanFitError]=peakfit(demodata,0,0,4,1,0,0,0,0,0,0);

% Compute errors by subtracting ActualPeaks
NumPeaks= max(P(:,1));
if NumPeaks==4;
  FindpeaksError=100.*(P(1:4,2:5)-ActualPeaks(1:4,2:5))./ActualPeaks(1:4,2:5);
end
plot(x,y,'.r')
hold on
for peak=1:NumPeaks
    text(P(peak,2),P(peak,3),['Peak ' num2str(P(peak,1))])
end
hold off
PeakfitError=100.*(FitResults(1:4,2:5)-ActualPeaks(1:4,2:5))./ActualPeaks(1:4,2:5);
disp('Average percent errors for all peaks')
disp('      Position       Height      Width        Area')
PercentFindpeaksError=mean(abs(FindpeaksError))
PercentPeakfitError=mean(abs(PeakfitError))
RatioOfFindpeaksToPeakfitErrors=mean(PercentFindpeaksError./PercentPeakfitError)
disp('  ')

% To plot actual peak heights vs measured heights by findpeaks:
% plotfit(ActualPeaks(:,3),FitResults(:,3),1)
% To plot actual peak heights vs measured heights by peakfit:
% plotfit(ActualPeaks(:,3),P(:,3),1)

function OKpeakdemo


increment=1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 1 1 1];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400];   % Positions of the peaks (Change if desired)
wid=[10 30 50 70];   % Widths of the peaks (Change if desired)
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2);  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints/2); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=OKpeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('Detection and measurement of 4 peaks with the same heights but ')
disp('different widths. This demonstrates the effect of SlopeThreshold')
disp('and SmoothWidth on peak detection. Increasing SlopeThreshold (S key)')
disp('will discriminate against the broader peaks. Increasing SmoothWidth')
disp('(D key) will discriminate against the narrower peaks and noise.')
disp('')
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks

if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks;
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end

disp('End of demo.')

% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

function OKdemo1

increment=1;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=40.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=0.5; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.5*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2)+2;  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW,
% with AUTOZERO=1 (ON) and the top window showing the first peak. Once
% iPeak is running, you can change any of these using 
% keystroke commands (press "K" to see a list of commands).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),120,1);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp(['Signal to noise ratio of smallest peak = ' num2str(min(ActualPeaks(:,3))./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks(:,1:5)

SizeResults=size(MeasuredPeaks);
merror=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(1:5,2),MeasuredPeaks(n,2));
    merror(n,:)=100.*abs(ActualPeaks(indexn,:)-MeasuredPeaks(n,1:5))./ActualPeaks(indexn,:);
    merror(n,1)=MeasuredPeaks(n,1);
    % MeasuredPeaks(n,1)=indexn;
end
AveragePercentError=mean(merror)
disp('Demonstration of autozero background correction, for separated, narrow ')
disp('peaks on a large baseline. Hint: Turn on the Autozero mode (T key),')
disp('adjust the zoom setting so that the peaks are shown one at a time in')
disp('the upper window, then press the P key to display the peak table.')
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')

disp('End of demo.')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);



function OKdemo2

increment=1;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=60.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=0.05; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
  end; 
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2);  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints/2); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200,1);


% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp(['Signal to noise ratio of smallest peak = ' num2str(min(ActualPeaks(:,3))./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks(:,1:5)

SizeResults=size(MeasuredPeaks);
merror=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(1:5,2),MeasuredPeaks(n,2));
    merror(n,:)=100.*abs(ActualPeaks(indexn,:)-MeasuredPeaks(n,1:5))./ActualPeaks(indexn,:);
    merror(n,1)=MeasuredPeaks(n,1);
    % MeasuredPeaks(n,1)=indexn;
end
AveragePercentError=mean(merror)

disp('Demonstration of error caused by overlapping peaks on a large offset')
disp('baseline. Hint: Use the B key and click on the baseline points,')
disp('then press the P key to display the peak table. Or turn on the Autozero mode ')
disp('(T key) and use the Normal curve fit (N key) or Multiple curve fit (M key).') 
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')
disp('End of demo.')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);





function OKdemo3

increment=1;
x=[1:increment:1200];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4 5 6 7 8 9 10];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400 500 600 700 800 900 1000];   % Positions of the peaks (Change if desired)
wid=[10 15 20 25 30 35 40 45 50 55];   % Widths of the peaks (Change if desired)
Noise=.1; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      % A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2);  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints/2); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);


% Compare MeasuredPeaks to ActualPeakst
disp('-----------------------------------------------------------------')
disp(['Signal to noise ratio of smallest peak=' num2str(min(amp)./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks(:,1:5)

SizeResults=size(MeasuredPeaks);
merror=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(1:5,2),MeasuredPeaks(n,2));
    merror(n,:)=100.*abs(ActualPeaks(indexn,:)-MeasuredPeaks(n,1:5))./ActualPeaks(indexn,:);
    merror(n,1)=MeasuredPeaks(n,1);
    % MeasuredPeaks(n,1)=indexn;
end
AveragePercentError=mean(merror)
disp('Demonstration of overlapping Lorentzians peaks, without a background.')
disp('Overlap of peaks causes significant errors in peak height, ')
disp('width, and area. Hint: turn OFF the Autozero mode (T key) and use the Normal ')
disp('curve fit (N key) or Multiple curve fit (M key) with peak shape 2 (Lorentzian).') 
disp('End of demo.')

% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

function OKdemo4

increment=1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400];   % Positions of the peaks (Change if desired)
wid=[30 30 30 30];   % Widths of the peaks (Change if desired)
Noise=0.5; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positio ns
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=31;  % Larger-than-usual values of SmoothWidth and FitWidth 
FitWidth=40; % will help deal with the very poor signal-to-noise ratio.

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('Detection and measurement of 4 peaks in a very noisy signal')
disp(['Signal-to-noise ratio of first peak = ' num2str(amp(1)/Noise) ] )
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks

if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end
disp('The peak at x=100 is usually detected, but the accuracy of peak ')
disp('parameter measurement is poor because of the low signal-to-noise ratio.')
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')
disp(' ')
disp('An even better way to measure this signal is to use the peakfit')
disp('function, now that we know the number of peaks and their approximate')
disp('positions and widths. Zoom out until all 4 peaks are shown, press')
disp('the N key and press enter the peak shape (1) and number ot trials. ')
disp('Compare to the the iPeak results obtained by pressing the P key. ')
disp(' ')
disp('End of demo.')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);





function OKdemo5
increment=1;
Noise=2;
x=1:increment:501;
% For each simulated peak, specify the amplitude, position, and width
% and the overall random noise level.
amp=[100 200 300 400];  % Amplitudes of the peaks  (Change if desired)
pos=[175 225 275 325];   % Positions of the peaks (Change if desired)
wid=[50 50 50 50];   % Widths of the peaks (Change if desired)
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % alternatively, A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix
disp('-----------------------------------------------------------------')
disp('In this demo the peaks are so highly overlapped that only one or two of')
disp('the highest peaks yield distinct maxima that are detected by iPeak. The') 
disp('height, width, and area estimates are highly inaccurate because of this')
disp('overlap. The normal peak fit function (''N'' key) would be useful in this')
disp('case, but it depends on iPeak for the number of peaks and for the initial')
disp('guesses, and so it would fit only the peaks that were found and numbered.')  
disp(' ')
disp(['Signal to noise ratio of smallest peak = ' num2str(min(ActualPeaks(:,3))./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
iPeakResults=ipeak(demodata,0,100,8.2069e-007,15,15,248.5,353,0);
SizeResults=size(iPeakResults);
error=zeros(SizeResults(1),6);

disp(' ')
disp('To help in this case, pressing the ''H'' key will activate the peak ')
disp('sharpen function that decreases peak width and increases peak height')
disp('of all the peaks, making it easier for Findpeaks to detect and number ')
disp('them for use by the peakfit function. Note: peakfit fits the original')
disp('unmodified peaks; the sharpening is used only to help locate the peaks.')
disp(' ')
disp('Peak shape (1-8): 1 ')
disp('Number of trials: 1')
figure(2)
PeakfitResults=peakfit(demodata,251,352,4,1,1,1,[178 60 226 56 276 57 324 46],0); 
SizeResults=size(PeakfitResults);
error=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(:,2),PeakfitResults(n,2));
    error(n,:)=100.*abs(ActualPeaks(indexn,:)-PeakfitResults(n,:))./ActualPeaks(indexn,:);
    error(n,1)=indexn;
    PeakfitResults(n,1)=indexn;
end
if SizeResults(1)==1,
    AveragePercentErrors=error;
else
    AveragePercentErrors=mean(error);
end
PeakfitResults
AveragePercentErrors(1)=0;
Peakfit_Average_Percent_Errors=AveragePercentErrors
disp(' ')
disp('The plot of the peakfit results are shown in Figure 2')
disp('The built-in peakfit function (''N'') key gives much more accurate')
disp('values for the height, width, and area of these peaks')
disp(' ')
disp('End of demo. ')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);


function OKdemo6

% integers.
increment=1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 1];  % Amplitudes of the peaks  (Change if desired)
pos=[100 300];   % Positions of the peaks (Change if desired)
wid=[30 60];   % Widths of the peaks (Change if desired)
Noise=0.2; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positio ns
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*propnoise(z)+Noise.*pinknoise(length(z));  % Optionally adds proportional random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=33;  % Larger-than-usual values of SmoothWidth and FitWidth 
FitWidth=40; % will help deal with the very poor signal-to-noise ratio.

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('Detection and measurement of two peaks in a very noisy signal')
disp('You can change the peak heights, positions, and widths in lines 14-16,')
disp('the noise level in line 17, and the noise type in line 33.')
disp(['Signal-to-noise ratio of first peak = ' num2str(amp(1)/Noise) ] )
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks

if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end
disp('The accuracy of peak parameter measurement ')
disp('is poor because of the low signal-to-noise ratio.')
disp(' ')
disp('An even better way to measure this signal is to use the peakfit')
disp('function, now that we know the number of peaks and their')
disp('approximate, positions and widths. Press the M key and ')
disp('enter the peak shape (1) and number of trials. Compare to ')
disp('the iPeak results obtained by pressing the P key. ')
disp(' ')
disp('End of script.')

function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);

function y=propnoise(x)
% Random noise whose amplitude is proportional to the amplitude of signal
% zwith white power spectrum with mean zero and unit standard deviation,
% equal in length to x.
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*propnoise(model);
% plot(model)
% 
z=x.*randn(size(x));
sz=std(z);
y=z./sz;

function y=whitenoise(x)
% Random noise with white power spectrum with mean zero 
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*whitenoise(model);
% plot(model)
% 
y=randn(size(x));

function ry=pinknoise(n)

% model=gaussian([1:1000],500,200);
% model=model+.1.*pinknoise(length(model));
% plot(model)
x=[1:n];
y=randn(size(x));  % Random normally-distributed white noise
% Fourier filter 
fy=fft(y); % Compute Fourier transform of signal y
% Compute filter shape
lft1=[1:(length(fy)/2)+1];
lft2=[(length(fy)/2):length(fy)];
ffilter1=ones(size(lft1))./(sqrt(lft1));
ffilter2=ones(size(lft2))./(sqrt(lft2));
ffilter=[ffilter1,ffilter2];
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Inverse transform to recover filtered signal 'ry'
ry=((ry-mean(ry))./std(ry)); % Normalize to unit standard deviation

function y=sqrtnoise(x)
% Random noise whose amplitude is proportional to the square root of the
% amplitude of signal with white power spectrum with mean zero and unit
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*sqrtnoise(model);
% plot(model)
% 
z=sqrt(abs(x)).*randn(size(x));
sz=std(z);
y=z./sz;



function [P,FitResults,LowestError,BestStart,xi,yi]=findpeaksfit(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype,peakshape,extra,NumTrials,autozero,fixedparameters,plots)
% function [P,FitResults,LowestError,BestStart,xi,yi]=
% findpeaksfit(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,
% smoothtype,peakshape,extra,NumTrials,autozero,fixedparameters,plots)

% Function to locate and fit the positive peaks in a noisy x-y time series
% data set to a selected peak model. Detects peaks using 'findpeaks' (see 
.
%
% Returns the findpeaks peak table in P, and the peakfit results table in 
% FitResults,LowestError,BestStart,xi,yi
%
% SlopeThreshold - Slope of the smoothed first-derivative that is taken to
% indicate a peak. This discriminates on the basis of peak width. Larger
% values of this parameter will neglect broad features of the signal. A
% reasonable initial value for Gaussian peaks is 0.7*WidthPoints^-2, where
% WidthPoints is the number of data points in the half-width of the peak. 
%
% AmpThreshold - Discriminates on the basis of peak height. Any peaks with
% height less than this value are ignored. 
%
% SmoothWidth - Width of the smooth function that is applied to data before
% the slope is measured. Larger values of SmoothWidth will neglect small,
% sharp features. A reasonable value is typically about equal to 1/2 of the
% number of data points in the half-width of the peaks. 
%
% FitWidth - The number of points around the "top part" of the (unsmoothed)
% peak that are taken to estimate the peak heights, positions, and widths.
% A reasonable value is typically about equal to 1/2 of the number of data
% points in the half-width of the peaks. The minimum value is 3. 
%
% Smoothtype - determines the smoothing algorithm:
%     If smoothtype=1, rectangular (sliding-average or boxcar)
%     If smoothtype=2, triangular (2 passes of sliding-average)
%     If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
%    Basically, higher values yield  greater reduction in high-frequency
%    noise
%
% peakshape - Specifies the peak shape of the model: "peakshape" = 1-20.
% (1=Gaussian (default), 2=Lorentzian, 3=logistic, 4=Pearson,
% 5=exponentionally broadened Gaussian; 6=equal-width Gaussians;
% 7=Equal-width Lorentzians; 8=exponentionally broadened equal-width
% Gaussian, 9=exponential pulse, 10=sigmoid, 11=Fixed-width Gaussian,
% 12=Fixed-width Lorentzian; 13=Gaussian/Lorentzian blend; 14=BiGaussian,
% 15=BiLorentzian, 16=Fixed-position Gaussians; 17=Fixed-position
% Lorentzians; 18=exponentionally broadened Lorentzian; 19=alpha function;
% 20=Voigt profile.
%
% extra - Specifies the value of 'extra', used to fine-tune the peak shape
% fot the Pearson, exponentionally broadened Gaussian, Gaussina/:orentzian
% blend, and the bifurcated Gaussian and Lorentzian shapes. 
%
% NumTrials - number of trial fits and selects the best one (with lowest
% fitting error). NumTrials can be any positive integer (default is 1).
%
% autozero - sets autozero mode 0,1,2, or 3. autozero=0 (default) does not
% subtract baseline from data segment; autozero=1 interpolates a linear
% baseline from the edges of the data segment and subtracts it from the
% signal (assumes that the peak returns to the baseline at the edges of the
% signal); autozero=2,  like mode 1 except that it computes a quadratic
% curved baseline; autozero=3 compensates for a flat baseline without
% reference to the signal itself (best if the peak does not return to the
% baseline at the edges of the signal).
%
% fixedparameters - values (in vector format if more than one peak) of the
% peak positions or widths only for the fixed positions or width shapes
% (shapes 11, 12, 16, or 17). 
%
% plots - set plots=1 to plot the fitting results, othersize no plotting.
%
% Only the first 6 input arguments are required; the others have default
% values if omitted: smoothtype=3; peakshape=1; extra=0; NumTrials=1;
% autozero=0; fixedparameters=1; plots=0;
%
% Optional output parameters 
% 1. P: Peak table from findpeaks function (peak number, position, height,
%    width, area.
% 2. FitResults: a table of model peak parameters, one row for each peak,
%    listing Peak number, Peak position, Height, Width, and Peak area.
% 3. LowestError: The rms fitting error of the best trial fit.
% 4. BestStart: the starting guesses that gave the best fit.
% 5. xi: vector containing 100 interploated x-values for the model peaks. 
% 6. yi: matrix containing the y values of each model peak at each xi. 
%    Type plot(xi,yi(1,:)) to plot peak 1 or plot(xi,yi) to plot all peaks
% For more details, see


%
% Example 1: One noisy Gaussian peak at x=5, peak height=1.0, width=1.66,
% area=1.77. Peakfit gives better results, especially for width and area.
% x=[0:.05:10];y=exp(-(x-5).^2)+.3*randn(1,length(x));
% [findpeaksr,peakfitr]=findpeaksfit(x,y,0.0007,0.5,31,40,3,1,0,1,0,0,1)
%
% Example 2: Two Gaussian peaks at x=3 and 5, peak heights=0.5 and 1.0.
% x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.05*randn(1,length(x));
% [findpeaksr,peakfitr,error]=findpeaksfit(x,y,0.002,0.2,4,7,3,1,0,1,0,0,1)
%
% Example 3: Three Lorentzian peaks at x=30,40,50, all peak heights=1.0.
% x=1:.2:100;
% y=lorentzian(x,40,10)+lorentzian(x,50,10)+lorentzian(x,60,10)+.01.*randn(size(x));
% [findpeaksr,peakfitr]=findpeaksfit(x,y,8e-005,0.5,17,31,3,2,0,1,0,0,1)



format short g
format compact
warning off all

% Define values of any missing arguments for the peakfit function
switch nargin
     case 6
        smoothtype=3;
        peakshape=1;
        extra=0;
        NumTrials=1;
        autozero=0;
        fixedparameters=1;
        plots=0;
    case 7
        peakshape=1;
        extra=0;
        NumTrials=1;
        autozero=0;
        fixedparameters=1;
        plots=0;
    case 8
        extra=0;
        NumTrials=1;
        autozero=0;
        fixedparameters=1;
        plots=0;
    case 9
        NumTrials=1;
        autozero=0;
        fixedparameters=1;
        plots=1;
    case 10
        autozero=0;
        fixedparameters=1;
        plots=0;
    case 11
        fixedparameters=1;
        plots=0;
    case 12
        plots=0;
    otherwise
end % switch nargin

if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=y;
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>3,
                    [coef,S,MU]=polyfit(xx,log(abs(yy)),2);  % Fit parabola to log10 of sub-group with centering and scaling
                    c1=coef(3);c2=coef(2);c3=coef(1);
                    PeakX=-((MU(2).*c2/(2*c3))-MU(1));   % Compute peak position and height of fitted parabola
                    PeakY=exp(c1-c3*(c2/(2*c3))^2);
                    MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*c3)));
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak
                % detected, containing the peak number, peak
                % position (x-value) and peak height (y-value).
                % If peak measurements fails and results in NaN, skip this
                % peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip this peak
                else % Otherwiase count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end % isnan(PeakX)
            end  % y(j) > AmpTest,
        end  % d(j)-d(j+1) > SlopeThreshold
    end  % sign(d(j)) > sign (d(j+1)),
end % j=
sizeP=size(P);
NumPeaks=sizeP(1);
startvector=[];
for peaks=1:NumPeaks,
    startvector=[startvector [P(peaks,2) P(peaks,4)]];
end
if P(1)>0,
    [FitResults,LowestError,BestStart,xi,yi]=peakfit([x;y],0,0,NumPeaks,peakshape,extra,NumTrials,startvector,autozero,fixedparameters,plots);
else
    P=0;FitResults=0;LowestError=0;BestStart=0;xi=0;yi=0;
    disp('No peaks detected.')
end
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values

% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function d=deriv(a)
% First derivative of vector using 2-point central difference.

n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. 
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]

if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end
% ----------------------------------------------------------------------
function [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots)
% 
% Matlab function in a single m-file. Uses a non-linear optimization
% algorithm to decompose a complex, overlapping-peak signal into its
% component parts. The objective is to determine wpeakhether your signal
% can be represented as the sum of fundamental underlying peaks shapes.
% Accepts signals of any length, including those with non-integer and
% non-uniform x-values. Fits any number of peaks of Gaussian, Lorentzian,
% equal-width Gaussian and Lorentzian, fixed-width Gaussian and Lorentzian,
% biburfated Gaussian and Lorentzian, exponentially-broadened Gaussian,
% Pearson, Logistic, exponential pulse, sigmoid, Gaussian/Lorentzian blend,
% and Voigt shapes. This is a command line version, usable from a remote
% terminal. 
% peakfit(signal);       
% Performs an iterative least-squares fit of a single Gaussian  
% peak to the data matrix "signal", which has x values 
% in column 1 and Y values in column 2 (e.g. [x y])
%
% peakfit(signal,center,window);
% Fits a single Gaussian peak to a portion of the 
% matrix "signal". The portion is centered on the 
% x-value "center" and has width "window" (in x units).
% 
% peakfit(signal,center,window,NumPeaks);
% "NumPeaks" = number of peaks in the model (default is 1 if not
% specified). 
% 
% peakfit(signal,center,window,NumPeaks,peakshape); 
% Specifies the peak shape of the model: "peakshape" = 1-19. (1=Gaussian
% (default), 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
% 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
% 10=sigmoid, 11=Fixed-width Gaussian, 12=Fixed-width Lorentzian;
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=BiLorentzian,
% 16=Fixed-position Gaussians; 17=Fixed-position Lorentzians;
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile.
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra) Specifies the
% value of 'extra', used only in the Voigt, Pearson, exponentionally
% broadened Gaussian, Gaussina/Lorentzian blend, and bifurcated Gaussian
% and Lorentzian shapes to fine-tune the peak shape.
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials);
% Performs "NumTrials" trial fits and selects the best one (with lowest
% fitting error). NumTrials can be any positive integer (default is 1).
%
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start)
% Specifies the first guesses vector "firstguess" for the peak positions
% and widths. Must be expressed as a vector , in square brackets, e.g.
% start=[position1 width1 position2 width2 ...]
% 
% peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero) 
% Sets autozero mode in the last argument: autozero=0 (default) does not
% subtract baseline from data segment;. autozero=1 interpolates a linear
% baseline from the edges of the data segment and subtracts it from the
% signal (assumes that the peak returns to the baseline at the edges of the
% signal); autozero=2,  like mode 1 except that it computes a quadratic
% curved baseline; autozero=3 compensates for a flat baseline without
% reference to the signal itself (best if the peak does not return to the
% baseline at the edges of the signal).
%
% [FitResults,MeanFitError]=peakfit(signal,center,window...) Returns the
% FitResults vector in the order peak number, peak position, peak height,
% peak width, and peak area), and the MeanFitError (the percent RMS
% difference between the data and the model in the selected segment of that
% data) of the best fit.
%
% [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,...)
% Prints out parameter error estimates for each peak fit.
%
% Optional output parameters 
% 1. FitResults: a table of model peak parameters, one row for each peak,
%    listing Peak number, Peak position, Height, Width, and Peak area.
% 2. LowestError: The rms fitting error of the best trial fit.
% 3. BestStart: the starting guesses that gave the best fit.
% 4. xi: vector containing 100 interploated x-values for the model peaks. 
% 5. yi: matrix containing the y values of each model peak at each xi. 
%    Type plot(xi,yi(1,:)) to plot peak 1 or plot(xi,yi) to plot all peaks
% 6. BootResults: a table of bootstrap precision results for a each peak
%    and peak parameter.

%
global AA xxx PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO

format short g
format compact
warning off all

NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  Must be isignal(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % Must be isignal(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=0;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);
ShapeString='Gaussian';

% Define values of any missing arguments
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        plots=1;
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 4
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 5
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 8
        AUTOZERO=0;
        FIXEDWIDTH=0;
        plots=1;
    case 9
        AUTOZERO=autozero;
        FIXEDWIDTH=0;
        plots=1;
    case 10
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
        plots=1;
    case 11
        AUTOZERO=autozero;
        FIXEDWIDTH=fixedparameters;
    otherwise
end % switch nargin

% Default values for placeholder zeros
if NumTrials==0;NumTrials=1;end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
if FIXEDWIDTH==0, FIXEDWIDTH=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3,end
if AUTOZERO<0,AUTOZERO=0,end

% Remove linear baseline from data segment if AUTOZERO==1
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

% Assign ShapStrings
switch peakshape
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Sigmoid';
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='BiLorentzian';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt profile';
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
           fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
              fixedstart(pc)=fixedstart(pc)+.1*randn().*fixedstart(pc);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
%         plot(zxx,zyy);pause
        TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart);     
      case 19
          TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart);
      case 20
          % extra=extra+.2*randn().*extra;
          TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH);
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        otherwise
    end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)
%
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=sigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDWIDTH);
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDWIDTH);
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BiLorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
       otherwise
  end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
heightsize=size(height');
AAsize=size(AA);
if heightsize(2)==AAsize(1),
   mmodel=height'*AA+baseline;
else
    mmodel=height*AA+baseline;
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
for m=1:NumPeaks,
    if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
    yi(m,:)=height(m)*AA(m,:); % (NEW) Place y values of individual model peaks into matrix yi
end
xi=xxx+xoffset; % (NEW) Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed lines
    if peakshape==16||peakshape==17
    else
        for marker=1:NumPeaks,
            markx=BestStart((2*marker)-1);
            subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
        end % for
    end % if peakshape
    plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
    hold off;
    axis([min(xx)+xoffset max(xx)+xoffset min(yy) max(yy)]);
    switch AUTOZERO,
        case 0
            title('peakfit 4.2   No baseline correction')
        case 1
            title('peakfit 4.2   Linear baseline subtraction')
        case 2
            title('peakfit 4.2   Quadratic subtraction baseline')
        case 3
            title('peakfit 4.2   Flat baseline correction')
    end
 
    switch peakshape
    case {4,20}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Shape Constant = ' num2str(extra) ] )
    case {5,8,18}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Time Constant = ' num2str(extra) ] )
    case 13
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     % Gaussian = ' num2str(extra) ] )
    case {14,15}
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '%     Asymmetry = ' num2str(extra) ] )
    otherwise
        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(round(100*LowestError)/100) '% ' ] )
    end

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'b.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
end % if plots

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                else
                    FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)];
                end
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models only
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                if peakshape==16||peakshape==17, % Fixed-position shapes only
                    FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                else
                    FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
                end
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% Display Fit Results on upper graph
if plots,
    subplot(2,1,1);
    startx=min(xx)+xoffset+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=(max(yy)-min(yy))./10;
    starty=max(yy)-dyy;
    FigureSize=get(gcf,'Position');
    if peakshape==9||peakshape==10||peakshape==19,  % Pulse and sigmoid shapes only
        text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
    else
        text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area'] );
    end
    % Display FitResults using sprintf
    for peaknumber=1:NumPeaks,
        for column=1:5,
            itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
    xposition=startx;
    yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
    if AUTOZERO==3,
       text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
    end
end % if plots

if NumArgOut==6,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(5,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDWIDTH);
        for peak=1:NumPeaks,
            BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:5);
        BootstrapSTD=BootstrapSTD(2:5);
        BootstrapIQR=BootstrapIQR(2:5);
        PercentRSD=PercentRSD(2:5);
        PercentIQR=PercentIQR(2:5);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==6,
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)
% 
global PEAKHEIGHTS FIXEDWIDTH FIXEDPOSITIONS AUTOZERO
format short g
format compact
warning off all
FIXEDWIDTH=fixedparameters;
xoffset=0;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;

for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'Display','off','MaxFunEvals',200 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials, 
  switch peakshape
    case 1
        TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart);
    case 2
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart);
    case 3
        TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart);
    case 4
        TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart);
    case 5
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart);
    case 6
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart);
    case 7
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitlorentziancw(lambda,xx,yy)),cwnewstart);
    case 8
        cwnewstart(1)=newstart(1);
        for pc=2:NumPeaks,
            cwnewstart(pc)=newstart(2.*pc-1);  
        end
        cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
        TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart);
      case 9
          TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart);
      case 10
          TrialParameters=fminsearch(@(lambda)(fitsigmoid(lambda,xx,yy)),newstart);
      case 11
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart);
      case 12
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart);
      case 13
          TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart);
      case 14
          TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart);
      case 15
          TrialParameters=fminsearch(@(lambda)(fitBiLorentzian(lambda,xx,yy,extra)),newstart);
      case 16
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart);
      case 17
          fixedstart=[];
          for pc=1:NumPeaks,
              fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
          end
          TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstar);
      case 18
        zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
        zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
        TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,zxx,zyy,-extra)),newstart);
      case 19
          TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart);
      case 20
        % extra=extra+.2*randn().*extra;
        TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart);
      otherwise
  end % switch peakshape
  
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch peakshape
    case 1
        A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 3
        A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 4
        A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 6
        A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
    case 7
        A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));  
    case 8
        A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';    
    case 9
        A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));  
    case 10
        A(m,:)=sigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 11
        A(m,:)=gaussian(xx,TrialParameters(m),FIXEDWIDTH);
    case 12
        A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDWIDTH); 
    case 13
        A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 14
        A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 15
        A(m,:)=BiLorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);       
    case 16
        A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 17
        A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
    case 18
        A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    case 19
        A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m)); 
    case 20
        A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra); 
   end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+randn/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+randn/10);
    end
end % for

% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3,
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end

% Compare trial model to data segment and compute the fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
end % for k (NumTrials)

for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

for m=1:NumPeaks,
    if m==1,
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12,  % Fixed-width shapes only
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    else
        if peakshape==6||peakshape==7||peakshape==8, % equal-width peak models
            FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
        else
            if peakshape==11||peakshape==12, % Fixed-width shapes only
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDWIDTH area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]];
            end
        end % if peakshape
    end % m==1
end % for m=1:NumPeaks
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPOSITIONS
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPOSITIONS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDWIDTH
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDWIDTH)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitlorentziancw(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar) 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fitlognormal(lambda,t,y)
%	Fitting function for lognormal, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlognormal assumes a lognormal function 
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lognormal(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar) 
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band signal.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
% figure(2);plot(ey);figure(1);
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponental pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS =abs(A\y');
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar) 
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% sigmiods of the form Height./(1 + exp((t1 - t)/t2))
global PEAKHEIGHTS AUTOZERO
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g=sigmoid(x,t1,t2)
% g=1./(1 + exp((t1 - x)./t2))'; % old version of sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); % Modified sigmoid
% Bifurcated sigmoid
% lx=length(x);
% hx=val2ind(x,t1);
% g(1:hx)=1./(1 + exp((t1 - x(1:hx))./t2));
% g(hx+1:lx)=1./(1 + exp((t1 - x(hx+1:lx))./(1.3*t2)));
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.
global PEAKHEIGHTS AUTOZERO
A = zeros(length(t),round(length(lambda)/2));


% Assumes that OKpeakdata.mat has been loaded into the workspace. 
% EXAMPLE 1: 
y=cos(.1:.1:100);OKpeak(y);
% EXAMPLE 2: 
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];OLKKk(M);
% EXAMPLE 3: 
x=[0:.1:100];y=(x.*sin(x)).^2;OKpeak(x,y);
% EXAMPLE 4: 
x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));OKpeak(x,y,10);
ipeak([x;y],10);
ipeak(humps(0:.01:2),3);
x=[0:.1:10];y=exp(-(x-5).^2); OKpeak([x' y'],1);
% EXAMPLE 5: 
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];ipeak(M,0,0,.0001,20,20);
% EXAMPLE 6: 
load OKpeakdata.mat;ipeak(Sample1,0,110,0.06,3,4,249.7,0.4);
% EXAMPLE 7: 
OKpeak(Sample1,0,110,0.06,3,4,249.7,0.4,1);
% EXAMPLE 8: 
OKpeak(Sample1,0,100,0.05,3,6,296,5,0.1,Positions,Names);
% OKPeak demos
OKdemo
OKdemo1
OKdemo2
OKdemo3
OKdemo4
OKdemo5



