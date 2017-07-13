function gauss = gcoeff3m(debiasedhalldata,probepos);
if ~(nargin==2)
	error('Usage: gausscoeff = gcoeff3m(debiasedhalldata,probepos)');
end

probemask = ones(31,1);

%Modified: 20 Mar 2015 (adapted from getcoeff2014 - stored on lab server)

%revamped by Onur Kara, 4 Jun 2015 to allow for more coefficients (using the legendre polynomials of Matlab)

%Usage: getcoeff(debiasedhalldata,probepos,[probemask]), where 'debiasedhalldata is
%assumed to have 41 columns, the first one being time and the remaining
%40 to have the ordering as specified in the file
%'probeOrderinginMatlab.txt', in the folder data/sixtycm on wave.umd.edu;
%this data is assumed to have the probe offsets, due to both the probe
%biases as well as the applied magnetic field, removed, so that it 
%represents the actual signal of the induced field; if hallfile is of the
%form 'filename', this program attempts to load the appropriate file from
%disk, while if it is of the form filename, the program looks in the
%current workspace for the variable filename; the optional argument
%'probemask' is expected to be a vector with 30 entries, each entry being
%either a 1 or a 0, with 1 indicating the corresponding hallprobe is alive
%and 0 indicating it is dead; if 'probemask' is provided, the program uses
%it to determine which probe data to discard, and calculates matrix
%required to perform the least squares fit; if 'probemask' is not provided,
%the program assumes that all probes are alive; if no arguments are given,
%the program outputs usage information


%Purpose: This program calculates the gauss coefficients for hall data
%taken using the new set-up in Matlab, and outputs gauss coefficients up 
%through l=4, m=4




% argument checking and handling

if ((nargin < 1) || (nargin > 3))
    error('Usage: getcoeff(debiasedhalldata,[probemask])');
end;

% *** HCN ***
% number of Gauss coefficients to retrieve
 %lmax = input('Give maximum degree to invert for (l<=4) \n');
 %ncoef = lmax*(lmax+2); % (no monopole, l=0)
%
lmax = 4
ncoef = lmax*(lmax+2);


if(ischar(debiasedhalldata))
    datafile = load(debiasedhalldata);
else
    datafile = debiasedhalldata;
end;

[dumy nprobes] = size(datafile);
probemasksize = [nprobes 1];

if (nargin == 3) % if probemask provided, check that it has correct form
    if (~isequal((size(probemask)),probemasksize) && ...
            ~isequal((size(probemask')),probemasksize))
        error('Input: argument "probemask" must be a 30 element vector');
    elseif ((max(probemask) > 1) || (min(probemask) < 0))
        error('Input: argument "probemask" should consist of 1s and 0s');
    end;
elseif (nargin == 2)   % if probemask not specified, use default value
    probemask = ones(nprobes,1); 
end;

% convert nonzero elements of probemask to ones (in case they aren't), and
% then sum them up to get number of probes being used, put in variable
% 'nprobes'

nprobes = sum((probemask & ones(size(probemask,1),size(probemask,2))));

nsamples = size(datafile,1);

% load necessary files, define necessary parameters

scalefac = 839.5; % from a long-ago calibration?

% probe positions (copied from probepos2014_allprobes.dat, located on wave
% in /data/zero)

r(1:nprobes,1) = 0; % preallocate r, theta, and phi
theta(1:nprobes,1) = 0;
phi(1:nprobes,1) = 0;
i = 1; 
for j=1:nprobes
    if(probemask(j))
        r(i,1) = probepos(j,1);
        theta(i,1) = probepos(j,2);
        phi(i,1) = probepos(j,3);
        i = i + 1;
    end;
end;


% use r, theta, and phi to construct the basis functions

% *** HCN : Matlab Schmitt normalization x (-1)^m, it seems ****
% get Legendre associated polynomials for all l, calculated at the probes' co-latitudes
f=zeros(ncoef,nprobes); % f(coeffnum,probenum)
ii=0;
for ll=1:lmax
	Plm(1:ll+1,1:nprobes) = legendre(ll,transpose(cos(theta)),'sch'); % associated Legendre polynomials of degree l and order m=0 to l (warning : m=0 for index=1) (Schmidt normalization used by Sisan)
	for mm=1:2:ll
		Plm(mm+1,:) = -Plm(mm+1,:); % change sign of odd m (except 0)
	end
	for mm=0:ll
		b_r(:,1) = ll*(ll+1) * r.^(-(ll+2)) .* transpose(Plm(mm+1,:)); % radial component of the magnetic field at all probe positions for given l and m
		dPlm_dtheta_m(:,1) = mm * cot(theta) .* transpose(Plm(mm+1,:)); % theta derivative of Plm when l=m
		if (ll ~= mm)
			if (mm==0)
				coef = sqrt(ll*(ll+1)/2); % no minus sign because (-1)^m Schmidt
			else
				coef = sqrt((ll-mm)*(ll+mm+1)); % no minus sign because (-1)^m Schmidt
			end
			dPlm_dtheta_m(:,1) = dPlm_dtheta_m(:,1) + coef*transpose(Plm(mm+2,:)); % added term when l ~= m (with normalization for Schmidt*(-1)^m)
		end
		b_theta(:,1) = -r.^(-(ll+2)) .* dPlm_dtheta_m ; % orthoradial component of the magnetic field at all probe positions for given l and m
		ii = ii+1; % index of the 'f' function
		%f(ii,:) = b_r .* sin(theta) + b_theta .* cos(theta); % that's it for m=0 
		f(ii,:) = b_r; %trying just radial comp
        if (mm > 0)
			truc(1:nprobes,1) = f(ii,1:nprobes);
			f(ii,:) = truc .* cos(mm*phi); % cosine phi term
			ii = ii+1; % index of the 'f' function
			f(ii,:) = truc .* sin(mm*phi); % sine phi term
		end
	end
end
% ****
disp(['Identify Coefficients HCN : ' int2str(ii)])
% *****

% ****** check with previous expressions ******
if (0 == 1)
	f_ori=zeros(ncoef,nprobes); % f(coeffnum,probenum)
	f_ori(1,:)=3*r.^-3.*cos(theta).*sin(theta);
	f_ori(2,:)=r.^-3.*(3*(cos(theta).^2)-2).*abs(sin(theta))./sin(theta).*cos(phi);
	f_ori(3,:)=r.^-3.*(3*(cos(theta).^2)-2).*abs(sin(theta))./sin(theta).*sin(phi);
	if (lmax > 1)
		f_ori(4,:)=3/4*r.^-4.*sin(theta).*(3+5*cos(2*theta));
		f_ori(5,:)=sqrt(3)*r.^-4.*(5*(cos(theta).^2)-4).*cos(theta).* ...
		    abs(sin(theta))./sin(theta).*cos(phi);
		f_ori(6,:)=sqrt(3)*r.^-4.*(5*(cos(theta).^2)-4).*cos(theta).* ...
		    abs(sin(theta))./sin(theta).*sin(phi);
		f_ori(7,:)=sqrt(3)/2*r.^-4.*sin(theta).*(3-5*(cos(theta).^2)).*cos(2*phi);
		f_ori(8,:)=sqrt(3)/2*r.^-4.*sin(theta).*(3-5*(cos(theta).^2)).*sin(2*phi);
	end
	if (lmax > 2)
		f_ori(9,:)=5/16*r.^-5.*(2*sin(2*theta)+7*sin(4*theta));
		f_ori(10,:)=1/16*sqrt(3/2)*r.^-5.*sin(theta)./abs(sin(theta)).* ... % from final.nb
		    (35*cos(4*theta)-3).*cos(phi);
		f_ori(11,:)=1/16*sqrt(3/2)*r.^-5.*sin(theta)./abs(sin(theta)).* ... % from final.nb
		    (35*cos(4*theta)-3).*sin(phi);
		f_ori(12,:)=sqrt(15)/2*r.^-5.*sin(theta).*cos(theta).* ...
		    (5-7*(cos(theta).^2)).*cos(2*phi);
		f_ori(13,:)=sqrt(15)/2.*r.^-5.*sin(theta).*cos(theta).* ...
		    (5-7*(cos(theta).^2)).*sin(2*phi);
		f_ori(14,:)=1/4*sqrt(5/2)*r.^-5.*sin(theta).*abs(sin(theta)).* ...
		    (7*cos(2*theta)-1).*cos(3*phi);
		f_ori(15,:)=1/4*sqrt(5/2)*r.^-5.*sin(theta).*abs(sin(theta)).* ...
		    (7*cos(2*theta)-1).*sin(3*phi);
	end
	if (lmax > 3)
		f_ori(16,:)=15/128*r.^-6.*(2*sin(theta)+7*(sin(3*theta)+3*sin(5*theta)));
		f_ori(17,:)=1/2*sqrt(5/2)*r.^-6.*sin(theta)./abs(sin(theta)).*cos(theta).* ... % from final.nb
		    (18-77*(cos(theta).^2)+63*(cos(theta).^4)).*cos(phi);
		f_ori(18,:)=1/2*sqrt(5/2)*r.^-6.*sin(theta)./abs(sin(theta)).*cos(theta).* ... % from final.nb
		    (18-77*(cos(theta).^2)+63*(cos(theta).^4)).*sin(phi);
		f_ori(19,:)=-sqrt(5)/32*r.^-6.*sin(theta).*(5+28*cos(2*theta)+ ...
		    63*cos(4*theta)).*cos(2*phi);
		f_ori(20,:)=-sqrt(5)/32*r.^-6.*sin(theta).*(5+28*cos(2*theta)+ ...
		    63*cos(4*theta)).*sin(2*phi);
		f_ori(21,:)=3/2*sqrt(35/2)*r.^-6.*sin(theta).*cos(theta).*abs(sin(theta)).* ... % from final.nb
		    (3*(cos(theta).^2)-2).*cos(3*phi);
		f_ori(22,:)=3/2*sqrt(35/2)*r.^-6.*sin(theta).*cos(theta).*abs(sin(theta)).* ... % from final.nb
		    (3*(cos(theta).^2)-2).*sin(3*phi);
		f_ori(23,:)=-sqrt(35)/16*r.^-6.*(sin(theta)).^3.*(9*cos(2*theta)-1).* ...
		    cos(4*phi);
		f_ori(24,:)=-sqrt(35)/16*r.^-6.*(sin(theta)).^3.*(9*cos(2*theta)-1).* ...
		    sin(4*phi);
	end

% 
	figure('name','ori')
	imagesc(f_ori); caxis([-1.2 1.2]); colorbar;
	figure('name','HCN')
	imagesc(f); caxis([-1.2 1.2]); colorbar;
end
% ******************


% construct matrix from basis functions

M(1:ncoef,1:ncoef) = 0;

for i=1:ncoef
    for j=1:ncoef
        M(i,j)=(f(i,:)*f(j,:)');
    end;
end;

% grab gauss array data from input, put it in order (G1-1 ... G1-7, G2-1 -
% G2-8, G3-1 - G3-8, G4-1 - G4-7); also grab sample times

%

%
%t = datafile(:,1);
%
%data(1:nsamples,30) = 0;
%
%%ring 1
%for k=1:7
%    data(:,k) = datafile(:,k+25);
%end;
%
%%ring 2
%for k=1:8
%    data(:,k+7) = datafile(:,k+9);
%end;
%
%%ring 3
%for k=1:8
%    data(:,k+15) = datafile(:,k+17);
%end;
%
%%ring 4
%data(:,24) = datafile(:,41);
%for k=1:6
%    data(:,k+24) = datafile(:,k+34);
%end;
%
%% apply probemask to get rid of any data from a dead probe
%
%d(1:nsamples,1:nprobes) = 0;
%i = 1;
%for j=1:30
%    if(probemask(j))
%        d(:,i) = data(:,j);
%        i = i + 1;
%    end;
%end;

d = datafile;
d=d'; % put data in correct format for matrix operations

% calculate gauss coefficients; note: using Matlab's matrix division
% operation to solve system of equations Bf = M*g for g, instead of
% previous method of inverting M, and then calculating g = M^(-1)*Bf (this
% is what Matlab's documentation recommends: g = M\Bf should give faster,
% and numerically more accurate result; it uses Gaussian elimination)

Bf=f*d; % w(coeffnum,nsamples)=f(coeffnum,nprobes)*d(nprobes,nsamples)

g=zeros(size(Bf)); % preallocate storage for calculation result
for ii=1:nsamples % loop over samplenum
%     g(:,ii)=m*w(:,ii); % matrix on vector gives vector
      g(:,ii) = M\Bf(:,ii);
end


g=g'; 

% g2=g*scalefac;
gauss=g; %removed scalefac

% gauss(1:nsamples,ncoef+1) = 0;
% gauss(:,1) = t;
% gauss(:,2:ncoef+1)=g2;

% save the data as the variable 'gauss' in the file gaussfile 
    
% save(gaussfile,'gauss');  % save the gauss coefficients data in the
                          % variable 'gauss' in the file gaussfile

if (0 == 1)
%	size(g)
%	size(d)
%	size(f)
% ------------------------------
% essai CS avec inversion minimisant norme l1
	AA0 = g'; % initial vector for l1 iterations
	dumy = 0;
%	epsilon = std(d)/2.
	epsilon = std(d)/10.
%	g = l1qc_logbarrier(AA0, f', dumy, d, epsilon); % Candes'l1magic function (last arguments -> default)
%	g = LassoPrimalDualLogBarrier(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoActiveSet(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
	g = LassoBlockCoordinate(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoConstrained(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoGaussSeidel(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoGrafting(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoIteratedRidge(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoIterativeSoftThresholding(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoNonNegativeSquared(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoProjection(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoShooting(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoSignConstraints(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoSubGradient(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
%	g = LassoUnconstrainedApx(f', d, epsilon); % Schmidt's lasso function (last arguments -> default)
	gauss=g*scalefac;
% ------------------------------

end









