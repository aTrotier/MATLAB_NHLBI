function [xm0, ym0, g_info]=vds_M0(smax,gmax,dt,interleaves,FOV,kr,trapflag)
% [xm0,ym0]=vds_M0(smax,gmax,dt,interleaves,FOV,kr,<trapflag=0,1>)
% Default to triangular rewinder, otherwise use MCR's trapz (Aug 2018).
%
% ### DEPENDENT FUNCTIONS ###
% add_rewinders.m, balance_moments.m
%
% R Ramasawmy, NHLBI, Oct 2018.
if nargin == 0
    xm0.add_rewinders = @add_rewinders;
    xm0.balance_moments = @balance_moments;
else
    if nargin < 7
        trapflag = 0;
    end
    
    % --- Run VDS at Gradient raster time ---
    [~,g] = vds_GRT(smax,gmax,dt,interleaves,FOV,kr);
    xgrad = (real(g(:))*cos(0)+imag(g(:))*sin(0));
    ygrad = (-real(g(:))*sin(0)+imag(g(:))*cos(0));
    
    % --- use scanner gradient settings ---
    Tgsample = 1e-5; % GRT
    amppersamp = (1/sqrt(3.0)) * smax * Tgsample;
    
    [xm0, ym0, g_info]=add_rewinders(xgrad,ygrad,amppersamp, gmax, trapflag);
end
% --- preview traj ---

% figure, plot(cumsum(xm0),cumsum(ym0))

end


% Ramasawmy
% translation of vds.cpp
% Matlab specific: Added in "trapflag" in order to switch between triagular
% and trapezoidal M0 balancing (traps code: MCR, Aug 2018).

% bool add_rewinders(double **xgrad,
% 			double **ygrad,
% 			int *numgrad,
% 			int numgradmax,
% 			double amppersamp,
% 			double gradmax,
% 			int maxmoment)
% {
function [xgrad_m0, ygrad_m0, g_info] = add_rewinders(xgrad,ygrad,amppersamp,gradmax, trapflag)
%    int numx=*numgrad;
%    int numy=*numgrad;
%    int numorig=*numgrad;
%    int a;
%    double gend;
%    bool success;
%
%    std::cout << "ADD_REWINDERS:  numgrad = " << numx << ", amppersamp = " << amppersamp << " G/cm/sample" << "maxmoment = " << maxmoment << std::endl;
%
%
%    // Add rampdown on both axes
%    gend = std::max (  fabs((*xgrad)[numorig-1]) , fabs((*ygrad)[numorig-1]) );
%    int rampdown = ceil(gend/amppersamp);
%    std::cout << "***rbhpg gend = "<<gend << " amppersamp = " << amppersamp << " rampdown = " <<rampdown << std::endl;
%
%    for (a=0;a<rampdown;a++) {
%       (*xgrad)[numorig+a] = (*xgrad)[numorig-1] * (rampdown-1-a)/rampdown;
%       (*ygrad)[numorig+a] = (*ygrad)[numorig-1] * (rampdown-1-a)/rampdown;
%    }
%    numx += rampdown;
%    numy += rampdown;

% % ### RR ###
% % Tgsample = 1e-5;
% % slewmax = 14414.4; % fast
% % amppersamp = (1/sqrt(3.0)) * slewmax * Tgsample;
% % gmax = 2.4;
gend = max([abs(xgrad(end)),abs(ygrad(end))]);
rampdown = ceil(gend/amppersamp);

rampdown_x = zeros(rampdown,1);
rampdown_y = zeros(rampdown,1);
for a=1:rampdown
    rampdown_x(a) = (xgrad(end)) * (rampdown-1-(a-1))/rampdown;
    rampdown_y(a) = (ygrad(end)) * (rampdown-1-(a-1))/rampdown;
end

% ### RR ###
g_info.n_spiral = length(xgrad);
g_info.n_rampdown = rampdown;
% ###    ###

% % figure, plot([[xgrad; rampdown_x] [ygrad; rampdown_y]])

%    //now calc moment balancing starting from the end of the rampdown
%    int numspiralendrampdown=numx; //now same as numy
%    //The last argument 0 allows balance_moments to calculate what is needed for each spiral channel.
%    //Later we take the longest and design both based on this. So that the rephas grads have same timing.
%
%    int ret = balance_moments(ygrad,&numy,amppersamp,gradmax,maxmoment,numgradmax,0);
%    if(ret== -1) {
% 		std::cout << "   failed...exceeded gradient array size NUMGRADMAX" << std::endl;
% 		return false; //to signal fSeqPrep that this spiral with rewinding is now too long
% 		}
%
%    ret = balance_moments(xgrad,&numx,amppersamp,gradmax,maxmoment,numgradmax,0);
%    if(ret== -1) {
% 		std::cout << "   failed...exceeded gradient array size NUMGRADMAX" << std::endl;
% 		return false; //to signal fSeqPrep that this spiral with rewinding is now too long
% 		}

xgrad1 = [xgrad; rampdown_x];
ygrad1 = [ygrad; rampdown_y];

[gx1] = balance_moments(xgrad1,amppersamp,gradmax,0, trapflag);
[gy1] = balance_moments(ygrad1,amppersamp,gradmax,0, trapflag);


% 	//choose the slowest one and recalculate the rephase so both axes have same timing, it's just neater.
% 	*numgrad = std::max(numx,numy);
% 	int nramp = (*numgrad - numspiralendrampdown) / 2;
%
% 	numx = numspiralendrampdown; //numx and numy were modified by balance_moments above
% 	numy = numspiralendrampdown;
%
%    ret = balance_moments(ygrad,&numy,amppersamp,gradmax,maxmoment,numgradmax,nramp);
%    if(ret== -1) {
% 		std::cout << "   failed...exceeded gradient array size NUMGRADMAX" << std::endl;
% 		return false; //to signal fSeqPrep that this spiral with rewinding is now too long
% 		}
%
%    ret = balance_moments(xgrad,&numx,amppersamp,gradmax,maxmoment,numgradmax,nramp);
%    if(ret== -1) {
% 		std::cout << "   failed...exceeded gradient array size NUMGRADMAX" << std::endl;
% 		return false; //to signal fSeqPrep that this spiral with rewinding is now too long
% 		}

numgrad = max(length(gx1),length(gy1));
if ~trapflag
    nramp = (numgrad) / 2;
else
    nramp = numgrad;
end

[gx2] = balance_moments(xgrad1,amppersamp,gradmax,nramp, trapflag);
[gy2] = balance_moments(ygrad1,amppersamp,gradmax,nramp, trapflag);


%    if (numx != *numgrad) {
% 	   std::cout << "error in addrewinders numx = " << numx << " should be " << *numgrad << std::endl;
%    }
%
%    if (numy != *numgrad) {
% 	   std::cout << "error in addrewinders numy = " << numx << " should be " << *numgrad << std::endl;
%    }
% ### RR laziness

% ### RR ###

xgrad_m0 = [xgrad1(:); gx2(:)];
ygrad_m0 = [ygrad1(:); gy2(:)];

g_info.n_M0 = length(gx2);

figure,
% subplot(1,2,1);
plot(xgrad_m0), hold on, plot(ygrad_m0); legend({'G_x','G_y'});
% subplot(1,2,2);
% plot(diff([xgrad1(:); gx2(:)])), hold on, plot(diff([ygrad1(:); gy2(:)])); legend({'G_x','G_y'});

% ###    ###

end

% Ramasawmy
% translation of vds.cpp
% Matlab specific: Added in "trapflag" in order to switch between triagular
% and trapezoidal M0 balancing (traps code: MCR, Aug 2018).

% int balance_moments(double **g,
% 			int *num,
% 			double amppersamp,
% 			double gradmax,
% 			int maxmoment,
% 			int numgradmax,
% 			int nrampreq)
%(
function [g_bal] = balance_moments(g, amppersamp, gradmax, nrampreq, trapflag)

%
%     int onum = (*num);
% 	double m0;
% 	double nramp,sramp, gpeak;
% 	int i, nrampc;
%   int nflatc *MCR*
%
% 	std::cout << "Entering balance_moments: input waveform length = " << onum << std::endl;
% 	std::cout << "Entering balance_moments: input amppersamp      = " << amppersamp << std::endl;
num = length(g);
onum = num;
disp(['Entering balance_moments: input waveform length = ' num2str(onum)]);
disp(['Entering balance_moments: input amppersamp = ' num2str(amppersamp)]);

% 	m0 = 0;
% 	for (i=0;i<*num;i++) {
% 	   m0 += (*g)[i];
% 	}
% 	std::cout << "    initial M0 = " << m0 << std::endl;

m0 = sum(g);

disp(['initial M0 = ' num2str(m0) ]);

if ~trapflag
    %% ### TRIANGLE M0 ###
    %%
    disp('*** Using triangular balancers');
    % if (m0!=0.0) {
    % 			nramp = sqrt (abs(m0)/amppersamp);
    % 			gpeak = -(m0/fabs(m0)) * nramp * amppersamp;
    % 			// Recalc gpeak for next slowest multiple of 10us ramp
    % 			nrampc = ceil(nramp);
    % 			if (nrampc==0) {
    % 				gpeak = 0.0;
    % 				nrampc = 1;
    % 			}
    % 			gpeak = gpeak * (nramp/(double)nrampc);
    % 			sramp = gpeak/(double)nrampc;
    % 			std::cout << "***rbhpg zeroth moment rephasing nramp =" << nramp << " nrampc ="<< nrampc << " sramp =" << sramp << " gpeak ="<< gpeak << std::endl;
    % 	} else {
    % 			nramp = 0;
    % 			nrampc = 0;
    % 			sramp = 0.0;
    % 			gpeak = 0.0;
    % 	}
    
    if (m0 ~= 0.0)
        nramp = sqrt (abs(m0)/amppersamp);
        gpeak = -(m0/abs(m0)) * nramp * amppersamp;
        % Recalc gpeak for next slowest multiple of 10us ramp
        nrampc = ceil(nramp);
        if (nrampc==0)
            gpeak = 0.0;
            nrampc = 1;
        end
        gpeak = gpeak * (nramp/nrampc);
        sramp = gpeak/nrampc;
        disp(['*** zeroth moment rephasing nramp = '  num2str(nramp) ' nrampc = ' num2str(nrampc) ' sramp = ' num2str(sramp) ' gpeak = ' num2str(gpeak) ]);
    else
        nramp = 0;
        nrampc = 0;
        sramp = 0.0;
        gpeak = 0.0;
    end
    
    % 	// MCR update 20180823
    % 	if (fabs(gpeak)>gradmax) {
    % 			nramp = (double)nrampc*fabs((gpeak)/gradmax);
    % 			gpeak = fabs(gpeak)/gpeak*gradmax;
    % 			nrampc = ceil(nramp);
    % 			gpeak = gpeak* (nramp/(double)nrampc);
    % 			sramp = gpeak/(double)nrampc;
    % 	}
    % 	// MCR update 20180823 <<<<
    
    if (abs(gpeak)>gradmax)
        nramp = nrampc*abs((gpeak)/gradmax);
        gpeak = abs(gpeak)/gpeak*gradmax;
        nrampc = ceil(nramp);
        gpeak = gpeak* (nramp/nrampc);
        sramp = gpeak/nrampc;
    end
    
    % 	//the calling function uses this to slow down the rephase pulse of the faster axis so that they have same duration
    % 	if (nrampreq>0) {
    % 		if (nrampreq<nrampc) {
    % 			std::cout << "***rbhpg balance_moments in vdsMiniIRT.cpp Unexpected request to speed up rephase pulse, should be slowing it down, nrampreq =" << nrampreq << " nrampc ="<< nrampc << std::endl;
    % 			}
    % 		gpeak = gpeak * (double)nrampc / (double)nrampreq;
    % 		sramp = gpeak / (double)nrampreq;
    % 		nrampc=nrampreq;
    % 		std::cout << "***rbhpg balance_moments in vdsMiniIRT.cpp slowed to nrampreq =" << nrampreq << " sramp ="<< sramp << std::endl;
    % 	}
    if (nrampreq>0)
        if (nrampreq<nrampc)
            disp(['*** error']);
        end
        gpeak = gpeak * nrampc /nrampreq;
        sramp = gpeak / nrampreq;
        nrampc=nrampreq;
        disp(['*** nrampreq = ' num2str(nrampreq) ' sramp = ' num2str(sramp) ]);
    end
    % 		// rbhpg Prevent writing past end of input arrays return -1 value if so
    % 	int numgradrephased= *num + 2*nrampc;
    % 	if (numgradrephased > numgradmax) {
    % 			printf("From balance_moments  new gradient total points %d exceeds array allocation NUMGRADMAX %d\n",numgradrephased,numgradmax);
    % 			printf("From balance_moments  returning -1\n");
    % 			return (-1);
    % 	}
    % ### RR error checking, too lazy to include atm
    
    % 	// Write the rephasing pulse into the end of the arbitrary grad array
    % 	for (i=0;i<nrampc;i++) (*g)[*num+i] = sramp*(i+1);
    % 	*num += nrampc;
    % 	for (i=0;i<nrampc;i++) (*g)[*num+i] = sramp*(double)(nrampc-1-i);
    % 	*num += nrampc;
    
    g_bal = zeros(1,2*nrampc);
    for i=1:nrampc
        %     g(num+i) = sramp*(i);
        g_bal(i)  = sramp*(i);
    end
    num = num + nrampc;
    
    for i=1:nrampc
        %     g(num+i) = sramp*(nrampc-1-i);
        g_bal(nrampc+i)  = sramp*(nrampc-i);
    end
    num = num + nrampc;
    
else
    %% ### TRAPEZOID M0 ###
    %%
    disp('*** Using trapezoidal balancers');
    %     	nramp = 0;
    % 	if (m0!=0.0) {
    % 			//t1 = gradmax/amppersamp;
    % 			//t2 = m0/gradmax-t1;
    %
    % 			nrampc = ceil((gradmax/1.5)/amppersamp);
    % 			nflatc = ceil(fabs(m0)/(gradmax/1.5)-nrampc);
    % 			//gpeak = -(m0/fabs(m0)) * nramp * amppersamp;
    % 			// Recalc gpeak for next slowest multiple of 10us ramp
    % 			//nrampc = ceil(nramp);
    % 			if (nrampc==0) {
    % 				gpeak = 0.0;
    % 				nrampc = 1;
    % 			}
    % 			nflatc = std::max(nflatc,0);
    % 			gpeak = -m0/(nrampc+nflatc);
    % 			sramp = gpeak/(double)nrampc;
    % 			std::cout << "***rbhpg zeroth moment rephasing nflatc =" << nflatc << " nrampc ="<< nrampc << " sramp =" << sramp << " gpeak ="<< gpeak << std::endl;
    % 	} else {
    % 			nramp = 0;
    % 			nflatc = 0;
    % 			nrampc = 0;
    % 			sramp = 0.0;
    % 			gpeak = 0.0;
    % 	}
    if (m0 ~= 0.0)
        nrampc = ceil((gradmax/1.5)/amppersamp);
        nflatc = ceil(abs(m0)/(gradmax/1.5)-nrampc);
        
        if (nrampc==0)
            gpeak = 0.0;
            nrampc = 1;
        end
        nflatc = max(nflatc,0);
        gpeak = -m0/(nrampc+nflatc);
        sramp = gpeak/nrampc;
        
        disp(['*** zeroth moment rephasing nflatc = '  num2str(nflatc) ' nrampc = ' num2str(nrampc) ' sramp = ' num2str(sramp) ' gpeak = ' num2str(gpeak) ]);
    else
        nramp = 0;
        nflatc = 0;
        nrampc = 0;
        sramp = 0.0;
        gpeak = 0.0;
    end
    
    %     	//the calling function uses this to slow down the rephase pulse of the faster axis so that they have same duration
    % 	if (nrampreq>0) {
    % 		if (nrampreq<nrampc) {
    % 			std::cout << "***rbhpg balance_moments in vdsMiniIRT.cpp Unexpected request to speed up rephase pulse, should be slowing it down, nrampreq =" << nrampreq << " nrampc ="<< nrampc << std::endl;
    % 			}
    % 		nflatc = nrampreq-2*nrampc;
    % 		gpeak = -m0/(nrampc+nflatc);
    % 		sramp = gpeak/(double)nrampc;
    % 		/*gpeak = gpeak * (double)nrampc / (double)nrampreq;
    % 		sramp = gpeak / (double)nrampreq;
    % 		nrampc=nrampreq;*/
    % 		std::cout << "***rbhpg balance_moments in vdsMiniIRT.cpp slowed to nrampreq =" << nrampreq << " sramp ="<< sramp << std::endl;
    % 	}
    
    if (nrampreq>0)
        if (nrampreq<nrampc)
            disp(['*** error']);
        end
        nflatc = nrampreq-2*nrampc;
        gpeak = -m0/(nrampc+nflatc);
        sramp = gpeak/nrampc;
        disp(['*** nrampreq = ' num2str(nrampreq) ' sramp = ' num2str(sramp) ]);
    end
    
    
    % 	// rbhpg Prevent writing past end of input arrays return -1 value if so
    % 	int numgradrephased= *num + 2*nrampc + nflatc;
    % 	if (numgradrephased > numgradmax) {
    % 			printf("From balance_moments  new gradient total points %d exceeds array allocation NUMGRADMAX %d\n",numgradrephased,numgradmax);
    % 			printf("From balance_moments  returning -1\n");
    % 			return (-1);
    % 	}
    % ### RR error checking, too lazy to include atm
    %
    %
    % 	// Write the rephasing pulse into the end of the arbitrary grad array
    % 	for (i=0;i<nrampc;i++) (*g)[*num+i] = sramp*(i+1);
    % 	*num += nrampc;
    % 	for (i=0;i<nflatc;i++) (*g)[*num+i] = sramp*nrampc;
    % 	*num += nflatc;
    % 	for (i=0;i<nrampc;i++) (*g)[*num+i] = sramp*(double)(nrampc-1-i);
    % 	*num += nrampc;
    
    g_bal = zeros(1,2*nrampc + nflatc);
    for i=1:nrampc
        %     g(num+i) = sramp*(i);
        g_bal(i)  = sramp*(i);
    end
    num = num + nrampc;
    
    for i=1:nflatc
        %     g(num+i) = sramp*(i);
        g_bal(nrampc+ i)  = sramp*nrampc;
    end
    num = num + nflatc;
    
    for i=1:nrampc
        %     g(num+i) = sramp*(nrampc-1-i);
        g_bal(nrampc+nflatc+ i)  = sramp*(nrampc-i);
    end
    num = num + nrampc;
    
    
end

% 	// Confirm achieved M0 rephasing
%
% 	m0 = 0;
% 	for (i=0;i<*num;i++) {
% 	   m0 += (*g)[i];
% 	}
% 	std::cout << "    final M0 = " << m0 << std::endl;
% 	std::cout << "    new length = " << *num << std::endl;
%
% 	return (1);
% }

m0 = sum([g(:); g_bal(:)]);

disp('##########################');
disp(['final M0 = ' num2str(m0) ]);
disp(['new length = ' num2str(num) ]);
disp('##########################');

%% RR :: Plotting space
% figure, plot(g_bal);

end