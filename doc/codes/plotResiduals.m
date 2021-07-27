%
%  Plot residuals form WaveHoltz iterations
%
%  cd /Users/henshaw/Dropbox/research/cgwave/doc/codes
%
function plotResiduals(varargin)

  clearvars -except varargin;

  clearvars -except varargin;
  % clear global;

  % --- Clear all open figures ----
  clearOpenFigures(1:3);
  delete(findall(gcf,'type','annotation'));  % clear annotations

%   FigList = findall(groot, 'Type', 'figure');
%   for iFig = 1:numel(FigList)
%       try
%           clf(FigList(iFig));
%       catch
%           % Nothing to do
%       end
%   end

% Set defaults for plotting 
fontSize=16; lineWidth=2; markerSize=6; 
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);
set(0,'DefaultLegendFontSize',fontSize);
% size of matlab figure : default = [560,420] 
xwidth = 560;
ywidth = 540; % 570; % 420; % 570;  


caseOption = 'trigHelmholtzCIC';
caseOption = 'gaussianSquareOmega17p17';
caseOption = 'gaussianSquare128O4Omega15';
caseOption = 'gaussianAnnulusO4Omega9p777'; 
caseOption = 'gaussianAnnulusO4Omega18p4';
caseOption = 'gaussianDiskO4Omega8p1';
caseOption = 'gaussianBoxO4Omega8p1'; 
caseOption = 'gaussianSphereO2Omega8p1';

plotConvergenceRate=0;
domain='square'; % 'annulus'


  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};

    caseOption = getString( line,'-caseOption',caseOption );
 
    % ms                     = getString( line,'-ms',ms );
    % plotOption             = getInt( line,'-plotOption',plotOption );
    % numResolutions         = getInt( line,'-numResolutions',numResolutions );
    % N0                     = getInt( line,'-N0',N0 );
% 
    % numFixPointIterations  = getInt( line,'-numFixPointIterations',numFixPointIterations );
    % numFixPointIterations  = getInt( line,'-numFPI'               ,numFixPointIterations ); % short form
% 
    % numFreq                = getInt( line,'-numFreq',numFreq );
        % 
    % numPer1                = getInt( line,'-numPer1',numPer1 );    
    % useFixPoint            = getInt( line,'-useFixPoint',useFixPoint ); 
% 
    % omega1                 = getReal( line,'-omega1',omega1 );       
    % omega2                 = getReal( line,'-omega2',omega2 ); 
    % omega3                 = getReal( line,'-omega3',omega3 ); 
% 
    % kx1                    = getReal( line,'-kx1',kx1 );       
    % kx2                    = getReal( line,'-kx2',kx2 );   
    % kx3                    = getReal( line,'-kx3',kx3 );   
% 
    % c                      = getReal( line,'-c',c );   
% 
    % par.savePlots          = getInt( line,'-savePlots',par.savePlots );
    % par.plotName           = getString( line,'-plotName',par.plotName );        

  end




nd=2;  % Number of space dimensions
if( strcmp(caseOption,'trigHelmholtzCIC')==1 )

  % -- Sine helmholtz solution : CIC grid ----
  plotName = 'trigHelmholtzCIC';

  trigHelmholtzCIC4FixedPoint
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  trigHelmholtzCIC4Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  c=1;                % fix me -- add to .m file 

  myTitle = sprintf('CgWaveHoltz: CIC Trig, \\omega=%.4g, L_{2h}-res',omega);

  % semilogy( itv1,res1,'r-x', itv2,res2,'b-o');
  % grid on;
  % title(sprintf('CgWaveHoltz: CIC Trig, \\omega=%.4g, L_{2h}-res',omega));
  % legend('FP','GMRES');
  % xlabel('iterations'); 
% 
  % pause; pause;

  
elseif( strcmp(caseOption,'gaussianSquareOmega17p17')==1 )

  % -- Gaussian source----
  plotName = 'gaussianSquareOmega17p17';

  gaussianWaveHoltzSquare64Order4FixedPoint
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWaveHoltzSquare64Order4Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  % c=1;                % fix me -- add to .m file 

  plotConvergenceRate=1; % we can eval the true CR
  boxDims = [1,1,1]; % size of square is 1x1

  myTitle = sprintf('CgWaveHoltz: Gaussian, Square, \\omega=%.4g, L_{2h}-res',omega);

elseif( strcmp(caseOption,'gaussianSquare128O4Omega15')==1 )

  % -- Gaussian source----
  plotName = 'gaussianSquare128O4Omega15';

  gaussianWHSquare128O4Omega15FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHSquare128O4Omega15Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  boxDims = [1,1,1]; % size of square is 1x1


  myTitle = sprintf('CgWaveHoltz: Gaussian, Square, \\omega=%.4g, L_{2h}-res',omega);

elseif( strcmp(caseOption,'gaussianAnnulusO4Omega9p777')==1 )

  % -- ANNULUS : Gaussian source----
  plotName = 'gaussianAnnulusO4Omega9p777';

  gaussianWHAnnulusO4Omega9p777FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHAnnulusO4Omega9p777Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  domain='annulus'; 
  annulusDims = [.5,1]; % inner and outer radii

  myTitle = sprintf('CgWaveHoltz: Gaussian, Annulus, \\omega=%.4g, L_{2h}-res',omega);

elseif( strcmp(caseOption,'gaussianAnnulusO4Omega18p4')==1 )

  % -- ANNULUS : Gaussian source----
  plotName = 'gaussianAnnulusO4Omega18p4';

  gaussianWHAnnulusO4Omega18p4FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHAnnulusO4Omega18p4Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  domain='annulus'; 
  annulusDims = [.5,1]; % inner and outer radii

  myTitle = sprintf('CgWaveHoltz: Gaussian, Annulus, \\omega=%.4g, L_{2h}-res',omega);

elseif( strcmp(caseOption,'gaussianDiskO4Omega8p1')==1 )

  % -- DISK : Gaussian source----
  plotName = 'gaussianDiskO4Omega8p1';

  gaussianWHDiskO4Omega8p1FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHDiskO4Omega8p1Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  domain='disk'; 
  annulusDims = [.5,1]; % inner and outer radii

  myTitle = sprintf('CgWaveHoltz: Gaussian, Disk, \\omega=%.4g, L_{2h}-res',omega);


elseif( strcmp(caseOption,'gaussianBoxO4Omega8p1')==1 )

  % -- DISK : Gaussian source----
  plotName = 'gaussianBoxO4Omega8p1';

  gaussianWHBoxO4Omega8p1FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHBoxO4Omega8p1Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  nd=3; 
  domain='square'; 
  boxDims = [1,1,1]; % box dimensions

  myTitle = sprintf('CgWaveHoltz: Gaussian, Box3d, \\omega=%.4g, L_{2h}-res',omega);


elseif( strcmp(caseOption,'gaussianSphereO2Omega8p1')==1 )

  % -- DISK : Gaussian source----
  plotName = 'gaussianSphereO2Omega8p1';

  gaussianWHSphereO2Omega8p1FPI
  itv1 = itv; 
  res1 = res;
  cr1 = convergenceRate; 
  
  gaussianWHSphereO2Omega8p1Krylov
  itv2 = itv;
  res2 = res;
  cr2 = convergenceRate; 

  plotConvergenceRate=1; % we can eval the true CR
  nd=3; 
  domain='sphere'; 
  sphereDims = [1,1,1]; % sphere dimensions

  myTitle = sprintf('CgWaveHoltz: Gaussian, Sphere, \\omega=%.4g, L_{2h}-res',omega);





else
   fprintf('ERROR: unknown caseOption=%s\n',caseOption);
   error(1);

end 

if( plotConvergenceRate )
  % Get theoretical FPI convergence rate 
  omegav=omega;
  Tv = numPeriods*2*pi./omegav;
  par.savePlots=1;
  par.plotName = sprintf('%s',plotName); 

  if( strcmp(domain,'square') )
    % Eigenvalues for a square domain: *check me*
    %  lambda^2 = c^2*( (mx*pi)^2 + (my*pi)^2 )
    nx=32; ny=32; nz=32; 
    nfx=nx/2; nfy=ny/2;  % num frequencies 
    if( nd==2 ) nfz=1; else nfz=nz/2; end 

    numLambda = nfx*nfy; 
    lambdav = zeros( numLambda,1 );
    k=1; 
    for mx=1:nfx
      for my=1:nfy
        for mz=1:nfz

          lambdav(k) = c*sqrt( ((mx-1)*pi/boxDims(1) )^2 + ((my-1)*pi/boxDims(2))^2 + ((mz-1)*pi/boxDims(3))^2 ); k=k+1; 
        end
      end
    end
    lambdav = sort(lambdav);

  elseif( strcmp(domain,'annulus') )
    % fprintf('\n ***** FIX ME: we need the eigenvalues of an annulus ****\n\n');

    % File written by cg/ad/codes/annulusEigenvalues.maple 

    annulusEigenvaluesDirichlet

    % WARNING: Only valid to about numBesselOrder + 5 or so (increase )
    if( omega>numBesselOrder )
      fprintf('\n *** WARNING: Annulus eigenvalues may not all be there for omega=%g\n ***\n\n');
    end


    lambdav = sort(lambdav);

  elseif( strcmp(domain,'disk') )

    % fprintf('\n ***** FIX ME: we need the eigenvalues of a DISK ****\n\n');

    % File written by cg/ad/codes/annulusEigenvalues.maple 

    diskEigenvaluesDirichlet

    % WARNING: Only valid to about numBesselOrder + 5 or so (increase )
    if( omega>numBesselOrder )
      fprintf('\n *** WARNING: Annulus eigenvalues may not all be there for omega=%g\n ***\n\n');
    end


    lambdav = sort(lambdav)

  elseif( strcmp(domain,'sphere') )

    fprintf('\n ***** FIX ME: we need the eigenvalues of a SPHERE ****\n\n');

    % File written by cg/ad/codes/annulusEigenvalues.maple 

    diskEigenvaluesDirichlet

    % WARNING: Only valid to about numBesselOrder + 5 or so (increase )
    if( omega>numBesselOrder )
      fprintf('\n *** WARNING: Annulus eigenvalues may not all be there for omega=%g\n ***\n\n');
    end


    lambdav = sort(lambdav)    

  else
    fprintf('ERROR: unknown domain=%s\n',domain);
    error(1);
  end;


  rate = getWaveHoltzConvergenceRate( omegav,Tv, lambdav,par );
  fprintf('WaveHoltz fix-point convergence rate=%9.3e (theory)\n',rate);
end 


figure(1); 
semilogy( itv1,res1,'r-x', itv2,res2,'b-o'); hold on;
legendNames{1}=sprintf('FP         CR=%.3g',cr1);
legendNames{2}=sprintf('GMRES CR=%.3g',cr2);

if( plotConvergenceRate )

  numFixPointIterations=length(itv1); 
  scale = sum(res1)/numFixPointIterations;
  x = [1; numFixPointIterations];
  y = [1; rate^(numFixPointIterations-1)]*scale;
  semilogy( x,y,'k--'); 
  legendNames{3} = sprintf('FP (theory)=%.3g',rate); 

end;

% annotation('textbox',[.5 .9 .1 .2],'String','Text outside the axes','EdgeColor','none');
x=.2; y=.6; w=.4; h=.3; 
myText = sprintf('Order=%d',orderOfAccuracy);
annotation('textbox',[x y w h],'String',myText,'FitBoxToText','on','FontSize',fontSize);

grid on;
title(myTitle);
legend(legendNames);
xlabel('iterations'); 
hold off;



%   % make sure figure size is correct
pos = get(gcf,'position');  
pos(3) = xwidth;
pos(4) = ywidth;
set(gcf,'position',pos);

savePlotFile(plotName,'pdf'); 


% %% generate plots

% % ----------- Major Case ----------
% n=1; 
% for icase = 1:nCases   

%   fprintf('------- option=%s, icase = %d, n=%d ----------\n',caseOption,icase,n);
  
%   [h1,ev1,pv1,qv1]=getData(n,Cases,N0); n=n+1; 
%   if( numSubCases>1 ) [h2,ev2,pv2,qv2]=getData(n,Cases,N0); n=n+1; end
%   if( numSubCases>2 ) [h3,ev3,pv3,qv3]=getData(n,Cases,N0); n=n+1; end
%   if( numSubCases>3 ) [h4,ev4,pv4,qv4]=getData(n,Cases,N0); n=n+1; end

%   % figure
%   figure
  

%   % number of refinements
%   Ng = length(ev1);

%   %% --- get reference line3 ----
%   order2= 2; % expected order
%   [expected2] = getReferenceLine( ev1,pv1,order2,Ng );
%   order4 = 4; % expected order
%   [expected4] = getReferenceLine( ev2,pv2,order4,Ng );


%   purple   =[.6 .4 1];
%   medBlue  =[.3 .3 1];
%   lightBlue=[.8 .8 1];

%   medRed   =[1 .3 .3 ];
%   lightRed =[1 .8 .8 ];

%   %% plot
%   if( 1==1 )
%     h = ...
%     loglog(h1,ev1,'-ob','linewidth',lw,'markersize',ms,'MarkerFaceColor',medBlue  ,'MarkerEdgeColor','k');
%     hold on;
%     loglog(h1,pv1,'-sb','linewidth',lw,'markersize',ms,'MarkerFaceColor',medBlue  ,'MarkerEdgeColor','k');
%     loglog(h1,qv1,'-sb','linewidth',lw,'markersize',ms,'MarkerFaceColor',medBlue  ,'MarkerEdgeColor','k');
%     if( numSubCases>1 )
%       loglog(h2,ev2,'-dr','linewidth',lw,'markersize',ms,'MarkerFaceColor',medRed,'MarkerEdgeColor','k');
%       loglog(h2,pv2,'-^r','linewidth',lw,'markersize',ms,'MarkerFaceColor',medRed,'MarkerEdgeColor','k');
%       loglog(h2,qv2,'-^r','linewidth',lw,'markersize',ms,'MarkerFaceColor',medRed,'MarkerEdgeColor','k');
%     end
%     % loglog(h3,ev3,'-or','linewidth',lw,'markersize',ms,'MarkerFaceColor',medRed   ,'MarkerEdgeColor','k');
%     % loglog(h3,pv3,'-sr','linewidth',lw,'markersize',ms,'MarkerFaceColor',medRed   ,'MarkerEdgeColor','k');
%     % loglog(h4,ev4,'-dr','linewidth',lw,'markersize',ms,'MarkerFaceColor',lightRed ,'MarkerEdgeColor','k');
%     % loglog(h4,pv4,'-^r','linewidth',lw,'markersize',ms,'MarkerFaceColor',lightRed ,'MarkerEdgeColor','k');

%     % Plot rate lines 
%     loglog(h1,expected2,'k--',...
%            h2,expected4,'k-','linewidth',lw,'markersize',ms);
%     % loglog(h1,expected2,'k--',...
%     %        h3,expected4,'k-','linewidth',lw,'markersize',ms);
  
%   end; 


%   %% labels
%   title(Titles{icase},'interpreter','latex');
%   ylabel('Max error','interpreter','latex');
%   xlabel('$h$','interpreter','latex');

%   %% create legend
%   leg = legend('$\mathbf{E}$, O2','$\mathbf{P}$, O2','$N$, O2', ...
%                '$\mathbf{E}$, O4','$\mathbf{P}$, O4','$N$, O4', ...
%                '$h^2$ ref','$h^4$ ref', ...
%                'Location','SouthEast');
%   % leg = legend('$\mathbf{E}$, N O2','$\mathbf{P}$, N O2', ...
%   %              '$\mathbf{E}$, R O2','$\mathbf{P}$, R O2', ...
%   %              '$\mathbf{E}$, N O4','$\mathbf{P}$, N O4', ...
%   %              '$\mathbf{E}$, R O4','$\mathbf{P}$, R O4', ...
%   %              '$h^2$ ref','$h^4$ ref', ...
%   %              'Location','SouthEast');

%   set(leg,'interpreter','latex');

%   %% create xticks
%   %% set(gca,'XLim',[0.125,1]/N0(n));

%   xlim([min(h1),max(h1)]);

%   NN = round(1./h1);
%   xticks = {};
%   for i = 1:Ng
%     xticks{i} = sprintf('$1/%d$',NN(i));
%   end
%   set(gca,'xtick',h1(end:-1:1));
%   set(gca,'xticklabel',xticks(end:-1:1),'ticklabelinterpreter','latex');

%   set(gca,'ytick',10.^(-16:2:2));
  
%   % ylim([1.e-10,.1]);
  
%   %% grid
%   grid on;
%   set(gca,'xminorgrid','off');   set(gca,'xminortick','off');   set(gca,'yminorgrid','off');   set(gca,'yminortick','off');

%   set(gca,'fontsize',fs);

%   drawnow;
%   % make sure figure size is correct
%   pos = get(gcf,'position');  
%   pos(3) = xwidth;
%   pos(4) = ywidth;
%   set(gcf,'position',pos);

%   savePlotFile(plotName{icase},'pdf'); 

%    % name = sprintf('%s.eps',plotName{icase});
%    % print('-depsc',name);
%    % fprintf('Wrote file [%s]\n',name); 
%    % system(sprintf('/Users/henshaw/bin/pspdf %s',name),'-echo');


  
  
end


%-----------------------------------------------------------------------
% Return the errors from the conv.p matlab file
%-----------------------------------------------------------------------
function [h,ev,pv,qv]=getData( n,Cases,N0 )

  if length(Cases{n}) >= 63
    % my matlab doesn't run really long file names for some reason...
    newFile = [Cases{n}(1:60),'.m']; 
    copyfile([Cases{n},'.m'],newFile);
    run(newFile);
  else
    run(Cases{n});
  end
  
  h  = hh/N0(n); % grid spacing starts at 1/10 for this case
  
  
  % 2D : u0=Ex, u1=Ey, u2=Hz, u3=Ev u4=Pv u5=Qv
  % 3D : u0=Ex, u1=Ey, u2=Hz, u3=Ev u4=Pv u5=Qv   (check)
  ev = u3;    % | Ev |
  pv = u4;    % | Pv |
  qv = u5;    % | N |
  
end

%-----------------------------------------------------------------------
%% --- reference line ----
% place reference line inside largest gap between lines
%-----------------------------------------------------------------------
function [expected] = getReferenceLine( ev1,pv1,order,Ng )
  lastErrors = [ev1(end),pv1(end)];
  if ev1(end) < 1e-14
    lastErrors = [pv1(end)];
  end
  lastErrorsSorted = sort(lastErrors);
  [~,I] = max(diff(log(lastErrorsSorted)));
  avgErr = exp(mean(log(lastErrorsSorted(I:(I+1)))));

  expected = avgErr*(2^order).^((Ng-1):(-1):0);
end

