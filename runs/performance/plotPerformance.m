
%
% FINISH ME -- from ogmg/doc/tex/plotPerformance.m
%

%
% Plot performance comparisons
%
% Examples:
%   plotPerformance -case=solve -verticalBars=0 -savePlots=0
%   plotPerformance -case=solve -verticalBars=0 -relative=1 -savePlots=0
%   plotPerformance -case=solve -verticalBars=0 -cpuSpeed=2.7 -savePlots=0  % Cg6 is 2.7GHtz processor
%
%   plotPerformance -case=ARC -verticalBars=0 -relative=1 -savePlots=0   % time for advance rect/curv only
%   plotPerformance -case=ARC -verticalBars=0 -cpuSpeed=2.7 -savePlots=0   % time for advance rect/curv only
% 
%  relative=1 : relative to order2 results
%  -cpuSpeed=GHZ  (if set, plot cycles/step/pt )

function plotPerformance(varargin)

caseOption = 'solve';
savePlots=0;
verticalBars=0; % vertical or horizontal bars
yMax=0; yMin=1;
relative=0;     % show relative reults, normalized by Ogmg
cpuSpeed = -1;       % give processor speed in GHz

% --- read command line args ---
for i = 1 : nargin
  line = varargin{i};
  caseOption       = getString( line,'-caseOption',caseOption );
  caseOption       = getString( line,'-case',caseOption );

  relative         = getInt( line,'-relative',relative );
  savePlots        = getInt( line,'-savePlots',savePlots );
  verticalBars     = getInt( line,'-verticalBars',verticalBars );
  yMax             = getReal( line,'-yMax',yMax );
  yMin             = getReal( line,'-yMin',yMin );
  cpuSpeed         = getReal( line,'-cpuSpeed',cpuSpeed );
end 



% Set defaults for plotting 
fontSize=20; lineWidth=2; markerSize=8; 
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);
set(0,'DefaultLegendFontSize',fontSize);





% --- LOAD PERFORMANCE DATA -----

perfData;

% data


% relative CPU times
% data(cpuTotal,  2  ,1,square1024)=    0.58  ; % data(:,order,solver,gridType)= value

numSchemes=1;

numGrids=2; 
gridName{1}="square"; 
gridName{2}="nonSquare";    
% gridName{3}="shapes"; 

% quant=cpuTotal; 
if( strcmp(caseOption,'solve')==1 )
  quant=cpuSolve; quantName{1}="solve"; 
  if relative==1  
    myTitle=sprintf('Relative solve CPU');    
    axisLabel='CPU (relative to O2)';    
    plotName='relativeSolveTimes'; 
  else
    if cpuSpeed>0
      myTitle=sprintf('Solve Cycles');
      axisLabel='CPU cycles/step/pt';  
      plotName='cyclesSolveTimes'; 
    else  
      myTitle=sprintf('Solve CPU');
      axisLabel='CPU ns/step/pt';  
      plotName='solveTimes'; 
    end 
  end 
  % gMax(1)=180; 
  % gMax(2)=140; 
  % gMax(3)=55;

elseif( strcmp(caseOption,'ARC')==1 )
  quant=cpuARC; quantName{1}="advance"; 
  if relative==1  
    myTitle=sprintf('Relative interior update');    
    axisLabel='CPU (relative to O2)';    
    plotName='relativeAdvanceTimes'; 
  else
    if cpuSpeed>0
      myTitle=sprintf('Interior update cycles');
      axisLabel='CPU cycles/step/pt';  
      plotName='cyclesAdvanceTimes'; 
    else  
      myTitle=sprintf('Interior update CPU');
      axisLabel='CPU ns/step/pt';  
      plotName='advanceTimes'; 
    end 
  end 
  % gMax(1)=180; 
  % gMax(2)=140; 
  % gMax(3)=55;
elseif( strcmp(caseOption,'memory')==1 )
  quant=storage; quantName{1}="memory";  
  if relative==1 
    myTitle=sprintf('Relative memory usage'); axisLabel='Memory'; plotName='relativeMemoryUsage';
  else
    myTitle=sprintf('Memory usage'); axisLabel='Memory reals/point'; plotName='memoryUsage';
  end 
  gMax(1)=200; % 25; 
  gMax(2)=200; % 20; 
  gMax(3)=200; % 15;   
else
  fprintf('ERROR: unknown case=%s\n',caseOption); pause; pause;
end

gridType=1; % square
scheme=1;
order2=2; 
if cpuSpeed<0 clockSpeed=1; else clockSpeed=1/cpuSpeed; end

for group=1:3  % order of accuracy in each bar-group 
  ord = 2*group;
  if relative==1 
    cpu(group,1:numGrids) = data(quant,ord,scheme,1:numGrids)./data(quant,order2,scheme,1:numGrids);
  else
    cpu(group,1:numGrids) = data(quant,ord,scheme,1:numGrids)/clockSpeed;
  end
end

cpu

if verticalBars 
  hb = bar(cpu(1:3,1:numGrids));
else
  hb = barh(cpu(1:3,1:numGrids));
end
title(sprintf('%s, %s',myTitle,solverName{1}));
str = {'O2'; 'O4'; 'O6'};
if verticalBars 
  set(gca, 'XTickLabel',str, 'XTick',1:numel(str));
  ylabel(axisLabel);
else
  set(gca, 'YTickLabel',str, 'YTick',1:numel(str)); ytickangle(90);
  xlabel(axisLabel);
end
legend('square','nonSquare','Location','best');
% legend('Ogmg','AMG','BiCGSt','GMRES','Location','best');
grid on; 

  hb(1).FaceColor = [.2 .6 .5];
  hb(2).FaceColor = [.4 .4 1.];
%  hb(2).FaceColor = [1 0 0];
%   hb(3).FaceColor = [0 0 1];
%   hb(4).FaceColor = [.3 .2 .7]; % purple

% -------------------- LABEL TOP OF BARS -------------
yMax=max(max(cpu(1:3,1:2)));
relCpu = cpu(1:3,1:2);
xShift0=.3;    % shift in direction of bars
if cpuSpeed>0 
  yShift0=2;
else
  yShift0=yMax/40;
end 
labelSize=20;

labelBars( hb,relCpu,yMax,xShift0,yShift0,labelSize,verticalBars );


  if savePlots 
    if verticalBars 
      fullPlotName=sprintf('%s%sVerticalBars',plotName,gridName{gridType});
    else
      fullPlotName=sprintf('%sHorizontalBars%s',plotName,gridName{gridType});
    end
    savePlotFile(fullPlotName,'pdf');   
  end

return
end


% ----------------------------------------------------
% Function to put numbers on top of bar chart bars
% ----------------------------------------------------
function labelBars( hb,relCpu,yMax,xShift0,yShift0,labelSize,verticalBars )

  nd=2; 

  M=relCpu;
  numbersToAdd=relCpu;
  barWidth = hb.BarWidth;
  numCol = size(M,2);
  cnt = 0;
  for ii = numbersToAdd'
      cnt = cnt + 1;
      xPos = linspace(cnt - barWidth/2, cnt + barWidth / 2, numCol+1);
      idx = 1;
      group=0;
      for jj = xPos(1:end-1)
          group=group+1; 
          val = numbersToAdd(cnt,idx);
          y = min(yMax-0,M(cnt,idx));
          % xShift=(.08 - (idx-1)*.02)*barWidth; 
          xShift=(xShift0 - (idx-1)*.085)*barWidth; 
          % xShift=(.275 - (idx-1)*.085)*barWidth; 
          yShift=yShift0; 
          % fprintf('idx=%g, jj=%g xShift=%g num=%s\n',idx,jj,xShift,num2str(val,'%0.1f'));
          if verticalBars
            if 1==0 && (val==1 || group==1 )
              ht=text(jj+xShift, y + yShift, sprintf('%0.1f (best)',val),'fontSize',labelSize); set(ht,'Rotation',90);
            else
              ht=text(jj+xShift, y + yShift, num2str(val,'%0.1f'),'fontSize',labelSize); set(ht,'Rotation',90);
            end
          else
            % horizontal bars
            if 1==0 && (val==1 || group==1)
              ht=text(y + yShift, jj+xShift, sprintf('%0.1f (best)',val), 'fontSize',labelSize); set(ht,'Rotation',0);
            else
              ht=text(y + yShift, jj+xShift,  num2str(val,'%0.1f'),'fontSize',labelSize); set(ht,'Rotation',0);
            end            
          end
          idx = idx +1;
      end     
  end
  return
end
