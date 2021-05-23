%
%  Save a matlab plot as a pdf 
%
% plotName (input) :  plot-name (without suffix '.pdf')
% type (input) : 'pdf' (default)
%
function savePlotFile( plotName, type )

  if( nargin<2 ) type ='pdf'; end 

  if( strcmp(type,'pdf') )
    epsPlotName=sprintf('%s.eps',plotName); 

    %% print('-depsc2',epsPlotName); % save as a .eps file (temporary)
    % new way: Aug, 2020
    saveas(gcf,epsPlotName,'epsc2'); 

    system(sprintf('/usr/local/bin/ps2pdf -dEPSCrop %s.eps %s.pdf',plotName,plotName),'-echo');
    system(sprintf('rm %s',epsPlotName)); % now remove the .eps 
    fprintf('savePlotFile:saved plot [%s.pdf].\n',plotName);
  else
    fprintf('savePlotFile:ERROR: unknown file type=[%s]\n',type);
  end 

end 
