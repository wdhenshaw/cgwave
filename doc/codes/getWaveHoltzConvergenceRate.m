%
% Return the asymptotic convergence rate for the Wave-Holtz iteration
%
% lamv(i) : eigenvalues of c^2\Delta
% par.savePlots
%
% 
function rate = getMultiFrequencyConvergenceRate( omegav,Tv, lambdav,par )

  % General form of 
  %       beta(lambda; omega,T) = (2/T) int0^T (cos(omega*t)-.25)*cos(lambda*t) dt 
  % No assumptions on omega or T (i.e. we do NOT assume that T=2*pi/omega))
  betaGeneral = @(lambda,omega,T) ( mySinc((omega-lambda)*T)  + mySinc((omega+lambda)*T)- .5*mySinc(lambda*T) );  

  numLambda = length(lambdav);
  % lambdav
  
  mud = zeros(numLambda,1); % discrete mu's 
  
  for n=1:numLambda
    lam = lambdav(n);

    mud(n) = betaGeneral(lam,omegav,Tv);

  end
  % mud(1:8)
  % mud(1:100)
  % max(abs(mud))

  % convergence rate: 
  rate = max( max(abs(mud)) );

  if( 1==1 )

    % ---- plot mu curve and discrete mu corresponding to discrete lambda ----
    N=500;
    lambdaMin=0; 
    % lambdaMax=max(lambdav); 
    lambdaMax=100; % fix me 
    % plot a range of lambda
    numd=numLambda; 
    for( i=1:numLambda )
      if( lambdav(i)>lambdaMax )
        numd = i; break;
      end
    end
    L=1:numd; % plot these lambdav

    lamv = linspace(lambdaMin,lambdaMax,N); 
    muv = zeros(N,1); % discrete mu's 
  
    for n=1:N
      lam = lamv(n);
      muv(n) = betaGeneral(lam,omegav,Tv);
    
    end  

    figure(3)
    plot(lamv,abs(muv(:)),'-', lambdav(L),abs(mud(L)),'x' ); hold on;
    plot(lambdav(L),0*L,'k+' ); hold on;
    
    % plot vertical lines at  omega
    plot( [omegav;omegav], [0,1.25], 'k-','LineWidth',1);

    title(sprintf('WaveHoltz \\mu, \\omega=%g, rate=%.5g',omegav,rate));
    legend('\mu','\mu(\lambda_j)');
    grid on; xlabel('\lambda'); ylabel('|\mu|');  
    hold off; 

    % Hardcopy 
    if( par.savePlots )
      savePlotFile( sprintf('%sFixPointMuFunction',par.plotName),'pdf' );
    end       

  end
  
  return
