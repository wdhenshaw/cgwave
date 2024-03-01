// This file automatically generated from plot.bC with bpp.
//#define BOUNDS_CHECK
//#define OV_DEBUG

#include "CgWave.h"
#include "PlotStuff.h"
#include "GL_GraphicsInterface.h"
#include "DialogData.h"
// #include "CompositeGridOperators.h"
#include "ParallelUtility.h"
#include "display.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

static int numberOfPushButtons=0, numberOfTextBoxes=0;



// =============================================================================================
/// \brief Create a grid function that holds all the things we can plot.
/// \param current (input) : use this grid function
/// \param v (output) : save the augmented solution here.
/// \param t (input) : if t<0 then only fill the component names into the grid function v.
// =============================================================================================
realCompositeGridFunction& CgWave::
getAugmentedSolution(int current, realCompositeGridFunction & v, const real t)
{

    realCompositeGridFunction *& ucg = dbase.get<realCompositeGridFunction*>("ucg");

  // if( true ) return ucg[current]; // test 

    realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
  
    const int numberOfDimensions = cg.numberOfDimensions();
    
    int numberOfComponents=1;

    const Real & c                          = dbase.get<real>("c");
    const int & addForcing                  = dbase.get<int>("addForcing");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const aString & knownSolutionOption     = dbase.get<aString>("knownSolutionOption");

    const int & solveForScatteredField      = dbase.get<int>("solveForScatteredField");
    const int & plotScatteredField          = dbase.get<int>("plotScatteredField");

  // ** FIX ME: 
    int & computeErrors = dbase.get<int>("computeErrors");
    bool twilightZone = addForcing && forcingOption==twilightZoneForcing;
  // computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";

    int plotErrors= computeErrors;
  // plotErrors=0;

    const bool saveErrors = plotErrors;

    Range all;
    int numberToPlot=numberOfComponents;                  // save fields
    int nErr=numberToPlot;    numberToPlot += numberOfComponents*int(saveErrors);

    const int nScat =numberToPlot;
    if( plotScatteredField ) 
        numberToPlot += numberOfComponents;

  // we build a grid function with more components (errors, dissipation) for plotting
    v.updateToMatchGrid(cg,all,all,all,numberToPlot);

    v=0;
    v.setName("u",0);
    if( plotErrors )
        v.setName("err",1);

    if( plotScatteredField )
    {
        if( solveForScatteredField )
            v.setName("total",2); // plot total field 
        else
            v.setName("scat",2);
    }

    if( t<0. )
    {
    // in this case we only assign the component names and return 
        return v;
    }

    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];

        realMappedGridFunction & u = ucg[current][grid];
        realMappedGridFunction & vg = v[grid];
        OV_GET_SERIAL_ARRAY(real,u,uLocal);
        OV_GET_SERIAL_ARRAY(real,vg,vLocal);
        
        Range N=numberOfComponents;
        vLocal(all,all,all,N)=uLocal(all,all,all,N); // for now make a copy *** fix this **

        if( saveErrors )
        {
            OV_GET_SERIAL_ARRAY(real,error[grid],errLocal);
            vLocal(all,all,all,nErr)=errLocal(all,all,all);
        }

        if( plotScatteredField )
        {
      // -- plot the SCATTERED FIELD (if solveForScatteredField=0)
      // OR
      // -- plot the TOTAL FIELD (if solveForScatteredField=1)

      // Parameters for the plane wave defining the incident field
      //   sin( kx*x + ky*y +kz*z - omega*t + phi )
            const Real amp   = dbase.get<Real>("ampPlaneWave");
            const Real kx    = dbase.get<Real>("kxPlaneWave");
            const Real ky    = dbase.get<Real>("kyPlaneWave");
            const Real kz    = dbase.get<Real>("kzPlaneWave");
            const Real phi   = dbase.get<Real>("phiPlaneWave");
            const Real omega = dbase.get<Real>("omegaPlaneWave");

      // Add or subtract the plane wave from the solution
            const Real fieldSign = solveForScatteredField ? 1. : -1.;
      // printF("plotScatteredField: Using amp=%g, [kx,ky,kz]=[%g,%g,%g]*2*pi, phi=%g, fieldSign=%g\n",amp,kx/twoPi,ky/twoPi,kz/twoPi,phi,fieldSign);

            mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
            OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
            getIndex(mg.dimension(),I1,I2,I3);

            Real omegat = omega*t - phi; 
            if( numberOfDimensions==2 )
                vLocal(I1,I2,I3,nScat) = uLocal(I1,I2,I3,0) + (fieldSign)*amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) - omegat );
            else
                vLocal(I1,I2,I3,nScat) = uLocal(I1,I2,I3,0) + (fieldSign)*amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + kz*xLocal(I1,I2,I3,2) - omegat );

        }
    }

    if( plotScatteredField )
        { // interpolate scattered field 
        v.interpolate();
    }

    const int & useSuperGrid            = dbase.get<int>("useSuperGrid");
    const int & adjustPlotsForSuperGrid = dbase.get<int>("adjustPlotsForSuperGrid");    // set solution to zero in any superGridLayers
    if( adjustPlotsForSuperGrid && useSuperGrid )
    {
        printF("getAugmentedSolution: adjustPlotsForSuperGrid...\n");

        adjustSolutionForSuperGrid( v );
    }

    return v;
    
}



// =============================================================================================
/// \brief Build the run time dialog for CgWave
// =============================================================================================
int CgWave::
buildRunTimeDialog()
{
    GenericGraphicsInterface & ps = gi;
    GUIState *& runTimeDialog = dbase.get<GUIState*>("runTimeDialog");
    
    if( runTimeDialog==NULL )
    {
        runTimeDialog = new GUIState;       
        GUIState & dialog = *runTimeDialog;
        

        dialog.setWindowTitle("CgWave");
        dialog.setExitCommand("finish", "finish");

        aString cmds[] = {"break",
                                            "continue",
                                            "movie mode",
                                            "movie and save",
                                            "contour", 
                                            "grid",
                                              "erase",
                                            "plot options...",
                      //  "parameters...",
                                            "plot v",
                                            ""};

        numberOfPushButtons=9;  // number of entries in cmds
        int numRows=(numberOfPushButtons+1)/2;
        dialog.setPushButtons( cmds, cmds, numRows ); 

    // get any extra components such as errors for tz flow or the pressure for CNS.
        realCompositeGridFunction v;
        real t=-1; // this means only fill in the component names. 
        realCompositeGridFunction & u = getAugmentedSolution(0,v,t);

        const int numberOfComponents = u.getComponentBound(0)-u.getComponentBase(0)+1;
    // create a new menu with options for choosing a component.
        aString *cmd = new aString[numberOfComponents+1];
        aString *label = new aString[numberOfComponents+1];
        for( int n=0; n<numberOfComponents; n++ )
        {
            label[n]=u.getName(n);
            cmd[n]="plot:"+u.getName(n);

        }
        cmd[numberOfComponents]="";
        label[numberOfComponents]="";
        
        dialog.addOptionMenu("plot component:", cmd,label,0);
        delete [] cmd;
        delete [] label;

        const int numberOfTextStrings=8;
        aString textLabels[numberOfTextStrings];
        aString textStrings[numberOfTextStrings];

        const real & tFinal       = dbase.get<real>("tFinal");
        const real & tPlot        = dbase.get<real>("tPlot");
        const real & cfl          = dbase.get<real>("cfl");
        const int & debug         = dbase.get<int>("debug");
        const int & plotFrequency = dbase.get<int>("plotFrequency");

        int nt=0;
        textLabels[nt] = "final time";  sPrintF(textStrings[nt], "%g",tFinal);  nt++; 
        textLabels[nt] = "times to plot";  sPrintF(textStrings[nt], "%g",tPlot);  nt++; 
        textLabels[nt] = "cfl";  
        sPrintF(textStrings[nt], "%g",cfl);  nt++; 
        textLabels[nt] = "debug";  sPrintF(textStrings[nt], "%i",debug);  nt++; 
        textLabels[nt] = "plotFrequency";  sPrintF(textStrings[nt], "%i",plotFrequency);  nt++; 
  
       // null strings terminal list
        textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
        dialog.setTextBoxes(textLabels, textLabels, textStrings);
        numberOfTextBoxes=nt;
        

    }
    return 0;

}


static void
setSensitivity( GUIState & dialog, bool trueOrFalse )
{
    dialog.getOptionMenu(0).setSensitive(trueOrFalse);
    int n;
    for( n=1; n<numberOfPushButtons; n++ ) // leave first push button sensitive (=="break")
        dialog.setSensitive(trueOrFalse,DialogData::pushButtonWidget,n);
    
    for( n=0; n<numberOfTextBoxes; n++ )
        dialog.setSensitive(trueOrFalse,DialogData::textBoxWidget,n);
    
}


// ========================================================================================
/// \brief Return the method name: FD22, FD24s etc.
// ========================================================================================
aString CgWave::getMethodName() const
{
    const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
    const int & upwind                = dbase.get<int>("upwind");

    aString suffix="";
    if( upwind )
        suffix = suffix + "u";  // name = FD22u, FD44u etc. for sosup upwind dissipation

    aString methodName = sPrintF("FD%i%i%s",orderOfAccuracyInTime,orderOfAccuracy,(const char *)suffix);

    return methodName;
}

// ========================================================================================
/// \brief Construct the label that defines the time-stepping method and used in the
///    title for plots
// ========================================================================================
void CgWave::
getTimeSteppingLabel( real dt, aString & label ) const
{
    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

    aString buff;
    label=sPrintF(buff,"dt=%4.1e",dt);

    if( timeSteppingMethod==explicitTimeStepping )
        label+=" TS=ME";
    else
        label+=" TS=IM";


  // if( method==nfdtd  )
  // {
  //   if( artificialDissipation!=0. && artificialDissipation==artificialDissipationCurvilinear )
  //     label+=sPrintF(buff," ad%i=%4.2f",orderOfArtificialDissipation,artificialDissipation);
  //   else if( artificialDissipationCurvilinear!=0. )
  //     label+=sPrintF(buff," adr%i=%4.2f,adc%i=%4.2f",orderOfArtificialDissipation,artificialDissipation,
  //                    orderOfArtificialDissipation,artificialDissipationCurvilinear);
            
  //   if( applyFilter )
  //     label+=sPrintF(buff,", filter%i",orderOfFilter);

  //   if( divergenceDamping>0. )
  //     label+=sPrintF(buff," dd=%5.3f",divergenceDamping);
  // }
  // if( method==bamx )
  // {
  //   // BAMX currently just uses a filter stage
  //   const int filterOrder=orderOfAccuracyInSpace + 2;
        
  //   if( artificialDissipation!=0. )
  //     label+=sPrintF(buff," F%i=%3.1f",filterOrder,artificialDissipation);

  //   if( divergenceDamping>0. )
  //     label+=sPrintF(buff," dd=%5.3f",divergenceDamping);
  // }
    
}



int CgWave::
plot( int current, real t, real dt )
// ========================================================================================
/// /brief Plot the solution.
/// \param current (input) : index of solution to plot.
/// \param t (input) : time of solution to plot
/// \param dt (input) : current time-step
/// 
/// 
///  plotOptions :  0 = no plotting
///                 1 - plot and wait
///                 3 - do not wait for response after plotting
/// /Return values: 0=normal exit. 1=user has requested "finish".
// ========================================================================================
{
    real timep=getCPU();
    
    int & plotOptions = dbase.get<int>("plotOptions");
    int & plotChoices = dbase.get<int>("plotChoices");

    GenericGraphicsInterface & ps = gi;
    PlotStuffParameters & psp =  dbase.get<PlotStuffParameters>("psp");

  // we need to know if the graphics is open on any processor -- fix this in the GraphicsInterface.
    int graphicsIsOn = ps.isGraphicsWindowOpen();
    graphicsIsOn=ParallelUtility::getMaxValue(graphicsIsOn);
    int readingCommandFile = ps.readingFromCommandFile();
    readingCommandFile=ParallelUtility::getMaxValue(readingCommandFile);

    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

    if( plotOptions==noPlotting || !graphicsIsOn )  // ************ WDH***  March 24, 2023 **TEMP
        return 0;

  // printF(" $$$$$$$$ cgWave: plot: plotOptions=%d t=%9.3e\n",plotOptions,t);

  // real cpu0=getCPU();
    int returnValue=0;


    const int & myid                  = dbase.get<int>("myid");
    const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");

    const int orderOfAccuracyInSpace=orderOfAccuracy;


    int & movieFrame        = dbase.get<int>("movieFrame");
    aString & movieFileName = dbase.get<aString>("movieFileName");

    real & tFinal           = dbase.get<real>("tFinal");
    real & tPlot            = dbase.get<real>("tPlot");
    real & cfl              = dbase.get<real>("cfl");
    int & debug             = dbase.get<int>("debug");
    int & plotFrequency     = dbase.get<int>("plotFrequency");

    aString methodName      = getMethodName();  // FD22, or FD24s etc.

    char buff[100];
    psp.set(GI_TOP_LABEL,sPrintF(buff,"CgWave %s t=%6.2e ",(const char *)methodName,t));
    aString label;
    getTimeSteppingLabel( dt,label );
    
    psp.set(GI_TOP_LABEL_SUB_1,label);



  // DialogData &plotOptionsDialog = *pPlotOptionsDialog;
  // DialogData &parametersDialog = *pParametersDialog;

    assert( dbase.get<GUIState*>("runTimeDialog")!=NULL );
    GUIState *& runTimeDialog = dbase.get<GUIState*>("runTimeDialog");
    GUIState & dialog = *runTimeDialog;

    aString answer;

  // Build the augmented solution (includes errors)
    realCompositeGridFunction v;
    realCompositeGridFunction & u = getAugmentedSolution(current,v,t);  // u is either solution or v

    const int numberOfComponents = u.getComponentBound(0)-u.getComponentBase(0)+1;

    if( movieFrame>=0   )
    { // save a ppm file as part of a movie.
        
        psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::ppm);
        ps.outputString(sPrintF(buff,"Saving file %s%i.ppm",(const char*)movieFileName,movieFrame));
        ps.hardCopy(    sPrintF(buff,            "%s%i.ppm",(const char*)movieFileName,movieFrame),psp);
        psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::postScript);
        movieFrame++;
    }

    

    ps.erase();
    if( plotOptions != noPlotting )
    {

    // Plot all the the things that the user has previously plotted
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);

        if( plotChoices & 1 )
        {
            PlotIt::plot(ps,cg,psp);
        }
        if( plotChoices & 2 )
            PlotIt::contour(ps,u,psp);
        if( plotChoices & 4 )
            PlotIt::streamLines(ps,u,psp);


        bool programHalted=false;
        int checkForBreak=false;
        const int processorForGraphics = ps.getProcessorForGraphics();

        if( plotOptions==plotNoWait && !(ps.readingFromCommandFile()) &&
                myid==processorForGraphics && ps.isGraphicsWindowOpen() )
        { // we are running interactively and we should check for a "break" command:
            checkForBreak=true; 
        }
        broadCast(checkForBreak,processorForGraphics); // broadcast to all from the processor for graphics

        if( checkForBreak )
        {
      // movie mode ** check here if the user has hit break ***
      // ps.outputString(sPrintF(buff,"Check for break at t=%e\n",t));
            answer="";
            
            int menuItem = ps.getAnswerNoBlock(answer,"monitor>");
      // printf("answer = [%s]\n",(const char*)answer);
            
            if( answer=="break" )
            {
                programHalted=true;
            }
            if( t> tFinal-dt*.5 )
                programHalted=true;
        }

    // printF("\n $$$$$$$$ plot: plotOptions=%d, programHalted=%d, t=%9.3e\n\n",plotOptions,programHalted,t);

        if( plotOptions==plotAndWait || programHalted )
        {
      // -------------- Display Run Time Dialog and wait for answers ----------

            if( plotOptions==plotNoWait )
            {
                setSensitivity( dialog,true );
            }

            if( programHalted )
            {
                plotOptions=plotAndWait; // reset movie mode if set.
                movieFrame=-1;
            }
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);

//       DialogData & fileOutputDialog = dialog.getDialogSibling(1);
//       DialogData & pdeDialog = dialog.getDialogSibling(2);

            int len;
            bool replot=false;
            for(;;)
            {
        // int menuItem = ps.getMenuItem(menu,answer,"choose an option");
                real timew=getCPU();
                int menuItem = ps.getAnswer(answer,"");
                timing(timeForWaiting)+=getCPU()-timew;

                if( answer=="contour" )
                {
                    if(plotChoices & 2 )
                        ps.erase();
                    
                    PlotIt::contour(ps,u,psp);
                    if( psp.getObjectWasPlotted() ) 
                        plotChoices |= 2;
                    else
                        plotChoices &= ~2;

                    replot=true;
                }

                else if( answer=="grid" )
                {
                    PlotIt::plot(ps,cg,psp);

                    if( psp.getObjectWasPlotted() ) 
                        plotChoices |= 1;
                    else
                        plotChoices &= ~1;
                }
                else if( answer=="erase" )
                {
                    ps.erase();
                    plotChoices=0;
                }

                
                else if( answer=="continue" )
                {
                    if( t >= tFinal-dt/10. )
                    {
                        printF("WARNING: t=tFinal. Choose `finish' if you really want to end\n");
                        printF(" INFO: t=%9.3e, dt=%9.3e, tFinal=%9.3e\n",t,dt,tFinal);
                    }
                    else
                        break;
                }
                else if( answer=="movie mode" )
                {

                    if( ps.isGraphicsWindowOpen() )
                        plotOptions=plotNoWait;   // don't wait
                    else
                        plotOptions=noPlotting; 

                    setSensitivity( dialog,false );
                    break;
                }
                else if( answer=="movie and save" )
                {
                    ps.inputString(answer,"Enter basic name for the ppm files (default=plot)");
                    if( answer !="" && answer!=" ")
                        movieFileName=answer;
                    else
                        movieFileName="plot";
                    ps.outputString(sPrintF(buff,"pictures will be named %s0.ppm, %s1.ppm, ...",
                        (const char*)movieFileName,(const char*)movieFileName));

                    movieFrame=0;

                    if( ps.isGraphicsWindowOpen() )
                        plotOptions=plotNoWait;   // don't wait
                    else
                        plotOptions=noPlotting; 

                    setSensitivity( dialog,false );
                    break;
                }
                else if( answer=="finish" )
                {
                    printF("plot: finish chosen...\n");
                    tFinal=t;
                    returnValue=1;
                    break;
                }

                else if( dialog.getTextValue(answer,"cfl","%g",cfl) ){}//
                else if( dialog.getTextValue(answer,"final time","%g",tFinal) ){}//
                else if( dialog.getTextValue(answer,"times to plot","%g",tPlot) ){}//
                else if( dialog.getTextValue(answer,"debug","%i",debug) ){}//
                else if( dialog.getTextValue(answer,"plotFrequency","%i",plotFrequency) ){}//

                else if( len=answer.matches("plot:") )
                {
          // plot a new component
                    aString name = answer(len,answer.length()-1);
                    int component=-1;
                    for( int n=0; n<numberOfComponents; n++ )
                    {
                        if( name==u.getName(n) )
                        {
                            component=n;
                            break;
                        }
                    }
                    if( component==-1 )
                    {
                        printF("ERROR: unknown component name =[%s]\n",(const char*)name);
                        component=0;
                    }

                    ps.erase();

                    dialog.getOptionMenu(0).setCurrentChoice(component);
                    if( plotChoices & 2 )
                    {
                        psp.set(GI_COMPONENT_FOR_CONTOURS,component);
                        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);

                        if( plotChoices & 1 )
                        {
                            PlotIt::plot(ps,cg,psp);
                        }
                        PlotIt::contour(ps,u,psp);

                        if( plotChoices & 4 )
                            PlotIt::streamLines(ps,u,psp);

                        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
                    }
                }
                else if( answer=="plot v" )
                {
          // plot the Helmhotz integral v  = time integral of u 
                    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
                    v.setName("Helmholtz v",0);
                    
                    ps.erase();
                    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                    PlotIt::contour(ps,v,psp);

                }

            
                else if( answer=="break" )
                {
                }
                else
                {
                    printF("plot: Unknown response=[%s]\n",(const char*)answer);
                }
                if( replot )
                {
                    replot=false;
                    getAugmentedSolution(current,v,t);
                    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
                    ps.erase();


                    if( plotChoices & 1 )
                    {
                        PlotIt::plot(ps,cg,psp);
                    }
                    if( plotChoices & 2 )
                        PlotIt::contour(ps,u,psp);
                    if( plotChoices & 4 )
                        PlotIt::streamLines(ps,u,psp);
    
                    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                }//
            }
        }
    }
    
    
    if( plotOptions==plotNoWait  )
    {
        ps.redraw(TRUE);
        
    }
    
    timing(timeForPlotting)+=getCPU()-timep;
    return returnValue;
}
