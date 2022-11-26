#  
# Generate code for the modified equation time-stepping
#    VERSION 2:
#        HIERACHICAL SCHEME BASED on LAPLACIAN
#
# restart: currentdir("/Users/henshaw/Dropbox/research/cgwave"): currentdir("/Users/henshaw/Dropbox/research/cgwave/maple"): read "writeModifiedEquationCode.mpl";
#
# with(LinearAlgebra):
# with(plottools):
# with(plots):
with(StringTools):
# with(CodeGeneration):

# kernelopts(printbytes=false): # turn off memory used messages
# kernelopts('printlevel'=0): 
# interface(echo=0): 
# interface(quiet=true): 
interface( quiet=true,echo=0 ):

rectangular:=1: curvilinear:=2: # enumerators 
gridTypeName[rectangular]:="Rectangular"; 
gridTypeName[curvilinear]:="Curvilinear"; 

ndStart := 2: ndEnd:=3:     # number of dimensions
orderStart:=6: orderEnd:=6: # order of accuracy
# orderOfAccuracy  := 8:      # order of accuracy
# gridType         := rectangular: 
# gridType         := curvilinear;
gridTypeStart:=1: gridTypeEnd:=2: # gridType 

# --- difference weights from highOrderDiff.maple -----
# Weights in deriv 1
if 1=1 then
  dw[0,1]:=1; dw[1,1]:=-1/6; dw[2,1]:=1/30; dw[3,1]:=-1/140; dw[4,1]:=1/630; dw[5,1]:=-1/2772; 
  # Weights in deriv 2
  dw[0,2]:=1; dw[1,2]:=-1/12; dw[2,2]:=1/90; dw[3,2]:=-1/560; dw[4,2]:=1/3150; dw[5,2]:=-1/16632; 
  # Weights in deriv 3
  dw[0,3]:=1; dw[1,3]:=-1/4; dw[2,3]:=7/120; dw[3,3]:=-41/3024; dw[4,3]:=479/151200; dw[5,3]:=-59/79200; 
  # Weights in deriv 4
  dw[0,4]:=1; dw[1,4]:=-1/6; dw[2,4]:=7/240; dw[3,4]:=-41/7560; dw[4,4]:=479/453600; dw[5,4]:=-59/277200; 
  # Weights in deriv 5
  dw[0,5]:=1; dw[1,5]:=-1/3; dw[2,5]:=13/144; dw[3,5]:=-139/6048; dw[4,5]:=37/6480; dw[5,5]:=-4201/2993760; 
  # Weights in deriv 6
  dw[0,6]:=1; dw[1,6]:=-1/4; dw[2,6]:=13/240; dw[3,6]:=-139/12096; dw[4,6]:=37/15120; dw[5,6]:=-4201/7983360; 
end if:

# printf("Read definitions of chain rule derivatives...\n");
# read "chainRuleCoefficientsDefinitions2d.maple":
# read "chainRuleCoefficientsDefinitions3d.maple":
# printf("...done\n");

  # ---------- START PROCEDURE declareVar --------
  #
  # Auxillary function to declare a variable
  #
  declareVar :=proc( name );
    global declareCount, dfile;
    if name = finishDeclarations then
      fprintf(dfile,"\n");
    else
      # printf("declareVar: declareCount=%s\n",convert(declareCount,string)):
      if declareCount=0 then
        fprintf(dfile,"real %s",name);
      else
        fprintf(dfile,", %s",name);
      end if:
      declareCount:=declareCount+1;
      if declareCount > numDeclaredPerLine then
        fprintf(dfile,"\n");
        declareCount:=0;
      end if:
    end if:
  end:
  # ---------- END PROCEDURE declareVar --------


  # ---------- START PROCEDURE saveDpDm --------
  #
  # Auxillary function to print a second UNDIVIDED difference formulae
  # Example:
  #    saveDpDm( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveDpDm :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c,lhsDeclare;

    if source=rx or source=rsxy then
      c := "m1,m2"; # component name 
    else
      c := "0": 
    end if:
    d1m:=d1; d2m:=d2; d3m:=d3; 
    if d1>0  then
      # give preferences for differences in the first direction (better memory accesses)
      d1m := d1-2; 
    elif d2>0 then
      d2m:=d2-2;
    elif d3>0 then 
      d3m:=d3-2;   
    else
      printf("saveDpDm: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;    

    lhs := sprintf("%s%d%d%d(i1,i2,i3,%s)",convert(target,string),d1,d2,d3,c); # LHS name

    if shortTargetName=1 then
      # short LHS name, use "i" instead of  (i1,i2,i3,0)
      if source=rx or source=rsxy then
        lhs := sprintf("%s%d%d%di(m1,m2)",convert(target,string),d1,d2,d3); 
      else
        lhs := sprintf("%s%d%d%di",convert(target,string),d1,d2,d3);
      end if:       
      # lhs := sprintf("%s%d%d%di",convert(target,string),d1,d2,d3);  
    end if;
    if target=source and (d1m+d2m+d3m)>0 then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    lhsDeclare := SubstituteAll(lhs,"(i1,i2,i3,0)","(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
    declareVar(lhsDeclare): # ----- DECLARE THE VARIABLE

    # fprintf(file,"    d400(i1,i2,i3,0) = (d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0))/(dx(0)**2)\n",
    
    if d1>0  then
      # give preferences for differences in the first direction (better memory accesses)
      # d1m := d1-2; 
      fprintf(file,"    %s = %s(i1+1,i2,i3,%s) - 2*%s(i1,i2,i3,%s) + %s(i1-1,i2,i3,%s)\n",lhs,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1+1,i2,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1-1,i2,i3,%s))/(dx(0)**2)\n",lhs,rhs,c,rhs,c,rhs,c );
    elif d2>0 then
      # d2m:=d2-2;
      fprintf(file,"    %s = %s(i1,i2+1,i3,%s) - 2*%s(i1,i2,i3,%s) + %s(i1,i2-1,i3,%s)\n",lhs,rhs,c,rhs,c,rhs,c ); 
      # fprintf(file,"    %s = (%s(i1,i2+1,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1,i2-1,i3,%s))/(dx(1)**2)\n",lhs,rhs,c,rhs,c,rhs,c ); 
    elif d3>0 then 
      # d3m:=d3-2;   
      fprintf(file,"    %s = %s(i1,i2,i3+1,%s) - 2*%s(i1,i2,i3,%s) + %s(i1,i2,i3-1,%s)\n",lhs,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1,i2,i3+1,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1,i2,i3-1,%s))/(dx(2)**2)\n",lhs,rhs,c,rhs,c,rhs,c );
    else
      printf("saveDpDm: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;
    
    RETURN(lhs):
  end:
  # ---------- END PROCEDURE --------

  # ---------- START PROCEDURE saveD0 --------
  #
  # Auxillary function to print a FIRST UNDIVIDED difference formulae
  # Example:
  #    saveD0( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveD0 :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c,lhsDeclare;

    if source=rx or source=rsxy then
      c := "m1,m2"; # component name 
    else
      c := "0": 
    end if:

    d1m:=d1; d2m:=d2; d3m:=d3; # for rhs, one less derivative 
    if d1>0 and (d1 mod 2 =1) then # RHS must use even derivatives 
      # give preferences for differences in the first direction (better memory accesses)
      d1m := d1-1; # for rhs, one less derivative 
    elif d2>0 and (d2 mod 2 =1) then
      d2m:=d2-1;
    elif d3>0 and (d3 mod 2 =1) then 
      d3m:=d3-1;   
    else
      printf("saveD0: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;    

    lhs := sprintf("%s%d%d%d(i1,i2,i3,%s)",convert(target,string),d1,d2,d3,c); # LHS name

    if shortTargetName=1 then
      # short LHS name, use "i" instead of  (i1,i2,i3,0)
      if source=rx or source=rsxy then
        lhs := sprintf("%s%d%d%di(m1,m2)",convert(target,string),d1,d2,d3); 
      else
        lhs := sprintf("%s%d%d%di",convert(target,string),d1,d2,d3);
      end if:  
    end if;
    if target=source  and (d1m+d2m+d3m)>0 then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    lhsDeclare := SubstituteAll(lhs,"(i1,i2,i3,0)","(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
    declareVar(lhsDeclare): # ----- DECLARE THE VARIABLE    

    if d1>0 and (d1 mod 2 =1) then
      # give preferences for differences in the first direction (better memory accesses)
      fprintf(file,"    %s = %s(i1+1,i2,i3,%s) - %s(i1-1,i2,i3,%s)\n",lhs,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1+1,i2,i3,%s) - %s(i1-1,i2,i3,%s))/(2.*dx(0))\n",lhs,rhs,c,rhs,c );
    elif d2>0 and (d2 mod 2 =1) then
      fprintf(file,"    %s = %s(i1,i2+1,i3,%s) - %s(i1,i2-1,i3,%s)\n",lhs,rhs,c,rhs,c ); 
      # fprintf(file,"    %s = (%s(i1,i2+1,i3,%s) - %s(i1,i2-1,i3,%s))/(2.*dx(1))\n",lhs,rhs,c,rhs,c ); 
    elif d3>0 and (d3 mod 2 =1) then 
      fprintf(file,"    %s = %s(i1,i2,i3+1,%s) - %s(i1,i2,i3-1,%s)\n",lhs,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1,i2,i3+1,%s) - %s(i1,i2,i3-1,%s))/(2.*dx(2))\n",lhs,rhs,c,rhs,c );
    else
      printf("saveDpDm: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;
    
    RETURN(lhs):
  end:
  # ---------- END PROCEDURE save D0 --------

  # ---------- START PROCEDURE saveD0D0 --------
  #
  # Auxillary function to print a MIXED UNDIVIDED DERIVATIVE difference formulae
  # Example:
  #    saveD0D0( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveD0D0 :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c,lhsDeclare;

    if source=rx or source=rsxy then
      c := "m1,m2"; # component name 
    else
      c := "0": 
    end if:
    d1m:=d1; d2m:=d2; d3m:=d3; # for rhs, one less derivative 
    if d1>0 and (d1 mod 2 =1) then # RHS must use even derivatives 
      # give preferences for differences in the first direction (better memory accesses)
      d1m := d1-1; # for rhs, one less derivative 
    end if:
    if d2>0 and (d2 mod 2 =1) then
      d2m:=d2-1;
    end if:
    if d3>0 and (d3 mod 2 =1) then 
      d3m:=d3-1;   
    end if;    

    lhs := sprintf("%s%d%d%d(i1,i2,i3,%s)",convert(target,string),d1,d2,d3,c); # LHS name

    if shortTargetName=1 then
      # short LHS name, use "i" instead of  (i1,i2,i3,0)
      if source=rx or source=rsxy then
        lhs := sprintf("%s%d%d%di(m1,m2)",convert(target,string),d1,d2,d3); 
      else
        lhs := sprintf("%s%d%d%di",convert(target,string),d1,d2,d3);
      end if:        
    end if;
    if target=source and (d1m+d2m+d3m)>0 then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    lhsDeclare := SubstituteAll(lhs,"(i1,i2,i3,0)","(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
    declareVar(lhsDeclare): # ----- DECLARE THE VARIABLE    

    if d1>0 and (d1 mod 2 =1) and d2>0 and (d2 mod 2 = 1) then
      # give preferences for differences in the first direction (better memory accesses)
      fprintf(file,"    %s = %s(i1+1,i2+1,i3,%s) - %s(i1-1,i2+1,i3,%s) - %s(i1+1,i2-1,i3,%s) + %s(i1-1,i2-1,i3,%s)\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1+1,i2+1,i3,%s) - %s(i1-1,i2+1,i3,%s) - %s(i1+1,i2-1,i3,%s) + %s(i1-1,i2-1,i3,%s))/(4.*dx(0)*dx(1))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
    elif d1>0 and (d1 mod 2 =1) and d3>0 and (d3 mod 2 = 1) then
      fprintf(file,"    %s = %s(i1+1,i2,i3+1,%s) - %s(i1-1,i2,i3+1,%s) - %s(i1+1,i2,i3-1,%s) + %s(i1-1,i2,i3-1,%s)\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1+1,i2,i3+1,%s) - %s(i1-1,i2,i3+1,%s) - %s(i1+1,i2,i3-1,%s) + %s(i1-1,i2,i3-1,%s))/(4.*dx(0)*dx(2))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
    elif d2>0 and (d2 mod 2 =1) and d3>0 and (d3 mod 2 = 1) then
      fprintf(file,"    %s = %s(i1,i2+1,i3+1,%s) - %s(i1,i2-1,i3+1,%s) - %s(i1,i2+1,i3-1,%s) + %s(i1,i2-1,i3-1,%s)\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1,i2+1,i3+1,%s) - %s(i1,i2-1,i3+1,%s) - %s(i1,i2+1,i3-1,%s) + %s(i1,i2-1,i3-1,%s))/(4.*dx(1)*dx(2))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
    else
      printf("saveD0D0: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;
    
    RETURN(lhs):
  end:
  # ---------- END PROCEDURE save D0D0 --------

  # ---------- START PROCEDURE saveD0D0D0 --------
  #
  # Auxillary function to print a third MIXED DERIVATIVE difference formulae
  # Example:
  #    saveD0D0D0( d,d,1,1,1 )
  #    saveD0D0D0( d,d,3,1,1 )
  #    saveD0D0D0( d,d,3,1,3 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveD0D0D0 :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c;

    if source=rx or source=rsxy then
      c := "m1,m2"; # component name 
    else
      c := "0": 
    end if:
    d1m:=d1; d2m:=d2; d3m:=d3; # for rhs, one less derivative 
    if d1>0 and (d1 mod 2 =1) then # RHS must use even derivatives 
      d1m := d1-1; # for rhs, one less derivative 
    end if:
    if d2>0 and (d2 mod 2 =1) then
      d2m:=d2-1;
    end if:
    if d3>0 and (d3 mod 2 =1) then 
      d3m:=d3-1;   
    end if;    

    lhs := sprintf("%s%d%d%d(i1,i2,i3,%s)",convert(target,string),d1,d2,d3,c); # LHS name

    if shortTargetName=1 then
      # short LHS name, use "i" instead of  (i1,i2,i3,0)
      if source=rx or source=rsxy then
        lhs := sprintf("%s%d%d%di(m1,m2)",convert(target,string),d1,d2,d3); 
      else
        lhs := sprintf("%s%d%d%di",convert(target,string),d1,d2,d3);
      end if:        
    end if;
    if target=source then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    # D0xD0yD0z
    fprintf(file,"    %s = ( (%s(i1+1,i2+1,i3+1,%s) - %s(i1-1,i2+1,i3+1,%s) - %s(i1+1,i2-1,i3+1,%s) + %s(i1-1,i2-1,i3+1,%s)) ",lhs,rhs,c,rhs,c,rhs,c,rhs,c):
    fprintf(file,"          -(%s(i1+1,i2+1,i3-1,%s) - %s(i1-1,i2+1,i3-1,%s) - %s(i1+1,i2-1,i3-1,%s) + %s(i1-1,i2-1,i3-1,%s))  )/(8.*dx(0)*dx(1)*dx(2))\n",rhs,c,rhs,c,rhs,c,rhs,c ):


    RETURN(lhs):
  end:
  # ---------- END PROCEDURE saveD0D0D0 --------


# ---------- START PROCEDURE startLoops --------
#   startLoops( file ):
#
  startLoops :=proc( file );
   fprintf(file,"beginLoops3d()\n");
   fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):  
   fprintf(file,"  if( mask(i1,i2,i3).ne.0 )then\n");
   fprintf(file,"  #End \n"):
  end:
# ---------- END PROCEDURE startLoops --------

# ---------- START PROCEDURE endLoops --------
#   endLoops( file ):
#
  endLoops :=proc( file );
  fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
  fprintf(file,"    end if ! mask .ne. 0\n");
  fprintf(file,"  #End \n"):  
  fprintf(file,"  endLoops3d() \n");  
 end:
# ---------- END PROCEDURE endLoops --------

# ---------------------------------------------------------------------------
# ------------ START LOOP OVER GRIDTYPE, ORDER AND DIMENSIONS ---------------
# ---------------------------------------------------------------------------
for gridType from gridTypeStart to gridTypeEnd do 
for orderOfAccuracy from orderStart by 2 to orderEnd do 
for nd from ndStart to ndEnd do 

  # printf("--- START: orderOfAccuracy=%d, nd=%d\n",orderOfAccuracy,nd):

 # At eighth order we split loops to improve efficiency 
 if orderOfAccuracy=8 then splitLoops:=1: else splitLoops:=0: end if:

 declareCount:=0;         # global variable for counting declared variables
 numDeclaredPerLine:=10;  # write this many declared variables per line  

numGhost := orderOfAccuracy/2; 
rNames[0]:="r": rNames[1]:="s": rNames[2]:="t":
xNames[0]:="x": xNames[1]:="y": xNames[2]:="z":

# -- define file names: 
myFileName := sprintf("update%ddOrder%d%s.h",nd,orderOfAccuracy,gridTypeName[gridType]); 
declareFileName := sprintf("declare%ddOrder%d%s.h",nd,orderOfAccuracy,gridTypeName[gridType]);  # for declarations of variables 


file := open(myFileName, WRITE);
dfile := open(declareFileName, WRITE);  # for declarations of variables 

  # if gridType=curvilinear then
  #   fprintf(file,"! -----------------------------------------------------------------------\n"):
  #   fprintf(file,"! Macro to assign neighbours of metric terms when comuputing derivatives:\n"):
  #   fprintf(file,"!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)\n"):
  #   fprintf(file,"!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)\n"):
  #   fprintf(file,"!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)\n"):
  #   fprintf(file,"! Extrapolate when necessary\n"):
  #   fprintf(file,"! -----------------------------------------------------------------------\n"):
  #   fprintf(file,"#beginMacro getMetricNeighbours%ddOrder%d(m1,m2)\n",nd,orderOfAccuracy):
  #   fprintf(file,"  do ig=2,%d\n",numGhost):
  #   fprintf(file,"    if( i1-ig.ge.nd1a )then\n"):
  #   fprintf(file,"      rxi1g(-ig) = rsxy(i1-ig,i2,i3,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi1g(-ig) = extrapJac(i1-ig,i2,i3,m1,m2,1,0,0)\n");
  #   fprintf(file,"    end if\n"):
  #   fprintf(file,"    if( i1+ig.le.nd1b )then\n"):
  #   fprintf(file,"      rxi1g(+ig) = rsxy(i1+ig,i2,i3,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi1g(+ig) = extrapJac(i1+ig,i2,i3,m1,m2,-1,0,0)\n");
  #   fprintf(file,"    end if\n"):
  # 
  #   fprintf(file,"    if( i2-ig.ge.nd2a )then\n"):
  #   fprintf(file,"      rxi2g(-ig) = rsxy(i1,i2-ig,i3,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi2g(-ig) = extrapJac(i1,i2-ig,i3,m1,m2,0,1,0)\n");
  #   fprintf(file,"    end if\n"):
  #   fprintf(file,"    if( i2+ig.le.nd2b )then\n"):
  #   fprintf(file,"      rxi2g(+ig) = rsxy(i1,i2+ig,i3,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi2g(+ig) = extrapJac(i1,i2+ig,i3,m1,m2,0,-1,0)\n");
  #   fprintf(file,"    end if\n"):        
  #   if nd=3 then
  #   fprintf(file,"    if( i3-ig.ge.nd3a )then\n"):
  #   fprintf(file,"      rxi3g(-ig) = rsxy(i1,i2,i3-ig,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi3g(-ig) = extrapJac(i1,i2,i3-ig,m1,m2,0,0,1)\n");
  #   fprintf(file,"    end if\n"):
  #   fprintf(file,"    if( i3+ig.le.nd3b )then\n"):
  #   fprintf(file,"      rxi3g(+ig) = rsxy(i1,i2,i3+ig,m1,m2)\n");
  #   fprintf(file,"    else\n"):
  #   fprintf(file,"      rxi3g(+ig) = extrapJac(i1,i2,i3,m1,m2,0,0,-1)\n");
  #   fprintf(file,"    end if\n"):        
  #   end if: 
  #   fprintf(file,"    end do\n"):
  #   fprintf(file,"#endMacro ! end getMetricNeighbours(m1,m2)\n"):  
  #   fprintf(file,"\n"): 
  # end if:   

fprintf(file,"! ===========================================================================\n");
fprintf(file,"!   Modified Equation : order=%d, DIMENSIONS=%d, gridType=%s\n",orderOfAccuracy,nd,gridTypeName[gridType]);
fprintf(file,"!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl\n");
fprintf(file,"!  MASK : USEMASK or NOMASK \n"):
fprintf(file,"!  FORCING : NOFORCING, TZ, USEFORCING \n"):
fprintf(file,"! ===========================================================================\n");
fprintf(file,"#beginMacro update%ddOrder%d%s(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)\n",nd,orderOfAccuracy,gridTypeName[gridType]); 

fprintf(dfile,"! ===========================================================================\n");
fprintf(dfile,"!   Modified Equation : order=%d, DIMENSIONS=%d, gridType=%s\n",orderOfAccuracy,nd,gridTypeName[gridType]);
fprintf(dfile,"!  Macro to DECLARE VARIABLES\n");
fprintf(dfile,"!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl\n");
fprintf(dfile,"! ===========================================================================\n");
fprintf(dfile,"#beginMacro declare%ddOrder%d%s()\n",nd,orderOfAccuracy,gridTypeName[gridType]); 

fprintf(file,"#If #FORCING eq noForcing\n"):
fprintf(file,"#defineMacro FV(m) \n"):
fprintf(file,"#Else\n"):
fprintf(file,"#defineMacro FV(m) +dtSq*fv(m)\n"):
fprintf(file,"#End\n"):


fprintf(file,"\n! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----\n");
fprintf(file,"! Example: \n");
if gridType=rectangular then
  fprintf(file,"! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...\n");
else
  fprintf(file,"! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...\n");
end if:
maxDeriv := orderOfAccuracy:
numTerms := orderOfAccuracy/2-1:

cx[0]:="": cx[1]:="": cx[2]:="":
for d1 from 1 by 1 to maxDeriv-2 do
  if gridType=rectangular then
    cx[0] := cat(cx[0],"x");   cx[1] := cat(cx[1],"y");  cx[2] := cat(cx[2],"z");
  else
    cx[0] := cat(cx[0],"r");   cx[1] := cat(cx[1],"s");  cx[2] := cat(cx[2],"t");
  end if:
  for ax from 0 to nd-1 do
    fprintf(file,"c%s0 = 1.; ",cx[ax]);
    declareVar(sprintf("c%s0",cx[ax])):
    # fprintf(dfile,"real c%s0\n",cx[ax]);
    for i from 1 to numTerms do
      fprintf(file,"c%s%d = %s.; ",cx[ax],i,convert(dw[i,d1],string),ax,2*i);
      # fprintf(file,"c%s%d = (%s.)*dx(%d)**%d; ",cx[ax],i,convert(dw[i,d1],string),ax,2*i);
      declareVar(sprintf("c%s%d",cx[ax],i)):
      # fprintf(dfile,"real c%s%d\n",cx[ax],i);
    end do;
    fprintf(file,"\n"):
  end do:
end do:


if gridType=curvilinear then
  fprintf(file,"dr1=dr(0); dr1i=1./dr1;\n"); 
  fprintf(file,"dr2=dr(1); dr2i=1./dr2;\n"); 
  fprintf(file,"dr3=dr(2); dr3i=1./dr3;\n"); 
  declareVar("dr1");  declareVar("dr2"):  declareVar("dr3"):
  declareVar("dr1i"); declareVar("dr2i"): declareVar("dr3i"):
else
  fprintf(file,"cxx=1./dx(0)**2;\n");
  fprintf(file,"cyy=1./dx(1)**2;\n");
  fprintf(file,"czz=1./dx(2)**2;\n");
  declareVar("cxx");  declareVar("cyy"):  declareVar("czz"):
end if:


fprintf(file,"fv(m)=0.\n");

if gridType=curvilinear then
  fprintf(file,"\n"):
  if nd=2 then
    fprintf(file,"#defineMacro c200(i1,i2,i3) lapCoeff(i1,i2,i3,0)\n"):
    fprintf(file,"#defineMacro c020(i1,i2,i3) lapCoeff(i1,i2,i3,1)\n"):
    fprintf(file,"#defineMacro c110(i1,i2,i3) lapCoeff(i1,i2,i3,2)\n"):
    fprintf(file,"#defineMacro c100(i1,i2,i3) lapCoeff(i1,i2,i3,3)\n"):
    fprintf(file,"#defineMacro c010(i1,i2,i3) lapCoeff(i1,i2,i3,4)\n"):
  else
    fprintf(file,"#defineMacro c200(i1,i2,i3) lapCoeff(i1,i2,i3,0)\n"):
    fprintf(file,"#defineMacro c020(i1,i2,i3) lapCoeff(i1,i2,i3,1)\n"):
    fprintf(file,"#defineMacro c002(i1,i2,i3) lapCoeff(i1,i2,i3,2)\n"):
    fprintf(file,"#defineMacro c110(i1,i2,i3) lapCoeff(i1,i2,i3,3)\n"):
    fprintf(file,"#defineMacro c101(i1,i2,i3) lapCoeff(i1,i2,i3,4)\n"):
    fprintf(file,"#defineMacro c011(i1,i2,i3) lapCoeff(i1,i2,i3,5)\n"):
    fprintf(file,"#defineMacro c100(i1,i2,i3) lapCoeff(i1,i2,i3,6)\n"):
    fprintf(file,"#defineMacro c010(i1,i2,i3) lapCoeff(i1,i2,i3,7)\n"):   
    fprintf(file,"#defineMacro c001(i1,i2,i3) lapCoeff(i1,i2,i3,8)\n"):
  
  end if:
  #  fprintf(file,"\n! we may need to extrapolate some metrics to extra ghost points \n"):
  #  if orderOfAccuracy=2 then
  #    fprintf(file,"#defineMacro extrapJac(i1,i2,i3,m1,m2,is1,is2,is3) (3.*rsxy(i1+(is1),i2+(is2),i3+(is3),m1,m2)- 3.*rsxy(i1+2*(is1),i2+2*(is2),i3+2*(is3),m1,m2) +    rsxy(i1+3*(is1),i2+3*(is2),i3+3*(is3),m1,m2))\n"):
  #  elif orderOfAccuracy=4 then
  #    fprintf(file,"#defineMacro extrapJac(i1,i2,i3,m1,m2,is1,is2,is3) (4.*rsxy(i1+(is1),i2+(is2),i3+(is3),m1,m2)- 6.*rsxy(i1+2*(is1),i2+2*(is2),i3+2*(is3),m1,m2) + 4.*rsxy(i1+3*(is1),i2+3*(is2),i3+3*(is3),m1,m2)-    rsxy(i1+4*(is1),i2+4*(is2),i3+4*(is3),m1,m2))\n"):
  #  elif orderOfAccuracy=6 then
  #    fprintf(file,"#defineMacro extrapJac(i1,i2,i3,m1,m2,is1,is2,is3) (5.*rsxy(i1+(is1),i2+(is2),i3+(is3),m1,m2)-10.*rsxy(i1+2*(is1),i2+2*(is2),i3+2*(is3),m1,m2) +10.*rsxy(i1+3*(is1),i2+3*(is2),i3+3*(is3),m1,m2)- 5.*rsxy(i1+4*(is1),i2+4*(is2),i3+4*(is3),m1,m2) +   rsxy(i1+5*(is1),i2+5*(is2),i3+5*(is3),m1,m2))\n"):
  #  elif orderOfAccuracy=8 then
  #    fprintf(file,"#defineMacro extrapJac(i1,i2,i3,m1,m2,is1,is2,is3) (6.*rsxy(i1+(is1),i2+(is2),i3+(is3),m1,m2)-15.*rsxy(i1+2*(is1),i2+2*(is2),i3+2*(is3),m1,m2) +20.*rsxy(i1+3*(is1),i2+3*(is2),i3+3*(is3),m1,m2)-15.*rsxy(i1+4*(is1),i2+4*(is2),i3+4*(is3),m1,m2) + 6*rsxy(i1+5*(is1),i2+5*(is2),i3+5*(is3),m1,m2) - rsxy(i1+6*(is1),i2+6*(is2),i3+6*(is3),m1,m2))\n"):
  #  else
  #    error("orderOfAccuracy"):
  #  end if:
end if:


if gridType=curvilinear then

  fprintf(file,"\nif( c200(0,0,0).le.0. )then\n");
  fprintf(file,"\n  ! --- Evaluate and store coefficients in Laplacian ---\n");

  fprintf(file,"  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'\n");

  fprintf(file,"\n"):


  # $$$$$ NEW WAY $$$$$
  fprintf(file,"  numGhost1=%d;\n",numGhost-1);
  fprintf(file,"  n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);\n");
  fprintf(file,"  n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);\n");
  fprintf(file,"  n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);\n");

  startLoops( file ):

  # fprintf(file,"  beginLoops3d()\n");
  # fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):  
  # fprintf(file,"    if( mask(i1,i2,i3).ne.0 )then\n");
  # fprintf(file,"  #End \n"):

  for m1 from 0 to nd-1 do
  for m2 from 0 to nd-1 do
    fprintf(file,"    %s%s = rsxy(i1,i2,i3,%d,%d)\n",rNames[m1],xNames[m2],m1,m2):  # rx = rsxy(i1,i2,i3,0,0)
    declareVar(sprintf("%s%s",rNames[m1],xNames[m2])); 
  end do:
  end do:

  declareVar("diffOrder1"); 
  declareVar("diffOrder2"); 
  declareVar("diffOrder3"); 
  fprintf(file,"\n    ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---\n"):
  for dir from 1 to nd do 
    # printf(">> orderOfAccuracy=%d\n",orderOfAccuracy):

    diffOrder:=orderOfAccuracy: 
    fprintf(file,"    if( (i%d-%d).ge.nd%da .and. (i%d+%d).le.nd%db )then\n",dir,diffOrder/2,dir,dir,diffOrder/2,dir):
    fprintf(file,"      diffOrder%d=%d\n",dir,diffOrder):
    for diffOrder from orderOfAccuracy-2 by -2 to 2 do
      fprintf(file,"    elseif( (i%d-%d).ge.nd%da .and. (i%d+%d).le.nd%db )then\n",dir,diffOrder/2,dir,dir,diffOrder/2,dir):
      fprintf(file,"      diffOrder%d=%d\n",dir,diffOrder):
    end do:
    fprintf(file,"    else\n"):
    fprintf(file,"      stop 999\n"):
    fprintf(file,"    end if\n"):
  end do:

  # compute (r,s,t) derivatives of rx
  fprintf(file,"    if( diffOrder1.eq.2 )then\n"):
  for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do
    fprintf(file,"      %s%sr = (rsxy(i1+1,i2,i3,%d,%d)-rsxy(i1-1,i2,i3,%d,%d))*(.5*dr1i) \n",rNames[m1],xNames[m2],m1,m2,m1,m2):  # rxr = 
  end do: end do: 
  if orderOfAccuracy>=4 then
    fprintf(file,"    elseif( diffOrder1.eq.4 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do
    fprintf(file,"      %s%sr = ( 8*(rsxy(i1+1,i2,i3,%d,%d)-rsxy(i1-1,i2,i3,%d,%d)) -(rsxy(i1+2,i2,i3,%d,%d)-rsxy(i1-2,i2,i3,%d,%d)) )*(dr1i/12.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2): 
    end do: end do:
  end if:
  if orderOfAccuracy>=6 then
    fprintf(file,"    elseif( diffOrder1.eq.6 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%sr = ( 45.*(rsxy(i1+1,i2,i3,%d,%d)-rsxy(i1-1,i2,i3,%d,%d)) -9.*(rsxy(i1+2,i2,i3,%d,%d)-rsxy(i1-2,i2,i3,%d,%d)) +(rsxy(i1+3,i2,i3,%d,%d)-rsxy(i1-3,i2,i3,%d,%d)) )*(dr1i/60.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2):     
    end do: end do:
  end if:
  if orderOfAccuracy>=8 then
    fprintf(file,"    elseif( diffOrder1.eq.8 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%sr = ( 672.*(rsxy(i1+1,i2,i3,%d,%d)-rsxy(i1-1,i2,i3,%d,%d)) -168.*(rsxy(i1+2,i2,i3,%d,%d)-rsxy(i1-2,i2,i3,%d,%d)) +32*(rsxy(i1+3,i2,i3,%d,%d)-rsxy(i1-3,i2,i3,%d,%d)) -3.*(rsxy(i1+4,i2,i3,%d,%d)-rsxy(i1-4,i2,i3,%d,%d)) )*(dr1i/840.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2):     
    end do: end do:
    end if:
  fprintf(file,"    end if\n"):

  fprintf(file,"    if( diffOrder2.eq.2 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%ss = (rsxy(i1,i2+1,i3,%d,%d)-rsxy(i1,i2-1,i3,%d,%d))*(.5*dr2i) \n",rNames[m1],xNames[m2],m1,m2,m1,m2):  # rxs = 
    end do: end do:
  if orderOfAccuracy>=4 then
    fprintf(file,"    elseif( diffOrder2.eq.4 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%ss = ( 8*(rsxy(i1,i2+1,i3,%d,%d)-rsxy(i1,i2-1,i3,%d,%d)) -(rsxy(i1,i2+2,i3,%d,%d)-rsxy(i1,i2-2,i3,%d,%d)) )*(dr2i/12.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2): 
    end do: end do:
  end if:
  if orderOfAccuracy>=6 then
    fprintf(file,"    elseif( diffOrder2.eq.6 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%ss = ( 8*(rsxy(i1,i2+1,i3,%d,%d)-rsxy(i1,i2-1,i3,%d,%d)) -(rsxy(i1,i2+2,i3,%d,%d)-rsxy(i1,i2-2,i3,%d,%d)) )*(dr2i/12.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2): 
    end do: end do:
  end if:
  if orderOfAccuracy>=8 then
    fprintf(file,"    elseif( diffOrder2.eq.8 )then\n"):
    for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
    fprintf(file,"      %s%ss = ( 8*(rsxy(i1,i2+1,i3,%d,%d)-rsxy(i1,i2-1,i3,%d,%d)) -(rsxy(i1,i2+2,i3,%d,%d)-rsxy(i1,i2-2,i3,%d,%d)) )*(dr2i/12.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2): 
    end do: end do:
  end if:
  fprintf(file,"    end if\n"):

  if nd=3 then
    fprintf(file,"    if( diffOrder3.eq.2 )then\n"):
      for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
      fprintf(file,"      %s%st = (rsxy(i1,i2,i3+1,%d,%d)-rsxy(i1,i2,i3-1,%d,%d))*(.5*dr2i) \n",rNames[m1],xNames[m2],m1,m2,m1,m2):  # rxt =
      end do: end do:
    if orderOfAccuracy>=4 then
      fprintf(file,"    elseif( diffOrder3.eq.4 )then\n"):
      for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
      fprintf(file,"      %s%st = ( 8*(rsxy(i1,i2,i3+1,%d,%d)-rsxy(i1,i2,i3-1,%d,%d)) -(rsxy(i1,i2,i3+2,%d,%d)-rsxy(i1,i2,i3-2,%d,%d)) )*(dr3i/12.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2): 
      end do: end do:
    end if:
    if orderOfAccuracy>=6 then
      fprintf(file,"    elseif( diffOrder3.eq.6 )then\n"):
      for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
      fprintf(file,"      %s%st = ( 45.*(rsxy(i1,i2,i3+1,%d,%d)-rsxy(i1,i2,i3-1,%d,%d)) -9.*(rsxy(i1,i2,i3+2,%d,%d)-rsxy(i1,i2,i3-2,%d,%d)) +(rsxy(i1,i2,i3+3,%d,%d)-rsxy(i1,i2,i3-3,%d,%d)))*(dr3i/60.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2): 
      end do: end do:
    end if:
    if orderOfAccuracy>=8 then
      fprintf(file,"    elseif( diffOrder3.eq.8 )then\n"):
      for m1 from 0 to nd-1 do  for m2 from 0 to nd-1 do    
      fprintf(file,"      %s%st = ( 672.*(rsxy(i1,i2,i3+1,%d,%d)-rsxy(i1,i2,i3-1,%d,%d)) -168.*(rsxy(i1,i2,i3+2,%d,%d)-rsxy(i1,i2,i3-2,%d,%d)) +32*(rsxy(i1,i2,i3+3,%d,%d)-rsxy(i1,i2,i3-3,%d,%d)) -3.*(rsxy(i1,i2,i3+4,%d,%d)-rsxy(i1,i2,i3-4,%d,%d)) )*(dr3i/840.) \n",rNames[m1],xNames[m2],m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2,m1,m2): 
      end do: end do:
    end if:
    fprintf(file,"    end if\n"):      
  end if:

  for m1 from 0 to nd-1 do
  for m2 from 0 to nd-1 do
    declareVar(sprintf("%s%sr",rNames[m1],xNames[m2])); 
    declareVar(sprintf("%s%ss",rNames[m1],xNames[m2])); 
    if nd=3 then
      declareVar(sprintf("%s%st",rNames[m1],xNames[m2])); 
    end if: 
  end do:
  end do:

  # write(*,'("i1,i2=",2i3," rx,ry,rxi1g(-2),rxi1g(2),rxi2g(-2),rxi2g(2),ryr,rys=",8(1pe10.2))') i1,i2,rx,ry,rxi1g(-2),rxi1g(2),rxi2g(-2),rxi2g(2),ryr,rys

  # Compute some (x,y,z) derivatives of rx
  for m1 from 0 to nd-1 do
    if nd=2 then
    fprintf(file,"    %sxx = rx*%sxr + sx*%sxs \n",rNames[m1],rNames[m1],rNames[m1]):  # rxx = rx*rxr + sx*rxs
    fprintf(file,"    %syy = ry*%syr + sy*%sys \n",rNames[m1],rNames[m1],rNames[m1]):  # ryy = ry*ryr + sy*rys
    declareVar(sprintf("%sxx",rNames[m1])); 
    declareVar(sprintf("%syy",rNames[m1])); 
    else
    fprintf(file,"    %sxx = rx*%sxr + sx*%sxs + tx*%sxt\n",rNames[m1],rNames[m1],rNames[m1],rNames[m1]):  # rxx = rx*rxr + sx*rxs + tx*rxt
    fprintf(file,"    %syy = ry*%syr + sy*%sys + ty*%syt\n",rNames[m1],rNames[m1],rNames[m1],rNames[m1]):  # ryy = ry*ryr + sy*rys + ty*syt
    fprintf(file,"    %szz = rz*%szr + sz*%szs + tz*%szt\n",rNames[m1],rNames[m1],rNames[m1],rNames[m1]):  # rzz = rz*rzr + sz*rzs + tz*szt
    declareVar(sprintf("%sxx",rNames[m1])); 
    declareVar(sprintf("%syy",rNames[m1])); 
    declareVar(sprintf("%szz",rNames[m1])); 
    end if: 
  end do:    
  fprintf(file,"\n    ! -- Coefficients in the Laplacian (scaled)\n");
  if nd=2 then
  fprintf(file,"    c200(i1,i2,i3) = (rx**2 + ry**2   )*dr1i**2\n");
  fprintf(file,"    c110(i1,i2,i3) = 2.*(rx*sx + ry*sy)*dr1i*dr2i*.25\n");
  fprintf(file,"    c020(i1,i2,i3) = (sx**2 + sy**2   )*dr2i**2\n");
  fprintf(file,"    c100(i1,i2,i3) = (rxx + ryy       )*dr1i*.5\n");
  fprintf(file,"    c010(i1,i2,i3) = (sxx + syy       )*dr2i*.5 \n");  

  else

  fprintf(file,"    c200(i1,i2,i3) = (rx**2 + ry**2 + rz**2 )*dr1i**2\n");
  fprintf(file,"    c020(i1,i2,i3) = (sx**2 + sy**2 + sz**2 )*dr2i**2\n");
  fprintf(file,"    c002(i1,i2,i3) = (tx**2 + ty**2 + tz**2 )*dr3i**2\n");

  fprintf(file,"    c110(i1,i2,i3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25\n");
  fprintf(file,"    c101(i1,i2,i3) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25\n");
  fprintf(file,"    c011(i1,i2,i3) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25\n");

  fprintf(file,"    c100(i1,i2,i3) = (rxx + ryy + rzz)*dr1i*.5\n");
  fprintf(file,"    c010(i1,i2,i3) = (sxx + syy + tyy)*dr2i*.5 \n");          
  fprintf(file,"    c001(i1,i2,i3) = (txx + tyy + tzz)*dr3i*.5 \n"); 


  end if: 
  fprintf(file,"\n");

  endLoops( file ):
  # fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
  # fprintf(file,"    end if ! mask .ne. 0\n");
  # fprintf(file,"  #End \n"):  
  # fprintf(file,"  endLoops3d() \n");

  fprintf(file,"end if ! end assignLapCoeff\n");


  

end if: # end if curvilinear

for ord from 2 by 2 to orderOfAccuracy do
  if nd=2 then ord3:=0; else ord3:=ord; end if;


  if ord=orderOfAccuracy then
    fprintf(file,"\n! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============\n");
  end if:

  shortTargetName:=1:            
  if ord<orderOfAccuracy then 
    arrayNeededForNextLevel:=0:  # this mean array is needed so do NOT used the shortTargetName
    suffix := "(i1,i2,i3,0)":
  else 
    arrayNeededForNextLevel:=1: 
    suffix := "i": 
  end if:     

  fprintf(file,"\n"):
  fprintf(file,"numGhost1=%d;\n",numGhost-ord/2);
  fprintf(file,"n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);\n");
  fprintf(file,"n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);\n");
  fprintf(file,"n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);\n");

  startLoops( file ):
  # fprintf(file,"beginLoops3d()\n");
  # fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):  
  # fprintf(file,"  if( mask(i1,i2,i3).ne.0 )then\n");
  # fprintf(file,"  #End \n"):

  if ord=2 then
    # ----------------- ORDER=2 ---------------- 

    # if ord<orderOfAccuracy then 
    #   shortTargetName:=0: 
    #   suffix := "(i1,i2,i3,0)":
    # else 
    #   shortTargetName:=1: 
    #   suffix := "i": 
    # end if:  

    lhsName:=saveDpDm( d,u,2,0,0,arrayNeededForNextLevel );  # Eval d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
    lhsName:=saveDpDm( d,u,0,2,0,arrayNeededForNextLevel );  # Eval d020(i1,i2,i3,0) = ...
    if nd=3 then
    lhsName:=saveDpDm( d,u,0,0,2,arrayNeededForNextLevel );
    end if;

    if orderOfAccuracy=2 then
      # no need to store these:
      lap2hName := "lap2h": declareVar("lap2h"):
    else
      # need to store these
      lap2hName := "lap2h(i1,i2,i3,0)": declareVar("lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
    end if:

    if gridType=curvilinear then
      # eval: d100i, d010i, d110i 
      # lhsName:=saveD0D0( d,u,d1,d2,d3,shortTargetName );  # rhs uses "u"
      shortTargetName:=1: 
      lhsName:=saveD0( d,u,1,0,0,shortTargetName );    # rhs uses "u"
      lhsName:=saveD0( d,u,0,1,0,shortTargetName );    # rhs uses "u"
      lhsName:=saveD0D0( d,u,1,1,0,shortTargetName );  # rhs uses "u"
      if nd=3 then
      lhsName:=saveD0( d,u,0,0,1,shortTargetName );    # rhs uses "u"
      lhsName:=saveD0D0( d,u,1,0,1,shortTargetName );  # rhs uses "u"
      lhsName:=saveD0D0( d,u,0,1,1,shortTargetName );  # rhs uses "u"
      end if:

      if nd=2 then
        fprintf(file,"    %s = \\\n",lap2hName);
        fprintf(file,"         c200(i1,i2,i3)*d200%s + \\\n",suffix);
        fprintf(file,"         c110(i1,i2,i3)*d110i + \\\n");
        fprintf(file,"         c020(i1,i2,i3)*d020%s +\\\n",suffix);
        fprintf(file,"         c100(i1,i2,i3)*d100i + \\\n");
        fprintf(file,"         c010(i1,i2,i3)*d010i\n");
      else
        fprintf(file,"    %s = \\\n",lap2hName);
        fprintf(file,"         c200(i1,i2,i3)*d200%s +\\\n",suffix);
        fprintf(file,"         c020(i1,i2,i3)*d020%s +\\\n",suffix);
        fprintf(file,"         c002(i1,i2,i3)*d002%s +\\\n",suffix);
        fprintf(file,"         c110(i1,i2,i3)*d110i + \\\n");
        fprintf(file,"         c101(i1,i2,i3)*d101i + \\\n");
        fprintf(file,"         c011(i1,i2,i3)*d011i + \\\n");
        fprintf(file,"         c100(i1,i2,i3)*d100i + \\\n");
        fprintf(file,"         c010(i1,i2,i3)*d010i + \\\n");
        fprintf(file,"         c001(i1,i2,i3)*d001i\n");        
      end if:



    else
      if nd=2 then
        fprintf(file,"    %s = cxx*d200%s + cyy*d020%s \n",lap2hName,suffix,suffix);
      else
        fprintf(file,"    %s = cxx*d200%s + cyy*d020%s + czz*d002%s \n",lap2hName,suffix,suffix,suffix);
      end if:
    end if:

    if ord=2 and orderOfAccuracy=2 then
      fprintf(file,"    #If #FORCING ne \"NOFORCING\" \n"):
      fprintf(file,"      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) \n");
      fprintf(file,"    #End \n\n"):      
      fprintf(file,"\n    ! --- Modified equation space-time update ----\n"):
      fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \\\n");
      fprintf(file,"             + cdtsq*( %s )                         \\\n",lap2hName);
      fprintf(file,"             FV(m)                             \n"):
    end if:


  elif ord=4 then
    # ----------------- ORDER=4 ---------------- 
    # if ord<orderOfAccuracy then 
    #   shortTargetName:=0: 
    #   suffix := "(i1,i2,i3,0)":
    # else 
    #   shortTargetName:=1: 
    #   suffix := "i": 
    # end if:

    lhsName:= saveDpDm( d,d,4,0,0,arrayNeededForNextLevel );  # Eval d400(i1,i2,i3,0)
    lhsName:= saveDpDm( d,d,0,4,0,arrayNeededForNextLevel );  # Eval d040(i1,i2,i3,0)
    if nd=3 then
    lhsName:= saveDpDm( d,d,0,0,4,arrayNeededForNextLevel ); 
    end if:

    


    if orderOfAccuracy=4 then
      # no need to store these:
      lap4hName := "lap4h": declareVar("lap4h"):
      lap2hSqName := "lap2hSq": declareVar("lap2hSq"):
    else
      # need to store these
      lap4hName := "lap4h(i1,i2,i3,0)": declareVar("lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
      lap2hSqName := "lap2hSq(i1,i2,i3,0)": declareVar("lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):        
    end if:

    if gridType=curvilinear then

      if orderOfAccuracy>4 then
        lhsName:= saveDpDm( d,d,2,2,0 );  # d220 is needed for order 6
        if nd=3 then
        lhsName:= saveDpDm( d,d,2,0,2 );  # needed for order 6
        lhsName:= saveDpDm( d,d,0,2,2 );  # needed for order 6
        end if:
      end if:

      lhsName:=saveD0( d,d,3,0,0,shortTargetName );    
      lhsName:=saveD0( d,d,0,3,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( d,d,0,0,3,shortTargetName );
      end if:

      lhsName:=saveD0D0( d,d,3,1,0,shortTargetName );  
      lhsName:=saveD0D0( d,d,1,3,0,shortTargetName );  
      if nd=3 then
      lhsName:=saveD0D0( d,d,3,0,1,shortTargetName );  
      lhsName:=saveD0D0( d,d,1,0,3,shortTargetName );  
      lhsName:=saveD0D0( d,d,0,3,1,shortTargetName );          
      lhsName:=saveD0D0( d,d,0,1,3,shortTargetName );          
      end if:

      fprintf(file,"\n    ! --- Laplacian to order 4 = lap2h + corrections \n");
      if nd=2 then
        fprintf(file,"    %s = lap2h(i1,i2,i3,0) \\\n",lap4hName);
        fprintf(file,"         + c200(i1,i2,i3)*crr1*d400%s \\\n",suffix);
        fprintf(file,"         + c110(i1,i2,i3)*(cr1*d310i + cs1*d130i) \\\n");
        fprintf(file,"         + c020(i1,i2,i3)*css1*d040%s \\\n",suffix);
        fprintf(file,"         + c100(i1,i2,i3)*cr1 *d300i \\\n");
        fprintf(file,"         + c010(i1,i2,i3)*cs1 *d030i \n");
      else
        fprintf(file,"    %s = lap2h(i1,i2,i3,0) \\\n",lap4hName);
        fprintf(file,"         + c200(i1,i2,i3)*crr1*d400%s \\\n",suffix);
        fprintf(file,"         + c020(i1,i2,i3)*css1*d040%s \\\n",suffix);
        fprintf(file,"         + c002(i1,i2,i3)*ctt1*d004%s \\\n",suffix);
        fprintf(file,"         + c110(i1,i2,i3)*(cr1*d310i + cs1*d130i) \\\n");
        fprintf(file,"         + c101(i1,i2,i3)*(cr1*d301i + ct1*d103i) \\\n");
        fprintf(file,"         + c011(i1,i2,i3)*(cs1*d031i + ct1*d013i) \\\n");
        fprintf(file,"         + c100(i1,i2,i3)*cr1 *d300i \\\n");
        fprintf(file,"         + c010(i1,i2,i3)*cs1 *d030i \\\n");        
        fprintf(file,"         + c001(i1,i2,i3)*ct1 *d003i \n");        
      end if:

 
      fprintf(file,"\n    ! --- Laplacian squared to order 2:\n"); 
      lhsName:= saveDpDm( lap2h,lap2h,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap2h,lap2h,0,2,0,arrayNeededForNextLevel );  
      if nd=3 then 
      lhsName:= saveDpDm( lap2h,lap2h,0,0,2,arrayNeededForNextLevel ); 
      end if:

      lhsName:=saveD0( lap2h,lap2h,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap2h,lap2h,0,1,0,shortTargetName );    
      lhsName:=saveD0D0( lap2h,lap2h,1,1,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0( lap2h,lap2h,0,0,1,shortTargetName );    
      lhsName:=saveD0D0( lap2h,lap2h,1,0,1,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,0,1,1,shortTargetName );                     
      end if: 
      
      fprintf(file,"    %s =  \\\n",lap2hSqName);
      if nd=2 then
        fprintf(file,"             c200(i1,i2,i3)*lap2h200%s \\\n",suffix);
        fprintf(file,"           + c020(i1,i2,i3)*lap2h020%s \\\n",suffix);
        fprintf(file,"           + c110(i1,i2,i3)*lap2h110i  \\\n");
        fprintf(file,"           + c100(i1,i2,i3)*lap2h100i  \\\n");
        fprintf(file,"           + c010(i1,i2,i3)*lap2h010i    \n");
      else
        fprintf(file,"            c200(i1,i2,i3)*lap2h200%s \\\n",suffix);
        fprintf(file,"          + c020(i1,i2,i3)*lap2h020%s \\\n",suffix);
        fprintf(file,"          + c002(i1,i2,i3)*lap2h002%s \\\n",suffix);
        fprintf(file,"          + c110(i1,i2,i3)*lap2h110i  \\\n");
        fprintf(file,"          + c101(i1,i2,i3)*lap2h101i  \\\n");
        fprintf(file,"          + c011(i1,i2,i3)*lap2h011i  \\\n");
        fprintf(file,"          + c100(i1,i2,i3)*lap2h100i  \\\n");
        fprintf(file,"          + c010(i1,i2,i3)*lap2h010i  \\\n");
        fprintf(file,"          + c001(i1,i2,i3)*lap2h001i    \n");        
      end if:


    else

      # RECTANGULAR
      fprintf(file,"\n    ! --- Laplacian to order 4 = lap2h + corrections \n");
      fprintf(file,"    %s = lap2h(i1,i2,i3,0) \\\n",lap4hName);
      if nd=2 then
      fprintf(file,"         + cxx*cxx1*d400%s \\\n",suffix);
      fprintf(file,"         + cyy*cyy1*d040%s   \n",suffix);
      else
      fprintf(file,"         + cxx*cxx1*d400%s \\\n",suffix);
      fprintf(file,"         + cyy*cyy1*d040%s \\\n",suffix);
      fprintf(file,"         + czz*czz1*d004%s     ",suffix);
      end if:
      fprintf(file,"\n"):

      fprintf(file,"\n    ! --- Laplacian squared to order 2:\n"); 
      lhsName:= saveDpDm( lap2h,lap2h,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap2h,lap2h,0,2,0,arrayNeededForNextLevel ); 
      if nd=3 then 
      lhsName:= saveDpDm( lap2h,lap2h,0,0,2,arrayNeededForNextLevel ); 
      end if:
      
      fprintf(file,"    %s =                      \\\n",lap2hSqName);
      if nd=2 then
      fprintf(file,"             cxx*lap2h200%s   \\\n",suffix);
      fprintf(file,"           + cyy*lap2h020%s     \n",suffix);
      else
      fprintf(file,"             cxx*lap2h200%s   \\\n",suffix);
      fprintf(file,"           + cyy*lap2h020%s   \\\n",suffix);
      fprintf(file,"           + czz*lap2h002%s     \n",suffix);
      end if:

    end if:

    if ord=4 and orderOfAccuracy=4 then
      fprintf(file,"    #If #FORCING ne \"NOFORCING\" \n"):
      fprintf(file,"      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) \n");
      fprintf(file,"    #End \n\n"):       
      fprintf(file,"\n    ! --- Modified equation space-time update ----\n"):
      fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \\\n");
      fprintf(file,"             + cdtsq*( %s )                         \\\n",lap4hName);
      fprintf(file,"             + cdtPow4By12*( %s )                   \\\n",lap2hSqName);
      fprintf(file,"             FV(m)                                    \n"):
      # fprintf(file,"             + dtSq*fv(m)                             \n"):
    end if:    

  elif ord=6 then

    # ----------------- ORDER=6 ----------------
    # if ord<orderOfAccuracy then 
    #   shortTargetName:=0: 
    #   suffix := "(i1,i2,i3,0)":
    # else 
    #   shortTargetName:=1: 
    #   suffix := "i": 
    # end if: 

    if orderOfAccuracy=6 then
      # no need to store these:
      lap6hName := "lap6h": declareVar("lap6h"):
      lap4hSqName := "lap4hSq": declareVar("lap4hSq"):
      lap2hCubedName := "lap2hCubed": declareVar("lap2hCubed"):
    else
      # need to store these
      lap6hName := "lap6h(i1,i2,i3,0)": declareVar("lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):
      lap4hSqName := "lap4hSq(i1,i2,i3,0)": declareVar("lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):        
      lap2hCubedName := "lap2hCubed(i1,i2,i3,0)": declareVar("lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)"):        
    end if:      

    lhsName:= saveDpDm( d,d,6,0,0,arrayNeededForNextLevel ); 
    lhsName:= saveDpDm( d,d,0,6,0,arrayNeededForNextLevel ); 
    if nd=3 then
    lhsName:= saveDpDm( d,d,0,0,6,arrayNeededForNextLevel ); 
    end if:

    if gridType=curvilinear then

      if orderOfAccuracy>6 then
        lhsName:= saveDpDm( d,d,4,2,0,arrayNeededForNextLevel );  
        lhsName:= saveDpDm( d,d,2,4,0,arrayNeededForNextLevel );  
        if nd=3 then     
        lhsName:= saveDpDm( d,d,4,0,2,arrayNeededForNextLevel );  
        lhsName:= saveDpDm( d,d,2,0,4,arrayNeededForNextLevel );  
        lhsName:= saveDpDm( d,d,0,4,2,arrayNeededForNextLevel );  
        lhsName:= saveDpDm( d,d,0,2,4,arrayNeededForNextLevel );  
        end if:
      end if:      

      lhsName:=saveD0( d,d,5,0,0,shortTargetName );    
      lhsName:=saveD0( d,d,0,5,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( d,d,0,0,5,shortTargetName );
      end if:      

      lhsName:=saveD0D0( d,d,5,1,0,shortTargetName ); 
      lhsName:=saveD0D0( d,d,1,5,0,shortTargetName );  
      lhsName:=saveD0D0( d,d,3,3,0,shortTargetName );

      if nd=3 then
      lhsName:=saveD0D0( d,d,5,0,1,shortTargetName ); 
      lhsName:=saveD0D0( d,d,1,0,5,shortTargetName ); 
      lhsName:=saveD0D0( d,d,0,5,1,shortTargetName ); 
      lhsName:=saveD0D0( d,d,0,1,5,shortTargetName ); 

      lhsName:=saveD0D0( d,d,3,0,3,shortTargetName );   
      lhsName:=saveD0D0( d,d,0,3,3,shortTargetName );                  
      end if:      

      fprintf(file,"    ! --- Laplacian to order 6 = lap4h + corrections \n");
      fprintf(file,"    %s = lap4h(i1,i2,i3,0) \\\n",lap6hName);
      if nd=2 then
      fprintf(file,"       + c200(i1,i2,i3)*crr2*d600%s \\\n",suffix);
      fprintf(file,"       + c110(i1,i2,i3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) \\\n");
      fprintf(file,"       + c020(i1,i2,i3)*css2*d060%s \\\n",suffix);
      fprintf(file,"       + c100(i1,i2,i3)*cr2 *d500i \\\n");
      fprintf(file,"       + c010(i1,i2,i3)*cs2 *d050i \n");
      else
      fprintf(file,"       + c200(i1,i2,i3)*crr2*d600%s \\\n",suffix);
      fprintf(file,"       + c020(i1,i2,i3)*css2*d060%s \\\n",suffix);
      fprintf(file,"       + c002(i1,i2,i3)*ctt2*d006%s \\\n",suffix);
      fprintf(file,"       + c110(i1,i2,i3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) \\\n");
      fprintf(file,"       + c101(i1,i2,i3)*(cr2*d501i + ct2*d105i + cr1*ct1*d303i ) \\\n");
      fprintf(file,"       + c011(i1,i2,i3)*(cs2*d051i + ct2*d015i + cs1*ct1*d033i ) \\\n");
      fprintf(file,"       + c100(i1,i2,i3)*cr2 *d500i \\\n");
      fprintf(file,"       + c010(i1,i2,i3)*cs2 *d050i \\\n");
      fprintf(file,"       + c001(i1,i2,i3)*ct2 *d005i \n");        
      end if:

      # Slower: 
      # if splitLoops=1 then 
      #    endLoops( file ):
      #    fprintf(file," ! --- SPLIT LOOPS\n");
      #    startLoops( file ):
      #  end if:

      # ! --- Laplacian cubed order 2:
      # if ord<orderOfAccuracy then 
      #   shortTargetName:=0: 
      #   suffix := "(i1,i2,i3,0)":
      # else 
      #   shortTargetName:=1: 
      #   suffix := "i": 
      # end if:       
      lhsName:= saveDpDm( lap2hSq,lap2hSq,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,2,0,arrayNeededForNextLevel );  
      if nd=3 then
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,0,2,arrayNeededForNextLevel );
      end if:

      lhsName:=saveD0( lap2hSq,lap2hSq,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap2hSq,lap2hSq,0,1,0,shortTargetName ); 
      if nd=3 then   
      lhsName:=saveD0( lap2hSq,lap2hSq,0,0,1,shortTargetName ); 
      end if:
      lhsName:=saveD0D0( lap2hSq,lap2hSq,1,1,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0D0( lap2hSq,lap2hSq,1,0,1,shortTargetName );
      lhsName:=saveD0D0( lap2hSq,lap2hSq,0,1,1,shortTargetName );
      end if:

      fprintf(file,"    %s =  \\\n",lap2hCubedName);
      if nd=2 then
        fprintf(file,"       + c200(i1,i2,i3)*lap2hSq200%s \\\n",suffix);
        fprintf(file,"       + c110(i1,i2,i3)*lap2hSq110i  \\\n");
        fprintf(file,"       + c020(i1,i2,i3)*lap2hSq020%s \\\n",suffix);
        fprintf(file,"       + c100(i1,i2,i3)*lap2hSq100i  \\\n");
        fprintf(file,"       + c010(i1,i2,i3)*lap2hSq010i   \n"); 
      else
        fprintf(file,"       + c200(i1,i2,i3)*lap2hSq200%s \\\n",suffix);
        fprintf(file,"       + c020(i1,i2,i3)*lap2hSq020%s \\\n",suffix);
        fprintf(file,"       + c002(i1,i2,i3)*lap2hSq002%s \\\n",suffix);
        fprintf(file,"       + c110(i1,i2,i3)*lap2hSq110i  \\\n");
        fprintf(file,"       + c101(i1,i2,i3)*lap2hSq101i  \\\n");
        fprintf(file,"       + c011(i1,i2,i3)*lap2hSq011i  \\\n");
        fprintf(file,"       + c100(i1,i2,i3)*lap2hSq100i  \\\n");
        fprintf(file,"       + c010(i1,i2,i3)*lap2hSq010i  \\\n");
        fprintf(file,"       + c001(i1,i2,i3)*lap2hSq001i   \n");         
      end if:

      # Important split: 
      if splitLoops=1 then 
        endLoops( file ):
        fprintf(file," ! --- SPLIT LOOPS\n");
        startLoops( file ):
      end if:
      fprintf(file,"\n    ! --- Laplacian squared to order 4 = \n"):
      fprintf(file,"    !  lap2h*( lap4h ) + corrections*( Lap2h )\n"):

      # if ord<orderOfAccuracy then 
      #   shortTargetName:=0: 
      #   suffix := "(i1,i2,i3,0)":
      # else 
      #   shortTargetName:=1: 
      #   suffix := "i": 
      # end if:  
      lhsName:= saveDpDm( lap4h,lap4h,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap4h,lap4h,0,2,0,arrayNeededForNextLevel );  
      if nd=3 then
      lhsName:= saveDpDm( lap4h,lap4h,0,0,2,arrayNeededForNextLevel );  
      end if:
 
      lhsName:= saveDpDm( lap2h,lap2h,4,0,0,arrayNeededForNextLevel ); 
      lhsName:= saveDpDm( lap2h,lap2h,0,4,0,arrayNeededForNextLevel ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2h,lap2h,0,0,4,arrayNeededForNextLevel ); 
      end if:


      if orderOfAccuracy>6 then
        lhsName:= saveDpDm( lap2h,lap2h,2,2,0,arrayNeededForNextLevel ); 
        if nd=3 then
        lhsName:= saveDpDm( lap2h,lap2h,2,0,2,arrayNeededForNextLevel ); 
        lhsName:= saveDpDm( lap2h,lap2h,0,2,2,arrayNeededForNextLevel ); 
        end if:        
      end if:

      lhsName:=saveD0( lap4h,lap4h,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap4h,lap4h,0,1,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap4h,lap4h,0,0,1,shortTargetName );
      end if:

      lhsName:=saveD0D0( lap4h,lap4h,1,1,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0D0( lap4h,lap4h,1,0,1,shortTargetName ); 
      lhsName:=saveD0D0( lap4h,lap4h,0,1,1,shortTargetName ); 
      end if:

      lhsName:=saveD0( lap2h,lap2h,3,0,0,shortTargetName );    
      lhsName:=saveD0( lap2h,lap2h,0,3,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap2h,lap2h,0,0,3,shortTargetName );
      end if:
 
      lhsName:=saveD0D0( lap2h,lap2h,3,1,0,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,1,3,0,shortTargetName );  
      if nd=3 then
      lhsName:=saveD0D0( lap2h,lap2h,3,0,1,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,1,0,3,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,0,3,1,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,0,1,3,shortTargetName );  
      end if:

      fprintf(file,"    %s =     \\\n",lap4hSqName);
      if nd=2 then
        fprintf(file,"         c200(i1,i2,i3)*( lap4h200%s + crr1*lap2h400%s )    \\\n",suffix,suffix);
        fprintf(file,"       + c110(i1,i2,i3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) \\\n");
        fprintf(file,"       + c020(i1,i2,i3)*( lap4h020%s + css1*lap2h040%s )     \\\n",suffix,suffix);
        fprintf(file,"       + c100(i1,i2,i3)*( lap4h100i + cr1 *lap2h300i )    \\\n");
        fprintf(file,"       + c010(i1,i2,i3)*( lap4h010i + cs1 *lap2h030i )      \n");
      else
        fprintf(file,"         c200(i1,i2,i3)*( lap4h200%s + crr1*lap2h400%s )    \\\n",suffix,suffix);
        fprintf(file,"       + c020(i1,i2,i3)*( lap4h020%s + css1*lap2h040%s )     \\\n",suffix,suffix);
        fprintf(file,"       + c002(i1,i2,i3)*( lap4h002%s + ctt1*lap2h004%s )     \\\n",suffix,suffix);
        fprintf(file,"       + c110(i1,i2,i3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) \\\n");
        fprintf(file,"       + c101(i1,i2,i3)*( lap4h101i + cr1*lap2h301i + ct1*lap2h103i ) \\\n");
        fprintf(file,"       + c011(i1,i2,i3)*( lap4h011i + cs1*lap2h031i + ct1*lap2h013i ) \\\n");
        fprintf(file,"       + c100(i1,i2,i3)*( lap4h100i + cr1 *lap2h300i )    \\\n");
        fprintf(file,"       + c010(i1,i2,i3)*( lap4h010i + cs1 *lap2h030i )    \\\n");
        fprintf(file,"       + c001(i1,i2,i3)*( lap4h001i + ct1 *lap2h003i )      \n");        
      end if:



    else
      # RECTANGULAR
      fprintf(file,"    ! --- Laplacian to order 6 = lap4h + corrections \n");
      fprintf(file,"    %s = lap4h(i1,i2,i3,0) \\\n",lap6hName);
      if nd=2 then
        fprintf(file,"       + cxx*cxx2*d600%s \\\n",suffix);
        fprintf(file,"       + cyy*cyy2*d060%s   \n",suffix);
      else
        fprintf(file,"       + cxx*cxx2*d600%s \\\n",suffix);
        fprintf(file,"       + cyy*cyy2*d060%s \\\n",suffix);        
        fprintf(file,"       + czz*czz2*d006%s   \n",suffix);        
      end if:

      # ! --- Laplacian cubed order 2:
      lhsName:= saveDpDm( lap2hSq,lap2hSq,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,2,0,arrayNeededForNextLevel ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,0,2,arrayNeededForNextLevel ); 
      end if:

      fprintf(file,"    %s =                    \\\n",lap2hCubedName);
      if nd=2 then
        fprintf(file,"         cxx*lap2hSq200%s \\\n",suffix);
        fprintf(file,"       + cyy*lap2hSq020%s   \n",suffix);
      else
        fprintf(file,"         cxx*lap2hSq200%s  \\\n",suffix);
        fprintf(file,"       + cyy*lap2hSq020%s  \\\n",suffix);
        fprintf(file,"       + czz*lap2hSq020%s    \n",suffix);
      end if:

      fprintf(file,"\n    ! --- Laplacian squared to order 4 = \n"):
      fprintf(file,"    !  lap2h*( lap4h ) + corrections*( Lap2h )\n"):

      lhsName:= saveDpDm( lap4h,lap4h,2,0,0,arrayNeededForNextLevel );  
      lhsName:= saveDpDm( lap4h,lap4h,0,2,0,arrayNeededForNextLevel ); 
      if nd=3 then
      lhsName:= saveDpDm( lap4h,lap4h,0,0,2,arrayNeededForNextLevel ); 
      end if:       

      lhsName:= saveDpDm( lap2h,lap2h,4,0,0,arrayNeededForNextLevel ); 
      lhsName:= saveDpDm( lap2h,lap2h,0,4,0,arrayNeededForNextLevel );
      if nd=3 then
      lhsName:= saveDpDm( lap2h,lap2h,0,0,4,arrayNeededForNextLevel );
      end if:       
      # if orderOfAccuracy>6 then
      #   lhsName:= saveDpDm( lap2h,lap2h,2,2,0,shortTargetName ); 
      # end if: 

      fprintf(file,"    %s =                                       \\\n",lap4hSqName);
      if nd=2 then
        fprintf(file,"         cxx*( lap4h200%s + cxx1*lap2h400%s )  \\\n",suffix,suffix);
        fprintf(file,"       + cyy*( lap4h020%s + cyy1*lap2h040%s )    \n",suffix,suffix);
      else
        fprintf(file,"         cxx*( lap4h200%s + cxx1*lap2h400%s )  \\\n",suffix,suffix);
        fprintf(file,"       + cyy*( lap4h020%s + cyy1*lap2h040%s )  \\\n",suffix,suffix);
        fprintf(file,"       + czz*( lap4h002%s + czz1*lap2h004%s )    \n",suffix,suffix);
      end if:


    end if:

    if ord=6 and orderOfAccuracy=6 then
      fprintf(file,"    #If #FORCING ne \"NOFORCING\" \n"):
      fprintf(file,"      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) \n");
      fprintf(file,"    #End \n\n"): 
      fprintf(file,"\n    ! --- Modified equation space-time update ----\n"):

      fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  \\\n");
      fprintf(file,"             + cdtsq*( %s )               \\\n",lap6hName);
      fprintf(file,"             + cdtPow4By12*( %s )       \\\n",lap4hSqName);
      fprintf(file,"             + cdtPow6By360*( %s )   \\\n",lap2hCubedName):
      fprintf(file,"             FV(m)                    \n"):
      # fprintf(file,"             + dtSq*fv(m)                      \n"):        
    end if:


  elif ord=8 then
    # ----------------- ORDER=8 ---------------- 

    
    lhsName:= saveDpDm( d,d,8,0,0,shortTargetName ); 
    lhsName:= saveDpDm( d,d,0,8,0,shortTargetName ); 
    if nd=3 then
    lhsName:= saveDpDm( d,d,0,0,8,shortTargetName ); 
    end if:   

    if gridType=curvilinear then

      lhsName:=saveD0( d,d,7,0,0,shortTargetName );    
      lhsName:=saveD0( d,d,0,7,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( d,d,0,0,7,shortTargetName );
      end if:       

      lhsName:=saveD0D0( d,d,7,1,0,shortTargetName );  
      lhsName:=saveD0D0( d,d,1,7,0,shortTargetName );  
      if nd=3 then
      lhsName:=saveD0D0( d,d,7,0,1,shortTargetName );  
      lhsName:=saveD0D0( d,d,1,0,7,shortTargetName );  
      lhsName:=saveD0D0( d,d,0,7,1,shortTargetName );  
      lhsName:=saveD0D0( d,d,0,1,7,shortTargetName );  
      end if:       
      lhsName:=saveD0D0( d,d,5,3,0,shortTargetName ); 
      lhsName:=saveD0D0( d,d,3,5,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0D0( d,d,5,0,3,shortTargetName ); 
      lhsName:=saveD0D0( d,d,3,0,5,shortTargetName );
      lhsName:=saveD0D0( d,d,0,5,3,shortTargetName ); 
      lhsName:=saveD0D0( d,d,0,3,5,shortTargetName );
      end if:        

      declareVar("lap8h"):
      fprintf(file,"    ! --- Laplacian to order 8 = lap6h + corrections \n");
      fprintf(file,"    lap8h = lap6h(i1,i2,i3,0)                                                         \\\n");
      if nd=2 then
        fprintf(file,"            + c200(i1,i2,i3)*crr3*d800i                                               \\\n");
        fprintf(file,"            + c110(i1,i2,i3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) \\\n");
        fprintf(file,"            + c020(i1,i2,i3)*css3*d080i                                               \\\n");
        fprintf(file,"            + c100(i1,i2,i3)* cr3*d700i                                               \\\n");
        fprintf(file,"            + c010(i1,i2,i3)* cs3*d070i \n");
      else
        fprintf(file,"            + c200(i1,i2,i3)*crr3*d800i                                               \\\n");
        fprintf(file,"            + c020(i1,i2,i3)*css3*d080i                                               \\\n");
        fprintf(file,"            + c002(i1,i2,i3)*ctt3*d008i                                               \\\n");
        fprintf(file,"            + c110(i1,i2,i3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) \\\n");
        fprintf(file,"            + c101(i1,i2,i3)*(cr3*d701i + ct3*d107i + cr2*ct1*d503i + cr1*ct2*d305i ) \\\n");
        fprintf(file,"            + c011(i1,i2,i3)*(cs3*d071i + ct3*d017i + cs2*ct1*d053i + cs1*ct2*d035i ) \\\n");
        fprintf(file,"            + c100(i1,i2,i3)* cr3*d700i                                               \\\n");
        fprintf(file,"            + c010(i1,i2,i3)* cs3*d070i                                               \\\n");
        fprintf(file,"            + c001(i1,i2,i3)* ct3*d007i \n");        
      end if:


      fprintf(file,"\n    ! --- Laplacian^4 4p (4th power) order 2: \n"):

      if ord<orderOfAccuracy then shortTargetName:=0: else shortTargetName:=1: end if:

      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,0,2,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,0,0,2,shortTargetName ); 
      end if:        
      shortTargetName:=1: 
      lhsName:=saveD0( lap2hCubed,lap2hCubed,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap2hCubed,lap2hCubed,0,1,0,shortTargetName );   
      if nd=3 then
      lhsName:=saveD0( lap2hCubed,lap2hCubed,0,0,1,shortTargetName );   
      end if:        
      lhsName:=saveD0D0( lap2hCubed,lap2hCubed,1,1,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0D0( lap2hCubed,lap2hCubed,1,0,1,shortTargetName );
      lhsName:=saveD0D0( lap2hCubed,lap2hCubed,0,1,1,shortTargetName );
      end if:       

      declareVar("lap2h4p"):
      fprintf(file,"    lap2h4p  =                             \\\n"):
      if nd=2 then
        fprintf(file,"             + c200(i1,i2,i3)*lap2hCubed200i  \\\n"):
        fprintf(file,"             + c110(i1,i2,i3)*lap2hCubed110i  \\\n"):
        fprintf(file,"             + c020(i1,i2,i3)*lap2hCubed020i  \\\n"):
        fprintf(file,"             + c100(i1,i2,i3)*lap2hCubed100i  \\\n"):
        fprintf(file,"             + c010(i1,i2,i3)*lap2hCubed010i    \n"): 
      else
        fprintf(file,"             + c200(i1,i2,i3)*lap2hCubed200i  \\\n"):
        fprintf(file,"             + c020(i1,i2,i3)*lap2hCubed020i  \\\n"):
        fprintf(file,"             + c002(i1,i2,i3)*lap2hCubed002i  \\\n"):
        fprintf(file,"             + c110(i1,i2,i3)*lap2hCubed110i  \\\n"):
        fprintf(file,"             + c101(i1,i2,i3)*lap2hCubed101i  \\\n"):
        fprintf(file,"             + c011(i1,i2,i3)*lap2hCubed011i  \\\n"):
        fprintf(file,"             + c100(i1,i2,i3)*lap2hCubed100i  \\\n"):
        fprintf(file,"             + c010(i1,i2,i3)*lap2hCubed010i  \\\n"):
        fprintf(file,"             + c001(i1,i2,i3)*lap2hCubed001i    \n"):         
      end if: 


      fprintf(file,"    ! --- Laplacian squared to order 6 :\n"):
      fprintf(file,"    !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 \n"):
      fprintf(file,"    !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]\n"):
      fprintf(file,"    !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)\n"):


      lhsName:= saveDpDm( lap6h,lap6h,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap6h,lap6h,0,2,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap6h,lap6h,0,0,2,shortTargetName ); 
      end if: 

      lhsName:=saveD0( lap6h,lap6h,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap6h,lap6h,0,1,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap6h,lap6h,0,0,1,shortTargetName );
      end if:           
      lhsName:=saveD0D0( lap6h,lap6h,1,1,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0D0( lap6h,lap6h,1,0,1,shortTargetName ); 
      lhsName:=saveD0D0( lap6h,lap6h,0,1,1,shortTargetName ); 
      end if:        

      lhsName:= saveDpDm( lap4h,lap4h,4,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap4h,lap4h,0,4,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap4h,lap4h,0,0,4,shortTargetName ); 
      end if:       
      # lhsName:= saveDpDm( lap4h,lap4h,2,2,0,shortTargetName ); 

      lhsName:=saveD0( lap4h,lap4h,3,0,0,shortTargetName );    
      lhsName:=saveD0( lap4h,lap4h,0,3,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap4h,lap4h,0,0,3,shortTargetName );
      end if:       

      lhsName:=saveD0D0( lap4h,lap4h,3,1,0,shortTargetName );  
      lhsName:=saveD0D0( lap4h,lap4h,1,3,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0D0( lap4h,lap4h,3,0,1,shortTargetName ); 
      lhsName:=saveD0D0( lap4h,lap4h,1,0,3,shortTargetName ); 
      lhsName:=saveD0D0( lap4h,lap4h,0,3,1,shortTargetName ); 
      lhsName:=saveD0D0( lap4h,lap4h,0,1,3,shortTargetName ); 
      end if: 

      lhsName:= saveDpDm( lap2h,lap2h,6,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap2h,lap2h,0,6,0,shortTargetName );       
      if nd=3 then
      lhsName:= saveDpDm( lap2h,lap2h,0,0,6,shortTargetName ); 
      end if: 

      lhsName:=saveD0( lap2h,lap2h,5,0,0,shortTargetName );    
      lhsName:=saveD0( lap2h,lap2h,0,5,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap2h,lap2h,0,0,5,shortTargetName );
      end if: 

      lhsName:=saveD0D0( lap2h,lap2h,5,1,0,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,1,5,0,shortTargetName );       
      lhsName:=saveD0D0( lap2h,lap2h,3,3,0,shortTargetName );       
      if nd=3 then
      lhsName:=saveD0D0( lap2h,lap2h,5,0,1,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,1,0,5,shortTargetName );        
      lhsName:=saveD0D0( lap2h,lap2h,0,5,1,shortTargetName );  
      lhsName:=saveD0D0( lap2h,lap2h,0,1,5,shortTargetName ); 
      lhsName:=saveD0D0( lap2h,lap2h,3,0,3,shortTargetName );       
      lhsName:=saveD0D0( lap2h,lap2h,0,3,3,shortTargetName );       
      end if: 

      declareVar("lap6hSq"):
      fprintf(file,"    lap6hSq =                                                                                     \\\n"):
      if nd=2 then
        fprintf(file,"              c200(i1,i2,i3)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       \\\n"):
        fprintf(file,"            + c110(i1,i2,i3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         \\\n"):
        fprintf(file,"                                        +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   \\\n"):
        fprintf(file,"            + c020(i1,i2,i3)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       \\\n"):
        fprintf(file,"            + c100(i1,i2,i3)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       \\\n"):
        fprintf(file,"            + c010(i1,i2,i3)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                         \n"):
      else
        fprintf(file,"              c200(i1,i2,i3)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       \\\n"):
        fprintf(file,"            + c020(i1,i2,i3)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       \\\n"):
        fprintf(file,"            + c002(i1,i2,i3)*(lap6h002i + ctt1*lap4h004i + ctt2*lap2h006i )                       \\\n"):
        fprintf(file,"            + c110(i1,i2,i3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         \\\n"):
        fprintf(file,"                                        +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   \\\n"):
        fprintf(file,"            + c101(i1,i2,i3)*(lap6h101i +  cr1*lap4h301i +  cr2*lap2h501i                         \\\n"):
        fprintf(file,"                                        +  ct1*lap4h103i +  ct2*lap2h105i + cr1*ct1*lap2h303i )   \\\n"):
        fprintf(file,"            + c011(i1,i2,i3)*(lap6h011i +  cs1*lap4h031i +  cs2*lap2h051i                         \\\n"):
        fprintf(file,"                                        +  ct1*lap4h013i +  ct2*lap2h015i + cs1*ct1*lap2h033i )   \\\n"):
        fprintf(file,"            + c100(i1,i2,i3)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       \\\n"):
        fprintf(file,"            + c010(i1,i2,i3)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                       \\\n"):
        fprintf(file,"            + c001(i1,i2,i3)*(lap6h010i +  ct1*lap4h003i +  ct2*lap2h005i )                         \n"):        
      end if:


      fprintf(file,"\n    ! --- Laplacian CUBED to order 4 \n"):

      lhsName:= saveDpDm( lap4hSq,lap4hSq,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap4hSq,lap4hSq,0,2,0,shortTargetName );  
      if nd=3 then
      lhsName:= saveDpDm( lap4hSq,lap4hSq,0,0,2,shortTargetName );  
      end if:

      lhsName:=saveD0( lap4hSq,lap4hSq,1,0,0,shortTargetName );    
      lhsName:=saveD0( lap4hSq,lap4hSq,0,1,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap4hSq,lap4hSq,0,0,1,shortTargetName );
      end if:           
      lhsName:=saveD0D0( lap4hSq,lap4hSq,1,1,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0D0( lap4hSq,lap4hSq,1,0,1,shortTargetName ); 
      lhsName:=saveD0D0( lap4hSq,lap4hSq,0,1,1,shortTargetName ); 
      end if:        

      lhsName:= saveDpDm( lap2hSq,lap2hSq,4,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,4,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,0,4,shortTargetName ); 
      end if:       
      # lhsName:= saveDpDm( lap2hSq,lap2hSq,2,2,0,shortTargetName ); 

      lhsName:=saveD0( lap2hSq,lap2hSq,3,0,0,shortTargetName );    
      lhsName:=saveD0( lap2hSq,lap2hSq,0,3,0,shortTargetName );
      if nd=3 then
      lhsName:=saveD0( lap2hSq,lap2hSq,0,0,3,shortTargetName );
      end if:       

      lhsName:=saveD0D0( lap2hSq,lap2hSq,3,1,0,shortTargetName );  
      lhsName:=saveD0D0( lap2hSq,lap2hSq,1,3,0,shortTargetName ); 
      if nd=3 then
      lhsName:=saveD0D0( lap2hSq,lap2hSq,3,0,1,shortTargetName );  
      lhsName:=saveD0D0( lap2hSq,lap2hSq,1,0,3,shortTargetName );  
      lhsName:=saveD0D0( lap2hSq,lap2hSq,0,3,1,shortTargetName );  
      lhsName:=saveD0D0( lap2hSq,lap2hSq,0,1,3,shortTargetName );  
      end if:       

      declareVar("lap4hCubed"):
      fprintf(file,"    lap4hCubed =                                                    \\\n"):
      if nd=2 then
        fprintf(file,"                 c200(i1,i2,i3)*(lap4hSq200i + crr1*lap2hSq400i )   \\\n"):
        fprintf(file,"               + c110(i1,i2,i3)*(lap4hSq110i +  cr1*lap2hSq310i     \\\n"):
        fprintf(file,"                                             +  cs1*lap2hSq130i )   \\\n"):
        fprintf(file,"               + c020(i1,i2,i3)*(lap4hSq020i + css1*lap2hSq040i )   \\\n"):
        fprintf(file,"               + c100(i1,i2,i3)*(lap4hSq100i + cr1 *lap2hSq300i )   \\\n"):
        fprintf(file,"               + c010(i1,i2,i3)*(lap4hSq010i + cs1 *lap2hSq030i )     \n"):
      else
        fprintf(file,"                 c200(i1,i2,i3)*(lap4hSq200i + crr1*lap2hSq400i )   \\\n"):
        fprintf(file,"               + c020(i1,i2,i3)*(lap4hSq020i + css1*lap2hSq040i )   \\\n"):
        fprintf(file,"               + c002(i1,i2,i3)*(lap4hSq002i + ctt1*lap2hSq004i )   \\\n"):
        fprintf(file,"               + c110(i1,i2,i3)*(lap4hSq110i +  cr1*lap2hSq310i     \\\n"):
        fprintf(file,"                                             +  cs1*lap2hSq130i )   \\\n"):
        fprintf(file,"               + c101(i1,i2,i3)*(lap4hSq101i +  cr1*lap2hSq301i     \\\n"):
        fprintf(file,"                                             +  ct1*lap2hSq103i )   \\\n"):
        fprintf(file,"               + c011(i1,i2,i3)*(lap4hSq011i +  cs1*lap2hSq031i     \\\n"):
        fprintf(file,"                                             +  ct1*lap2hSq013i )   \\\n"):
        fprintf(file,"               + c100(i1,i2,i3)*(lap4hSq100i + cr1 *lap2hSq300i )   \\\n"):
        fprintf(file,"               + c010(i1,i2,i3)*(lap4hSq010i + cs1 *lap2hSq030i )   \\\n"):
        fprintf(file,"               + c001(i1,i2,i3)*(lap4hSq001i + ct1 *lap2hSq003i )     \n"):        
      end if:



    else

       # RECTANGULAR
      declareVar("lap8h"):
      fprintf(file,"    ! --- Laplacian to order 8 = lap6h + corrections \n");
      fprintf(file,"    lap8h = lap6h(i1,i2,i3,0)  \\\n");
      if nd=2 then
        fprintf(file,"            + cxx*cxx3*d800i  \\\n");
        fprintf(file,"            + cyy*cyy3*d080i    \n");
      else
        fprintf(file,"            + cxx*cxx3*d800i  \\\n");
        fprintf(file,"            + cyy*cyy3*d080i  \\\n");
        fprintf(file,"            + czz*czz3*d008i    \n");
      end if:


      fprintf(file,"\n    ! --- Laplacian^4 4p (4th power) order 2: \n"):

      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,0,2,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2hCubed,lap2hCubed,0,0,2,shortTargetName ); 
      end if:        

      declareVar("lap2h4p"):
      fprintf(file,"    lap2h4p =                   \\\n"):
      if nd=2 then
        fprintf(file,"            + cxx*lap2hCubed200i   \\\n"):
        fprintf(file,"            + cyy*lap2hCubed020i     \n"):
      else
        fprintf(file,"            + cxx*lap2hCubed200i   \\\n"):
        fprintf(file,"            + cyy*lap2hCubed020i   \\\n"):
        fprintf(file,"            + czz*lap2hCubed002i     \n"):
      end if:

      fprintf(file,"    ! --- Laplacian squared to order 6 :\n"):
      fprintf(file,"    !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 \n"):
      fprintf(file,"    !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]\n"):
      fprintf(file,"    !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)\n"):


      lhsName:= saveDpDm( lap6h,lap6h,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap6h,lap6h,0,2,0,shortTargetName );  
      if nd=3 then
      lhsName:= saveDpDm( lap6h,lap6h,0,0,2,shortTargetName );  
      end if:       

      lhsName:= saveDpDm( lap4h,lap4h,4,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap4h,lap4h,0,4,0,shortTargetName );
      if nd=3 then
      lhsName:= saveDpDm( lap4h,lap4h,0,0,4,shortTargetName );
      end if:        
      # lhsName:= saveDpDm( lap4h,lap4h,2,2,0,shortTargetName ); 

      lhsName:= saveDpDm( lap2h,lap2h,6,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap2h,lap2h,0,6,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2h,lap2h,0,0,6,shortTargetName ); 
      end if:       

      declareVar("lap6hSq"):
      fprintf(file,"    lap6hSq =                                                    \\\n"):
      if nd=2 then
        fprintf(file,"              cxx*(lap6h200i + cxx1*lap4h400i + cxx2*lap2h600i ) \\\n"):
        fprintf(file,"            + cyy*(lap6h020i + cyy1*lap4h040i + cyy2*lap2h060i )   \n"):
      else
        fprintf(file,"              cxx*(lap6h200i + cxx1*lap4h400i + cxx2*lap2h600i ) \\\n"):
        fprintf(file,"            + cyy*(lap6h020i + cyy1*lap4h040i + cyy2*lap2h060i ) \\\n"):
        fprintf(file,"            + czz*(lap6h002i + czz1*lap4h004i + czz2*lap2h006i )   \n"):
      end if:


      fprintf(file,"\n    ! --- Laplacian CUBED to order 4 \n"):

      lhsName:= saveDpDm( lap4hSq,lap4hSq,2,0,0,shortTargetName );  
      lhsName:= saveDpDm( lap4hSq,lap4hSq,0,2,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap4hSq,lap4hSq,0,0,2,shortTargetName ); 
      end if:        

      lhsName:= saveDpDm( lap2hSq,lap2hSq,4,0,0,shortTargetName ); 
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,4,0,shortTargetName ); 
      if nd=3 then
      lhsName:= saveDpDm( lap2hSq,lap2hSq,0,0,4,shortTargetName ); 
      end if:       
      # lhsName:= saveDpDm( lap2hSq,lap2hSq,2,2,0,shortTargetName ); 

      declareVar("lap4hCubed"):
      fprintf(file,"    lap4hCubed =                                         \\\n"):
      if nd=2 then
        fprintf(file,"                 cxx*(lap4hSq200i + cxx1*lap2hSq400i )   \\\n"):
        fprintf(file,"               + cyy*(lap4hSq020i + cyy1*lap2hSq040i )     \n"):
      else
        fprintf(file,"                 cxx*(lap4hSq200i + cxx1*lap2hSq400i )   \\\n"):
        fprintf(file,"               + cyy*(lap4hSq020i + cyy1*lap2hSq040i )   \\\n"):
        fprintf(file,"               + czz*(lap4hSq002i + czz1*lap2hSq004i )     \n"):
      end if:

    end if:

    if ord=8 and orderOfAccuracy=8 then

      fprintf(file,"    #If #FORCING ne \"NOFORCING\" \n"):
      fprintf(file,"      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) \n");
      fprintf(file,"    #End \n\n"): 

      fprintf(file,"\n    ! --- Modified equation space-time update ----\n"):
      fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  \\\n");
      fprintf(file,"             + cdtsq*( lap8h )               \\\n");
      fprintf(file,"             + cdtPow4By12*( lap6hSq )       \\\n");
      fprintf(file,"             + cdtPow6By360*( lap4hCubed )   \\\n"):
      fprintf(file,"             + cdtPow8By20160*( lap2h4p )    \\\n"):
      fprintf(file,"             FV(m)                    \n"):
      # fprintf(file,"             + dtSq*fv(m)                      \n"):

    end if:    

  end if:

  endLoops( file ):
  # fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
  # fprintf(file,"  end if ! mask .ne. 0\n");
  # fprintf(file,"  #End \n"):
  # fprintf(file,"endLoops3d() \n");
  
end do:








  fprintf(file,"#endMacro"); 
  declareVar(finishDeclarations); # finish declarations
  fprintf(dfile,"#endMacro"); 
  fclose(file);
  fclose(dfile);
  printf("Wrote file=%s\n", myFileName);
  printf("Wrote file=%s\n", declareFileName);
  

end do:
end do: # end orderOfAccuracy
end do: # end gridType
# -------------------------------------------------------------------------------
# -------------- END LOOP OVER ORDER AND DIMENSIONS AND GRIDTYPE ---------------
# -------------------------------------------------------------------------------
