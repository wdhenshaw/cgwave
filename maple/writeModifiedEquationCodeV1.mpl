#  
# Generate code for the modified equation time-stepping
#    VERSION 1 USING HIERACHICAL 
#      -- CURVLINEAR VERSION IS VERY ***SLOW**
#
# restart: currentdir("/Users/henshaw/Dropbox/research/cgwave"): currentdir("/Users/henshaw/Dropbox/research/cgwave/maple"): read "writeModifiedEquationCodeV1.mpl";
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
# orderOfAccuracy  := 6:      # order of accuracy
orderStart:=2: orderEnd:=8: 
gridType         := rectangular: 
# gridType         := curvilinear;

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
  # Weights in deriv 7 and 8
  dw[0,7]:=1;   dw[1,7]:=0; dw[2,7]:=0; dw[3,7]:=0; dw[4,7]:=0; dw[5,7]:=0;   # *** FIX ME IF WE GO TO HIGHER ODRER THAN 8
  dw[0,8]:=1;   dw[1,8]:=0; dw[2,8]:=0; dw[3,8]:=0; dw[4,8]:=0; dw[5,8]:=0;   # *** FIX ME IF WE GO TO HIGHER ODRER THAN 8
end if:

printf("Read definitions of chain rule derivatives...\n");
read "chainRuleCoefficientsDefinitions2d.maple":
read "chainRuleCoefficientsDefinitions3d.maple":
printf("...done\n");

  # ---------- START PROCEDURE saveDpDm --------
  #
  # Auxillary function to print a second difference formulae
  # Example:
  #    saveDpDm( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveDpDm :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c;

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
    if target=source then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    # fprintf(file,"    d400(i1,i2,i3,0) = (d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0))/(dx(0)**2)\n",
    
    if d1>0  then
      # give preferences for differences in the first direction (better memory accesses)
      # d1m := d1-2; 
      fprintf(file,"    %s = %s(i1+1,i2,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1-1,i2,i3,%s)\n",lhs,rhs,c,rhs,c,rhs,c );
      # fprintf(file,"    %s = (%s(i1+1,i2,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1-1,i2,i3,%s))/(dx(0)**2)\n",lhs,rhs,c,rhs,c,rhs,c );
    elif d2>0 then
      # d2m:=d2-2;
      fprintf(file,"    %s = %s(i1,i2+1,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1,i2-1,i3,%s)\n",lhs,rhs,c,rhs,c,rhs,c ); 
      # fprintf(file,"    %s = (%s(i1,i2+1,i3,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1,i2-1,i3,%s))/(dx(1)**2)\n",lhs,rhs,c,rhs,c,rhs,c ); 
    elif d3>0 then 
      # d3m:=d3-2;   
      fprintf(file,"    %s = %s(i1,i2,i3+1,%s) - 2.*%s(i1,i2,i3,%s) + %s(i1,i2,i3-1,%s)\n",lhs,rhs,c,rhs,c,rhs,c );
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
  # Auxillary function to print a FIRST difference formulae
  # Example:
  #    saveD0( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveD0 :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c;

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
    if target=source then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    if d1>0 and (d1 mod 2 =1) then
      # give preferences for differences in the first direction (better memory accesses)
      fprintf(file,"    %s = (%s(i1+1,i2,i3,%s) - %s(i1-1,i2,i3,%s))/(2.*dx(0))\n",lhs,rhs,c,rhs,c );
    elif d2>0 and (d2 mod 2 =1) then
      fprintf(file,"    %s = (%s(i1,i2+1,i3,%s) - %s(i1,i2-1,i3,%s))/(2.*dx(1))\n",lhs,rhs,c,rhs,c ); 
    elif d3>0 and (d3 mod 2 =1) then 
      fprintf(file,"    %s = (%s(i1,i2,i3+1,%s) - %s(i1,i2,i3-1,%s))/(2.*dx(2))\n",lhs,rhs,c,rhs,c );
    else
      printf("saveDpDm: ERROR: invalid d1,d2,d3\n");
      error("stop here");    
    end if;
    
    RETURN(lhs):
  end:
  # ---------- END PROCEDURE save D0 --------

  # ---------- START PROCEDURE saveD0D0 --------
  #
  # Auxillary function to print a MIXED DERIVATIVE difference formulae
  # Example:
  #    saveD0D0( d,d,4,2,0 )
  #
  # target,source (input) : names of lhs and rhs variables
  # d1,d2,d3  (Input)     : number of derivatives in x,y,z or r,s,t
  #
  saveD0D0 :=proc( target,source,d1,d2,d3, shortTargetName:=0 )
    local lhs,rhs,d1m,d2m,d3m,c;

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
    if target=source then
      rhs := sprintf("%s%d%d%d",convert(source,string),d1m,d2m,d3m); # RHS name
    else
      rhs := convert(source,string);
    end if;

    if d1>0 and (d1 mod 2 =1) and d2>0 and (d2 mod 2 = 1) then
      # give preferences for differences in the first direction (better memory accesses)
      fprintf(file,"    %s = (%s(i1+1,i2+1,i3,%s) - %s(i1-1,i2+1,i3,%s) - %s(i1+1,i2-1,i3,%s) + %s(i1-1,i2-1,i3,%s))/(4.*dx(0)*dx(1))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
    elif d1>0 and (d1 mod 2 =1) and d3>0 and (d3 mod 2 = 1) then
      fprintf(file,"    %s = (%s(i1+1,i2,i3+1,%s) - %s(i1-1,i2,i3+1,%s) - %s(i1+1,i2,i3-1,%s) + %s(i1-1,i2,i3-1,%s))/(4.*dx(0)*dx(2))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
    elif d2>0 and (d2 mod 2 =1) and d3>0 and (d3 mod 2 = 1) then
      fprintf(file,"    %s = (%s(i1,i2+1,i3+1,%s) - %s(i1,i2-1,i3+1,%s) - %s(i1,i2+1,i3-1,%s) + %s(i1,i2-1,i3-1,%s))/(4.*dx(1)*dx(2))\n",lhs,rhs,c,rhs,c,rhs,c,rhs,c );
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
  # ---------- END PROCEDURE save D0D0D0 --------



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

  # ---------- START PROCEDURE defineDifferences --------
  #
  # Auxillary function to define even and odd differences of u or rsxy in the main loop
  #
  defineDifferences :=proc( d,u );
    global file,nd,gridType,orderOfAccuracy,rectangular,curvilinear:
    local d1,d2,d3,shortTargetName,numOdd,dSum,lhsName,ord,ord3,isOdd1,isOdd2,isOdd3:

    shortTargetName:=1; # use a short name for LHS
    ord := orderOfAccuracy:
    if nd=2 then ord3:=0; else ord3:=ord; end if;

    for d3 from 0 by 1 to ord3 do
    for d2 from 0 by 1 to ord do
    for d1 from 0 by 1 to ord do
      if d1 mod 2 =1 then isOdd1:=1; else isOdd1:=0; end if:
      if d2 mod 2 =1 then isOdd2:=1; else isOdd2:=0; end if:
      if d3 mod 2 =1 then isOdd3:=1; else isOdd3:=0; end if:  
      numOdd:= isOdd1 + isOdd2 + isOdd3:
      dSum := d1 + d2 + d3:
      if numOdd=0 then
        # all even derivatives
        if dSum=ord then
          lhsName := saveDpDm( d,d,d1,d2,d3,shortTargetName );                     # higest even derivative computed with D+D-
          lhsName := SubstituteAll(lhsName,"(m1,m2)",sprintf("(0:%d,0:%d)",nd-1,nd-1)):
          declareVar(lhsName); 
        elif dSum>0 and dSum<ord then
          # lower even derivatives have already been computed
          if convert(d,string)="rx" then
            fprintf(file,"    %s%d%d%di(m1,m2) = %s%d%d%d(i1,i2,i3,m1,m2)\n",d,d1,d2,d3,d,d1,d2,d3);  # set d020i = d020(i1,i2,i3,0)
            declareVar(sprintf("%s%d%d%di(0:%d,0:%d)",d,d1,d2,d3,nd-1,nd-1));  
           else 
            fprintf(file,"    %s%d%d%di = %s%d%d%d(i1,i2,i3,0)\n",d,d1,d2,d3,d,d1,d2,d3);  # set d020i = d020(i1,i2,i3,0)
            declareVar(sprintf("%s%d%d%di",d,d1,d2,d3));             
          end if   
        end if:

      elif numOdd=1 and dSum<=ord and gridType=curvilinear  then
        # one odd derivative
        if dSum=1 then
          lhsName:=saveD0( d,u,d1,d2,d3,shortTargetName );  # rhs uses "u"
        else
          lhsName:=saveD0( d,d,d1,d2,d3,shortTargetName );
        end if: 
        lhsName := SubstituteAll(lhsName,"(m1,m2)",sprintf("(0:%d,0:%d)",nd-1,nd-1)):
        declareVar(lhsName);  

      elif numOdd=2 and dSum<=ord and gridType=curvilinear then
        # two odd derivatives
        if dSum=2 then
          lhsName:=saveD0D0( d,u,d1,d2,d3,shortTargetName );  # rhs uses "u"
        else
          lhsName:=saveD0D0( d,d,d1,d2,d3,shortTargetName );
        end if: 
        lhsName := SubstituteAll(lhsName,"(m1,m2)",sprintf("(0:%d,0:%d)",nd-1,nd-1)):
        declareVar(lhsName);        

      elif numOdd=3 and dSum<=ord and gridType=curvilinear then
        # three odd derivatives
        if dSum=3 then
          lhsName:=saveD0D0D0( d,u,d1,d2,d3,shortTargetName );  # rhs uses "u"
        else
          lhsName:=saveD0D0D0( d,d,d1,d2,d3,shortTargetName );
        end if: 
        lhsName := SubstituteAll(lhsName,"(m1,m2)",sprintf("(0:%d,0:%d)",nd-1,nd-1)):
        declareVar(lhsName);  
      elif  dSum<ord and gridType=curvilinear then
        printf("ERROR: this case should not happen");
        error();     
      end if:
    
    end do:
    end do:
    end do:

  end:
  # ---------- END PROCEDURE defineDifferences --------

# -----------------------------------------------------------------
# ------------ START LOOP OVER ORDER AND DIMENSIONS ---------------
# -----------------------------------------------------------------
for orderOfAccuracy from orderStart by 2 to orderEnd do 
for nd from ndStart to ndEnd do 

 declareCount:=0;         # global variable for counting declared variables
 numDeclaredPerLine:=10;  # write this many declared variables per line  


# -- define file names: 
myFileName := sprintf("update%ddOrder%d%sV1.h",nd,orderOfAccuracy,gridTypeName[gridType]); 
declareFileName := sprintf("declare%ddOrder%d%sV1.h",nd,orderOfAccuracy,gridTypeName[gridType]);  # for declarations of variables 

#else
#  myFileName := sprintf("update%ddOrder%dCurvilinear.h",nd,orderOfAccuracy,gridType); 
#  declareFileName := sprintf("declare%ddOrder%dCurvilinear.h",nd,orderOfAccuracy,gridType); 
#end if:

file := open(myFileName, WRITE);
dfile := open(declareFileName, WRITE);  # for declarations of variables 

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

if gridType = rectangular then
  fprintf(file,"cxx=1./dx(0)**2;\n");
  fprintf(file,"cyy=1./dx(1)**2;\n");
  fprintf(file,"czz=1./dx(2)**2;\n");
  declareVar("cxx");  declareVar("cyy"):  declareVar("czz"):
end if:

numGhost := orderOfAccuracy/2; 
rNames[0]:="r": rNames[1]:="s": rNames[2]:="t":
xNames[0]:="x": xNames[1]:="y": xNames[2]:="z":

for ord from 2 by 2 to orderOfAccuracy-2 do
  if nd=2 then ord3:=0; else ord3:=ord; end if;

  fprintf(file,"numGhost1=%d;\n",numGhost-ord/2);
  fprintf(file,"n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);\n");
  fprintf(file,"n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);\n");
  fprintf(file,"n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);\n");

  fprintf(file,"beginLoops3d()\n");
  fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
  fprintf(file,"  if( mask(i1,i2,i3).ne.0 )then\n");
  fprintf(file,"  #End \n"):
  if ord=2 then
    # special case at order 2 : we difference u and later d
    lhsName:= saveDpDm( d,u,2,0,0 );
    declareVar(sprintf("d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)")):
    lhsName:=saveDpDm( d,u,0,2,0 );
    declareVar(sprintf("d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)")):
    if nd=3 then
    lhsName:=saveDpDm( d,u,0,0,2 );
    declareVar(sprintf("d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)")):
    end if;

  else
    for d3 from 0 by 2 to ord3 do
    for d2 from 0 by 2 to ord do
    for d1 from 0 by 2 to ord do
      if d1+d2+d3 = ord then
        saveDpDm( d,d,d1,d2,d3 );
        declareVar(sprintf("d%d%d%d(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)",d1,d2,d3)):
      end if;
    end do;
    end do;
    end do;
  end if:
  if gridType=curvilinear then
    fprintf(file,"    do m1=0,%d\n",nd-1);
    fprintf(file,"    do m2=0,%d\n",nd-1);
    # fprintf(file,"    do m2=m1,nd-1\n");
    if ord=2 then
      # special case at order 2 : we difference u and later d
      saveDpDm( rx,rsxy,2,0,0 );
      declareVar(sprintf("rx200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:%d,0:%d)",nd-1,nd-1)):
      saveDpDm( rx,rsxy,0,2,0 );
      declareVar(sprintf("rx020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:%d,0:%d)",nd-1,nd-1)):
      if nd=3 then
      saveDpDm( rx,rsxy,0,0,2 ); 
      declareVar(sprintf("rx002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:%d,0:%d)",nd-1,nd-1)):
      end if:     
    else
        for d3 from 0 by 2 to ord3 do
        for d2 from 0 by 2 to ord do
        for d1 from 0 by 2 to ord do
          if d1+d2+d3 = ord then
            saveDpDm( rx,rx,d1,d2,d3 );
            declareVar(sprintf("rx%d%d%d(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:%d,0:%d)",d1,d2,d3,nd-1,nd-1)):
          end if;
        end do;
        end do;
        end do;      
    end if:
    fprintf(file,"    end do\n");
    fprintf(file,"    end do\n");  
  end if:

  fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
  fprintf(file,"  end if ! mask .ne. 0\n");
  fprintf(file,"  #End \n"):
  fprintf(file,"endLoops3d() \n");
  
end do:

fprintf(file,"\n! ------------ MAIN LOOP, %dD, ORDER %d ----------------\n",nd,orderOfAccuracy);

ord:=orderOfAccuracy;
if nd=2 then ord3:=0; else ord3:=ord; end if:
fprintf(file,"numGhost1=%d; \n",numGhost-ord/2);
fprintf(file,"n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);\n");
fprintf(file,"n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);\n");
fprintf(file,"n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);\n");


fprintf(file,"\n! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----\n");
fprintf(file,"! Example: \n");
if gridType=rectangular then
  fprintf(file,"! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...\n");
else
  fprintf(file,"! u.rr = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...\n");
end if:
maxDeriv := orderOfAccuracy:
numTerms := orderOfAccuracy/2-1:

cx[0]:="": cx[1]:="": cx[2]:="":
for d1 from 1 by 1 to maxDeriv do
  if gridType=rectangular then
    cx[0] := cat(cx[0],"x");   cx[1] := cat(cx[1],"y");  cx[2] := cat(cx[2],"z");
  else
    cx[0] := cat(cx[0],"r");   cx[1] := cat(cx[1],"s");  cx[2] := cat(cx[2],"t");
  end if:
  for ax from 0 to nd-1 do
    fprintf(file,"c%s0 = 1./(dx(%d)**%d); ",cx[ax],ax,d1);
    declareVar(sprintf("c%s0",cx[ax])):
    # fprintf(dfile,"real c%s0\n",cx[ax]);
    for i from 1 to numTerms do
      fprintf(file,"c%s%d = (%s.); ",cx[ax],i,convert(dw[i,d1],string));
      # fprintf(file,"c%s%d = (%s.)*dx(%d)**%d; ",cx[ax],i,convert(dw[i,d1],string),ax,2*i);
      declareVar(sprintf("c%s%d",cx[ax],i)):
      # fprintf(dfile,"real c%s%d\n",cx[ax],i);
    end do;
    fprintf(file,"\n"):
  end do:
end do:


fprintf(file,"\n");
fprintf(file,"fv(m)=0.\n");
fprintf(file,"beginLoops3d()\n");
fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
fprintf(file,"  if( mask(i1,i2,i3).ne.0 )then\n");
fprintf(file,"  #End \n"):
fprintf(file,"\n    ! -- Note: Highest differences do not need to be stored in an array \n");

ord :=orderOfAccuracy:

#  **NEW WAY : do even and odd derivatives in one loop

fprintf(file,"\n    ! ----- DIFFERENCES OF U -----\n");
defineDifferences( d,u ):

if gridType=curvilinear then
  fprintf(file,"\n    ! ----- DIFFERENCES OF METRICS -----\n");
  fprintf(file,"    do m1=0,%d\n",nd-1);
  fprintf(file,"    do m2=0,%d\n",nd-1);

  # fprintf(file,"    do m2=m1,nd-1\n");
  fprintf(file,"    rx000i(m1,m2) = rsxy(i1,i2,i3,m1,m2)\n");
  declareVar("rx000i(0:2,0:2)"):

  defineDifferences( rx,rsxy ):

  fprintf(file,"    end do\n");
  fprintf(file,"    end do\n");  

end if:
fprintf(file,"\n\n");


if gridType=rectangular then
  # ------------ RECTANGULAR DERIVATIVES ------------------
  for ord from 2 by 2 to orderOfAccuracy do
    if nd=2 then ord3:=0; else ord3:=ord; end if;

    fprintf(file,"\n    ! ----- spatial derivatives of order %d ----\n",ord);
    for d3 from 0 by 2 to ord3 do
    for d2 from 0 by 2 to ord do
    for d1 from 0 by 2 to ord do
      if d1+d2+d3 = ord then
        # saveDpDm( d,d,d1,d2,d3,shortTargetName );
        suffix:=""; xn:=""; yn:=""; zn:=""; 
        for i1 from 2 by 2 to d1 do xn := cat( xn, "xx"); suffix := cat( suffix, "xx"); end do;
        for i2 from 2 by 2 to d2 do yn := cat( yn, "yy"); suffix := cat( suffix, "yy"); end do;
        for i3 from 2 by 2 to d3 do zn := cat( zn, "zz"); suffix := cat( suffix, "zz"); end do;
        derivName := cat("u",suffix); 
        # fprintf(file,"%s = d%d%d%d(i1,i2,i3,0) + finish me\n",derivName,d1,d2,d3);
        fprintf(file,"    %s =",derivName,d1,d2,d3);
        declareVar(derivName):

        if d1>0 then d1Max:=orderOfAccuracy; else d1Max:=0; end if;
        if d2>0 then d2Max:=orderOfAccuracy; else d2Max:=0; end if;
        if d3>0 then d3Max:=orderOfAccuracy; else d3Max:=0; end if;
        for i3 from d3 by 2 to d3Max do
        for i2 from d2 by 2 to d2Max do
        for i1 from d1 by 2 to d1Max do
          if i1+i2+i3 <= orderOfAccuracy then
            if i1=d1 and i2=d2 and i3=d3 then
              if d1>0 then
                fprintf(file,"c%s%d*",xn,0);    # e.g. cxx0
              end if:
              if d2>0 then
                fprintf(file,"c%s%d*",yn,0); # e.g. cyy0
              end if:
              if d3>0 then
                fprintf(file,"c%s%d*",zn,0); # e.g. czz0
              end if:
              fprintf(file,"(",zn,0); # e.g. czz0
            else
              fprintf(file," + "); 
              if i1-d1>0 then
                fprintf(file,"c%s%d*",xn,(i1-d1)/2);    # e.g. cxx1, cxx2, cxxxx1, ...
              end if:
              if i2-d2>0 then
                fprintf(file,"c%s%d*",yn,(i2-d2)/2); # e.g. cyy1, cyy2, cyyyy1, ...
              end if:
              if i3-d3>0 then
                fprintf(file,"c%s%d*",zn,(i3-d3)/2); # e.g. czz1, czz2, czzzz1, ...
              end if:
            end if;          
            fprintf(file,"d%d%d%di",i1,i2,i3);   # e.g. d240i , d202i, ...
          end if;
        end do;
        end do;
        end do;
        fprintf(file,")\n",i1,d2,d3);




      end if;
    end do:
    end do:
    end do:  

  end do:




else

  # ------------ CURVILINEAR DERIVATIVES ------------------

  for ord from 1 by 1 to orderOfAccuracy do     # order of derivative 

    if nd=2 then ord3:=0; else ord3:=ord; end if;

    fprintf(file,"\n    ! ----- parametric derivatives of u of order %d ----\n",ord);
    # Example: 
    # ! ----- parametric derivatives of u of order 2 ----
    # urr =d200i + crr1*d400i + crr2*d600i + crr3*d800i
    # urs =d110i + cr1*d310i + cr2*d510i + cr3*d710i + cs1*d130i + cr1*cs1*d330i + cr2*cs1*d530i + cs2*d150i + cr1*cs2*d350i + cs3*d170i
    # uss =d020i + css1*d040i + css2*d060i + css3*d080i 

    for d3 from 0 by 1 to ord3 do
    for d2 from 0 by 1 to ord do
    for d1 from 0 by 1 to ord do
      if d1+d2+d3 = ord then
        # saveDpDm( d,d,d1,d2,d3,shortTargetName );
        suffix:=""; xn:=""; yn:=""; zn:=""; rn:=""; sn:=""; tn:="": 
        for i1 from 1 by 1 to d1 do xn := cat( xn, "x"); rn := cat( rn, "r"); suffix := cat( suffix, "x"); end do;
        for i2 from 1 by 1 to d2 do yn := cat( yn, "y"); sn := cat( sn, "s"); suffix := cat( suffix, "y"); end do;
        for i3 from 1 by 1 to d3 do zn := cat( zn, "z"); tn := cat( tn, "t"); suffix := cat( suffix, "z"); end do;
        derivName := cat("u",suffix); 
        rDerivName := cat("u",cat(rn,cat(sn,tn))):
        fprintf(file,"    %s =",rDerivName,d1,d2,d3);
        declareVar(rDerivName):

        if d1>0 then d1Max:=orderOfAccuracy; else d1Max:=0; end if;
        if d2>0 then d2Max:=orderOfAccuracy; else d2Max:=0; end if;
        if d3>0 then d3Max:=orderOfAccuracy; else d3Max:=0; end if;
        for i3 from d3 by 2 to d3Max do
        for i2 from d2 by 2 to d2Max do
        for i1 from d1 by 2 to d1Max do
          if i1+i2+i3 <= orderOfAccuracy then
            if i1=d1 and i2=d2 and i3=d3 then

            else
              fprintf(file," + "); 
              if i1-d1>0 then
                fprintf(file,"c%s%d*",rn,(i1-d1)/2);    # e.g. cxx1, cxx2, cxxxx1, ...
              end if:
              if i2-d2>0 then
                fprintf(file,"c%s%d*",sn,(i2-d2)/2); # e.g. cyy1, cyy2, cyyyy1, ...
              end if:
              if i3-d3>0 then
                fprintf(file,"c%s%d*",tn,(i3-d3)/2); # e.g. czz1, czz2, czzzz1, ...
              end if:
            end if;          
            fprintf(file,"d%d%d%di",i1,i2,i3);   # e.g. d240i , d202i, ...
          end if;
        end do;
        end do;
        end do;
        fprintf(file,"\n",i1,d2,d3);

      end if;
    end do:
    end do:
    end do:  

    fprintf(file,"\n    !-- get coefficients in spatial derivatives of order %d --\n",ord);
    if ord<orderOfAccuracy then
      fprintf(file,"    getDerivCoeff%dd(%d)\n",nd,ord);
    end if:

    if ord=1 then
      # declare
      #   rx = rsxy(i1,i2,i3,0,0)
      #   sx = rsxy(i1,i2,i3,1,0)
      # etc
      for i1 from 0 to nd-1 do
      for i2 from 0 to nd-1 do
        fprintf(file,"    %s%s = rsxy(i1,i2,i3,%d,%d)\n",rNames[i1],xNames[i2],i1,i2):
        declareVar(sprintf("%s%s",rNames[i1],xNames[i2])):
      end do:
    end do:
    end if:

    if ord<orderOfAccuracy then # no need to define these for the highest derivative
      fprintf(file,"\n    ! --- define parametric and spatial derivatives of the metrics of order %d (used by higher derivatives of u) ---\n",ord);
      fprintf(file,"    do m1=0,%d\n",nd-1);
      fprintf(file,"    do m2=0,%d\n",nd-1);
      # fprintf(file,"    do m2=m1,nd-1\n");
      for d3 from 0 by 1 to ord3 do
      for d2 from 0 by 1 to ord do
      for d1 from 0 by 1 to ord do
        if d1+d2+d3 = ord then
          rn:=""; sn:=""; tn:="": 
          for i1 from 1 by 1 to d1 do  rn := cat( rn, "r");  end do;
          for i2 from 1 by 1 to d2 do  sn := cat( sn, "s");  end do;
          for i3 from 1 by 1 to d3 do  tn := cat( tn, "t");  end do;
          rDerivName := cat("rx",cat(rn,cat(sn,tn))):
          fprintf(file,"      %sa(m1,m2) =",rDerivName);
          declareVar(sprintf("%sa(0:2,0:2)",rDerivName)):

          if d1>0 then d1Max:=orderOfAccuracy; else d1Max:=0; end if;
          if d2>0 then d2Max:=orderOfAccuracy; else d2Max:=0; end if;
          if d3>0 then d3Max:=orderOfAccuracy; else d3Max:=0; end if;
          for i3 from d3 by 2 to d3Max do
          for i2 from d2 by 2 to d2Max do
          for i1 from d1 by 2 to d1Max do
            if i1+i2+i3 <= orderOfAccuracy then
              if i1=d1 and i2=d2 and i3=d3 then

              else
                fprintf(file," + "); 
                if i1-d1>0 then
                  fprintf(file,"c%s%d*",rn,(i1-d1)/2);    # e.g. crr1, crr, ...
                end if:
                if i2-d2>0 then
                  fprintf(file,"c%s%d*",sn,(i2-d2)/2); # e.g. cs1, css2, ...
                end if:
                if i3-d3>0 then
                  fprintf(file,"c%s%d*",tn,(i3-d3)/2); # e.g. ct1, ctt2, ...
                end if:
              end if;          
              fprintf(file,"rx%d%d%di(m1,m2)",i1,i2,i3);   # e.g. rx240i(m1,m2) , rx202i(m1,m2), ...
              # fprintf(file,"rx%d%d%d(i1,i2,i3,m1,m2)",i1,i2,i3);   # e.g. rx240 , rx202, ...
            end if;
          end do;
          end do;
          end do;
          fprintf(file,"\n",i1,d2,d3);

        end if;
      end do:
      end do:
      end do: 

      # Example: 
      # rxxxa(m1,m2) = cuxx100*rxr(m1,m2,0)+cuxx200*rxrra(m1,m2)+cuxx010*rxr(m1,m2,1)+cuxx110*rxrsa(m1,m2)+cuxx020*rxssa(m1,m2)
      for d3 from 0 by 1 to ord3 do
      for d2 from 0 by 1 to ord do
      for d1 from 0 by 1 to ord do 
        if d1+d2+d3 = ord then

          suffix:=""; xn:=""; yn:=""; zn:=""; 
          for i1 from 1 by 1 to d1 do xn := cat( xn, "x"); end do;
          for i2 from 1 by 1 to d2 do yn := cat( yn, "y"); end do;
          for i3 from 1 by 1 to d3 do zn := cat( zn, "z"); end do;
          derivName := cat("rx",cat(xn,cat(yn,zn)));
          coeffName := cat("u",cat(xn,cat(yn,zn)));

          fprintf(file,"      %sa(m1,m2) =",derivName);  # e.g. ux=, uy=, uxx=, uxy=, ...
          declareVar(sprintf("%sa(0:2,0:2)",derivName)):
          for i3 from 0 by 1 to ord3 do
          for i2 from 0 by 1 to ord do
          for i1 from 0 by 1 to ord do
            if i1+i2+i3>0 and i1+i2+i3 <= ord then
             rn:=""; sn:=""; tn:="": 
             for j1 from 1 by 1 to i1 do rn := cat( rn, "r"); end do;
             for j2 from 1 by 1 to i2 do sn := cat( sn, "s"); end do;
             for j3 from 1 by 1 to i3 do tn := cat( tn, "t"); end do;

             fprintf(file," +c%s%d%d%d*rx%s%s%sa(m1,m2)",coeffName,i1,i2,i3,rn,sn,tn);
            end if:
          end do;
          end do;
          end do;     

          fprintf(file,"\n");        
        end if:
      end do:
      end do:
      end do: 

      fprintf(file,"    end do\n");
      fprintf(file,"    end do\n");
      fprintf(file,"    ! --- define short form names metric derivatives ---\n");

      for d3 from 0 by 1 to ord3 do
      for d2 from 0 by 1 to ord do
      for d1 from 0 by 1 to ord do 
        if d1+d2+d3 = ord then

          suffix:=""; xn:=""; yn:=""; zn:=""; 
          for i1 from 1 by 1 to d1 do xn := cat( xn, "x"); end do;
          for i2 from 1 by 1 to d2 do yn := cat( yn, "y"); end do;
          for i3 from 1 by 1 to d3 do zn := cat( zn, "z"); end do;

          derivName := cat("rx",cat(xn,cat(yn,zn)));

          fprintf(file,"   ");  
          for j1 from 0 to nd-1 do
          for j2 from 0 to nd-1 do
            rxName := cat(rNames[j1],xNames[j2]); 
            fprintf(file," %s%s%s%s=%sa(%d,%d);",rxName,xn,yn,zn, derivName,j1,j2);
            declareVar(sprintf(" %s%s%s%s",rxName,xn,yn,zn)):
          end do:
          end do:
          fprintf(file,"\n");
        end if:
      end do:
      end do:
      end do: 
    end if:

    #  do m1=0,nd-1
    #  do m2=0,nd-1
    #    ! ---- 2nd parameteric derivatives of the metrics ----
    #    rxrra(m1,m2) = rx200(i1,i2,i3,m1,m2) - (dx(0)**2/12.)*rx400(i1,i2,i3,m1,m2) + (dx(0)**4/90.)*rx600(i1,i2,i3,m1,m2) - (dx(0)**6/560.)*rx800i(m1,m2)
    #    rxssa(m1,m2) = rx020(i1,i2,i3,m1,m2) - (dx(1)**2/12.)*rx040(i1,i2,i3,m1,m2) + (dx(1)**4/90.)*rx060(i1,i2,i3,m1,m2) - (dx(1)**6/560.)*rx080i(m1,m2) 
    #    rxrsa(m1,m2) = rx110i(m1,m2) - (dx(0)**2/6.)*rx310i(m1,m2) - (dx(1)**2/6.)*rx130i(m1,m2) + (dx(0)**4/30.)*rx510i(m1,m2) + (dx(1)**4/30.)*rx150i(m1,m2) + (dx(0)**2 * dx(1)**2/36. )*rx330i(m1,m2)      \
    #          - ( (1./140.)*dx(0)**6 )*rx710i(m1,m2) - ( (1./140.)*dx(1)**6 )*rx170i(m1,m2) - (dx(0)**4*dx(1)**2/(6.*30.))*rx530i(m1,m2) - (dx(0)**2*dx(1)**4/(6.*30.))*rx350i(m1,m2)   
    #  
    #    ! ---- 2nd spatial derivatives of the metrics ----
    #    rxxxa(m1,m2) = cuxx100*rxr(m1,m2,0)+cuxx200*rxrra(m1,m2)+cuxx010*rxr(m1,m2,1)+cuxx110*rxrsa(m1,m2)+cuxx020*rxssa(m1,m2)
    #    rxxya(m1,m2) = cuxy100*rxr(m1,m2,0)+cuxy200*rxrra(m1,m2)+cuxy010*rxr(m1,m2,1)+cuxy110*rxrsa(m1,m2)+cuxy020*rxssa(m1,m2)
    #    rxyya(m1,m2) = cuyy100*rxr(m1,m2,0)+cuyy200*rxrra(m1,m2)+cuyy010*rxr(m1,m2,1)+cuyy110*rxrsa(m1,m2)+cuyy020*rxssa(m1,m2)
    #  end do
    #  end do
    #  rxxx=rxxxa(0,0); ryxx=rxxxa(0,1); sxxx=rxxxa(1,0); syxx=rxxxa(1,1); 
    #  rxxy=rxxya(0,0); ryxy=rxxya(0,1); sxxy=rxxya(1,0); syxy=rxxya(1,1); 
    #  rxyy=rxyya(0,0); ryyy=rxyya(0,1); sxyy=rxyya(1,0); syyy=rxyya(1,1); 

    if ord mod 2 = 0 then
      fprintf(file,"\n    ! (even) spatial derivatives of u\n");
      # uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
      for d3 from 0 by 1 to ord3 do
      for d2 from 0 by 1 to ord do
      for d1 from 0 by 1 to ord do 
        if d1+d2+d3 = ord and (d1 mod 2 = 0) and (d2 mod 2 = 0) and (d3 mod 2 = 0) then

          suffix:=""; xn:=""; yn:=""; zn:=""; 
          for i1 from 1 by 1 to d1 do xn := cat( xn, "x"); end do;
          for i2 from 1 by 1 to d2 do yn := cat( yn, "y"); end do;
          for i3 from 1 by 1 to d3 do zn := cat( zn, "z"); end do;
          derivName := cat("u",cat(xn,cat(yn,zn)));

          # eval selected coeff in chain rule at highest order only: 
          if ord=orderOfAccuracy then itStart:=0: else itStart:=1: end if: 
          for it from itStart to 1 do

            if it=1 then 
              fprintf(file,"    %s =",derivName);   # e.g. ux=, uy=, uxx=, uxy=, ...
              declareVar(derivName):
            end if:
            
            for i3 from 0 by 1 to ord3 do
            for i2 from 0 by 1 to ord do
            for i1 from 0 by 1 to ord do
              if i1+i2+i3>0 and i1+i2+i3 <= ord then
               rn:=""; sn:=""; tn:="": 
               for j1 from 1 by 1 to i1 do rn := cat( rn, "r"); end do;
               for j2 from 1 by 1 to i2 do sn := cat( sn, "s"); end do;
               for j3 from 1 by 1 to i3 do tn := cat( tn, "t"); end do;

               derivCoeffName:= sprintf("c%s%d%d%d",derivName,i1,i2,i3);
               derivCoeffVar:= cat(c,cat(derivName,cat(i1,cat(i2,i3))));

               if it=0 then fprintf(file,"    %s=%s\n",derivCoeffName,convert(derivCoeffVar,string)); end if:
               if it=1 then fprintf(file," +c%s%d%d%d*u%s%s%s",derivName,i1,i2,i3,rn,sn,tn); end if:
              end if:
            end do;
            end do;
            end do;     
          end do: # for it 
          fprintf(file,"\n");        
        end if:
      end do:
      end do:
      end do: 

    end if:

    #  uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
    #  ! uxxxy = cuxxxy100*ur+cuxxxy200*urr+cuxxxy300*urrr+cuxxxy400*urrrr+cuxxxy010*us+cuxxxy110*urs+cuxxxy210*urrs+cuxxxy310*urrrs+cuxxxy020*uss+cuxxxy120*urss+cuxxxy220*urrss+cuxxxy030*usss+cuxxxy130*ursss+cuxxxy040*ussss
    #  uxxyy = cuxxyy100*ur+cuxxyy200*urr+cuxxyy300*urrr+cuxxyy400*urrrr+cuxxyy010*us+cuxxyy110*urs+cuxxyy210*urrs+cuxxyy310*urrrs+cuxxyy020*uss+cuxxyy120*urss+cuxxyy220*urrss+cuxxyy030*usss+cuxxyy130*ursss+cuxxyy040*ussss
    #  ! uxyyy = cuxyyy100*ur+cuxyyy200*urr+cuxyyy300*urrr+cuxyyy400*urrrr+cuxyyy010*us+cuxyyy110*urs+cuxyyy210*urrs+cuxyyy310*urrrs+cuxyyy020*uss+cuxyyy120*urss+cuxyyy220*urrss+cuxyyy030*usss+cuxyyy130*ursss+cuxyyy040*ussss
    #  uyyyy = cuyyyy100*ur+cuyyyy200*urr+cuyyyy300*urrr+cuyyyy400*urrrr+cuyyyy010*us+cuyyyy110*urs+cuyyyy210*urrs+cuyyyy310*urrrs+cuyyyy020*uss+cuyyyy120*urss+cuyyyy220*urrss+cuyyyy030*usss+cuyyyy130*ursss+cuyyyy040*ussss


  end do: # end order of derivative

end if:

fprintf(file,"\n");
fprintf(file,"  #If #FORCING ne \"NOFORCING\" \n"):
fprintf(file,"    getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) \n");
fprintf(file,"  #End \n\n"):

fprintf(file,"    ! ---- Modified Equation UPDATE ----- \n");
if nd=2 then
  fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \\\n");
  fprintf(file,"                     + cdtsq*( uxx + uyy ) \\\n");
  if orderOfAccuracy>2 then
  fprintf(file,"                     + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  \\\n");
  end if:
  if orderOfAccuracy>4 then
  fprintf(file,"                     + cdtPow6By360*( uxxxxxx + uyyyyyy\\\n");
  fprintf(file,"                               + 3.*( uxxxxyy + uxxyyyy) ) \\\n");
  end if:
  if orderOfAccuracy>6 then
  fprintf(file,"                     + cdtPow8By20160*( uxxxxxxxx + uyyyyyyyy \\\n");
  fprintf(file,"                                 + 4.*( uxxxxxxyy + uxxyyyyyy )\\\n");
  fprintf(file,"                                 + 6.*( uxxxxyyyy ) ) \\\n");
  end if:
  # fprintf(file,"                     + dtSq*fv(m)  \n");  
  fprintf(file,"     FV(m)  \n"): # optionally add forcing   
  fprintf(file,"\n");
else
  fprintf(file,"    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \\\n");
  fprintf(file,"                     + cdtsq*( uxx + uyy +uzz ) \\\n");
  if orderOfAccuracy>2 then
  fprintf(file,"                     + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*(uxxyy+uxxzz+uyyzz) )  \\\n");
 end if:
  if orderOfAccuracy>4 then
  fprintf(file,"                     + cdtPow6By360*( uxxxxxx + uyyyyyy + uzzzzzz \\\n");
  fprintf(file,"                                + 3.*( uxxxxyy + uxxyyyy + uxxxxzz + uxxzzzz + uyyyyzz + uyyzzzz ) \\\n");
  fprintf(file,"                                + 6.*(  uxxyyzz  ) ) \\\n");
  end if:
  if orderOfAccuracy>6 then  
  fprintf(file,"                     + cdtPow8By20160*( uxxxxxxxx + uyyyyyyyy + uzzzzzzzz \\\n");
  fprintf(file,"                                 + 4.*( uxxxxxxyy + uxxyyyyyy + uxxxxxxzz + uxxzzzzzz + uyyyyyyzz + uyyzzzzzz ) \\\n");
  fprintf(file,"                                 + 6.*( uxxxxyyyy + uxxxxzzzz + uyyyyzzzz ) \\\n");
  fprintf(file,"                                 +12.*( uxxxxyyzz + uxxyyyyzz + uxxyyzzzz ) )\\\n");
  end if:
  # fprintf(file,"                     + dtSq*fv(m)  \n");  
  fprintf(file,"                       FV(m)  \n"): # optionally add forcing 
  fprintf(file,"\n");

end if;


fprintf(file,"  #If #MASK eq \"USEMASK\" \n"):
fprintf(file,"  end if ! mask .ne. 0\n");
fprintf(file,"  #End \n"):
fprintf(file,"endLoops3d() \n");


  fprintf(file,"#endMacro"); 
  declareVar(finishDeclarations); # finish declarations
  fprintf(dfile,"#endMacro"); 
  fclose(file);
  fclose(dfile);
  printf("Wrote file=%s\n", myFileName);
  printf("Wrote file=%s\n", declareFileName);
  

end do:
end do:
# -----------------------------------------------------------------
# -------------- END LOOP OVER ORDER AND DIMENSIONS ---------------
# -----------------------------------------------------------------
