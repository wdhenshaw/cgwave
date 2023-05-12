#  
# Generate stencil include files for modified equation time stepping
#
# restart: currentdir("/Users/henshaw/Dropbox/research/cgwave/maple"): currentdir("/Users/henshaw/Dropbox/research/cgwave/maple"): read "writeStencilFiles.mpl";
#
# with(LinearAlgebra):
# with(plottools):
# with(plots):
with(StringTools):
with(CodeGeneration):

# kernelopts(printbytes=false): # turn off memory used messages
# kernelopts('printlevel'=0): 
# interface(echo=0): 
# interface(quiet=true): 
interface( quiet=true,echo=0 ):

rectangular:=1: curvilinear:=2: # enumerators 
gridTypeName[rectangular]:="Rectangular"; 
gridTypeName[curvilinear]:="Curvilinear"; 

ndStart := 2: ndEnd:=2:     # number of dimensions
orderStart:=2: orderEnd:=2: # order of accuracy
# orderOfAccuracy  := 8:      # order of accuracy
# gridType         := rectangular: 
# gridType         := curvilinear;

gridTypeStart:=1: gridTypeEnd:=2: # gridType 

includeDir := "../include"; # put include files here


# ********* SET CASE HERE ********
ndStart:=3; icase:=8; 

if icase=1 then
  nd:=ndStart; orderOfAccuracy:=2; gridType:=rectangular;
elif icase=2 then
  nd:=ndStart; orderOfAccuracy:=2; gridType:=curvilinear;
elif icase=3 then
  nd:=ndStart; orderOfAccuracy:=4; gridType:=rectangular;
elif icase=4 then
  nd:=ndStart; orderOfAccuracy:=4; gridType:=curvilinear;     
elif icase=5 then
  nd:=ndStart; orderOfAccuracy:=6; gridType:=rectangular;
elif icase=6 then
  nd:=ndStart; orderOfAccuracy:=6; gridType:=curvilinear; 
elif icase=7 then
  nd:=ndStart; orderOfAccuracy:=8; gridType:=rectangular;
elif icase=8 then
  nd:=ndStart; orderOfAccuracy:=8; gridType:=curvilinear;  
else
end if: 



# files came from AMP
#   Dropbox/AMP/fame/codes/fame_stencilCoefGenerators/stencil_2D/generatedStencilCoefs_160423/
#   Dropbox/AMP/fame/codes/fame_stencilCoefGenerators/stencil_3D/generatedStencilCoefs_160423> 
#


fileName := sprintf("stencil%dD_%d.txt",nd,orderOfAccuracy);
printf("Read file=%s\n", fileName);



stencilWidth:=orderOfAccuracy+1;
hw := orderOfAccuracy/2;          # stencil half-width
if nd=3 then
  hw3:=hw;
else
  hw3:=0;
end if: 
numStencilCoeff := stencilWidth^nd;

extra := orderOfAccuracy/2;

printf("OPEN to read file=%s\n", fileName);
file := fopen(fileName, READ):
# numLines:=5; 
numLines:= numStencilCoeff;

sc := array(1..numLines):
# b := array(1..numLines,1..3,1..3):
# c10 := array(1..3,1..3);
# c01 := array(1..3,1..3);
# c20 := array(1..3,1..3);
# c11 := array(1..3,1..3);
# c02 := array(1..3,1..3);
cpu0 := time();
j:=1; numZero:=0; 
for i from 1 to numLines+extra do
  str := readline(file);
  # printf("str=[%s]\n",str);
  if str=0 then
    break;
  end if;

  if str[1]="b" then

    str := SubstituteAll(str,"=",":="):
    str := SubstituteAll(str,"I","i1"):
    str := SubstituteAll(str,"J","i2"):
    str := SubstituteAll(str,"K","i3"):

    str := SubstituteAll(str,"i1-1","i1m1"):
    str := SubstituteAll(str,"i1+1","i1p1"):
    str := SubstituteAll(str,"i1-2","i1m2"):
    str := SubstituteAll(str,"i1+2","i1p2"): 
    str := SubstituteAll(str,"i1-3","i1m3"):
    str := SubstituteAll(str,"i1+3","i1p3"):     

    str := SubstituteAll(str,"i2-1","i2m1"):
    str := SubstituteAll(str,"i2+1","i2p1"):
    str := SubstituteAll(str,"i2-2","i2m2"):
    str := SubstituteAll(str,"i2+2","i2p2"): 
    str := SubstituteAll(str,"i2-3","i2m3"):
    str := SubstituteAll(str,"i2+3","i2p3"):     

    if nd=3 then
      str := SubstituteAll(str,"i3-1","i3m1"):
      str := SubstituteAll(str,"i3+1","i3p1"):
      str := SubstituteAll(str,"i3-2","i3m2"):
      str := SubstituteAll(str,"i3+2","i3p2"): 
      str := SubstituteAll(str,"i3-3","i3m3"):
      str := SubstituteAll(str,"i3+3","i3p3"):       
    end if:       

    # printf("str=[%s]\n",str);

    sc[j] := parse(str);
    if sc[j]=0 then
      numZero:=numZero+1;
      printf("** read zero coeff: sc[%d]=0\n",j);
    end if;

    j:=j+1;

  end if;
end do:
if numZero>0 then
  printf("$$$ there were %d zero coefficients out of %d $$$\n",numZero,numStencilCoeff);
end if;


fclose(file);
cpu := time()-cpu0;
printf("CLOSE file=%s, time to read=%9.2e (s)\n",fileName,cpu);




  if gridType=rectangular then

    # --- RECTANGULAR (CARTESIAN) GRID VERSION ----

    printf("RECTANGULAR (CARTESIAN) GRID VERSION...\n");
    if nd=2 then
      c20 := (i1,i2,i3) -> cdx2i;  # c^2/dx^2
      c02 := (i1,i2,i3) -> cdy2i; 
      c11 := (i1,i2,i3) -> 0; 
      c10 := (i1,i2,i3) -> 0; 
      c01 := (i1,i2,i3) -> 0;
    else
      c200 := (i1,i2,i3) -> cdx2i;
      c020 := (i1,i2,i3) -> cdy2i; 
      c002 := (i1,i2,i3) -> cdz2i; 
      c110 := (i1,i2,i3) -> 0; 
      c101 := (i1,i2,i3) -> 0; 
      c011 := (i1,i2,i3) -> 0; 
      c100 := (i1,i2,i3) -> 0; 
      c010 := (i1,i2,i3) -> 0;      
      c001 := (i1,i2,i3) -> 0;      
    end if: 
    # printf(".. CALL Matlab codegen...\n");
    # printf(".. sc[1]=%s\n",convert(sc[1],string));
    myCode := Matlab( sc,output=string ,optimize,declare=[S::numeric,i1::integer,i2::integer]):
    # printf(".. Done Matlab codegen...\n");

    # results:
    #   scr(1) = 
    #   scr(2) = 
    #   scr(3) = 
    myCode := RegSubs("sc\\(([0-9]*)\\)" = "scr(\\1)", myCode):

    myCode := RegSubs("\\^" = "**", myCode):  # power to **

    # printf("myCode=%s\n",myCode);



  else
    printf("CURVILINEAR GRID VERSION, call CodeGen...\n");

    cpu0 := time();
    myCode := Matlab( sc,output=string ,optimize,declare=[S::numeric,i1::integer,i2::integer]):
    cpu := time()-cpu0;
    printf("... time for CodeGen = %9.2e (s)\n",cpu);

    if nd=2 then
      myCode := RegSubs("sc\\(([0-9]*)\\)" = "sc(\\1,i1,i2)", myCode):
    else
      myCode := RegSubs("sc\\(([0-9]*)\\)" = "sc(\\1,i1,i2,i3)", myCode):
    end if:
    # myCode;

  end if:


  # myCode := RegSubs("sc\(1\)"="SC\(1,i2,i3\)",myCode ):
  # myCode := RegSubs("sc\(([1..9]*)"="SC\(\\1",myCode ):
  # myCode := RegSubs( "sc\((.*)\)"="sc(\\1,i2,i3)" , myCode );

  # myCode := Fortran( S ,optimize,declare=[S::numeric,i1::integer,i2::integer]):
  # myCode := Fortran( S ,optimize,declare=[S::float,c20::float,c11::float,c20::float,c01::float,c10::float,i1::integer,i2::integer]):
  # myCode := fortran( S ,filename="prog.f",optimize,declare=[S::float,c20::float,c11::float,c20::float,c01::float,c10::float,i1::integer,i2::integer]):
  # myCode := C( S ,output=string,optimize,declare=[S::float,c20::float,c11::float,c20::float,c01::float,c10::float,i1::integer,i2::integer]):
  # myCode := C([S[1],S[2],S[3]],output=string,optimize,declare=[c20::float,c11::float,c20::float,c01::float,c10::float,i1::integer,i2::integer]):
  # myCode := C([S[1],S[2]],output=string,resultname="XX",declare=[k::float,beta::float,A1::float,B1::float,A2::float,B2::float,A3::float,B3::float]):
  # myCode;

  # -- define file name: 
  myFileName := sprintf("%s/defineStencilVariables%ddOrder%d%s.h",includeDir,nd,orderOfAccuracy,gridTypeName[gridType]); 

  file := open(myFileName, WRITE);

  fprintf(file,"! Define variables to valuate stencil coefficients, dim=%d, order=%d, gridType=%s\n",nd,orderOfAccuracy,gridTypeName[gridType]);
  fprintf(file,"! File generated by cgWave/maple/writeStencilFiles.mpl\n");

  # declare temp variables

  fprintf(file,"integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3\n");
  fprintf(file,"integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3\n");
  if nd=3 then
    fprintf(file,"integer i3m3,i3m2,i3m1,i3p1,i3p2,i3p3\n");
  end if:

  # two-D order 8 has t20455
  printf("Write file with declarations of temporaries...\n");
  cpu0 := time();
  maxTemp:=50000;
  fprintf(file,"real t0"):
  j:=0;
  numPerLine:=30; 
  for i from 1 to maxTemp do
    if RegMatch( cat("t",i," ="), myCode ) then
      fprintf(file,",t%d",i);
      j:=j+1;
      if modp(j,numPerLine)=0 then
        fprintf(file,"\\\n    ");
      fi
    end if:
  end do;
  fprintf(file,";\n");

  fclose(file);
  cpu := time()-cpu0; 
  printf("Wrote file=%s, cpu=%9.2e(s)\n", myFileName,cpu);


  # -- define file name: 
  myFileName := sprintf("%s/evalStencilCoeff%ddOrder%d%s.h",includeDir,nd,orderOfAccuracy,gridTypeName[gridType]); 
  file := open(myFileName, WRITE);

  fprintf(file,"! Evaluate stencil coefficients, dim=%d, order=%d, gridType=%s\n",nd,orderOfAccuracy,gridTypeName[gridType]);
  fprintf(file,"! File generated by cgWave/maple/writeStencilFiles.mpl\n");

  if orderOfAccuracy>2 then
    fprintf(file,"i1m1=i1-1; i1p1=i1+1;\n");
    fprintf(file,"i2m1=i2-1; i2p1=i2+1;\n");
    if nd=3 then
      fprintf(file,"i3m1=i3-1; i3p1=i3+1;\n");
    end if:
  end if:
  if orderOfAccuracy>4 then
    fprintf(file,"i1m2=i1-2; i1p2=i1+2;\n");
    fprintf(file,"i2m2=i2-2; i2p2=i2+2;\n");
    if nd=3 then
      fprintf(file,"i3m2=i3-2; i3p2=i3+2;\n");
    end if:  
  end if:
  if orderOfAccuracy>6 then
    fprintf(file,"i1m3=i1-3; i1p3=i1+3;\n");
    fprintf(file,"i2m3=i2-3; i2p3=i2+3;\n");
    if nd=3 then
      fprintf(file,"i3m3=i3-3; i3p3=i3+3;\n");
    end if:  
  end if:

  fprintf(file,"%s",myCode);

  fclose(file);
  printf("Wrote file=%s\n", myFileName);



  # ---------- TEMP ----
  if 1=1 then

  # -- define file name: 
  myFileName := sprintf("%s/stencil%ddOrder%d%s.h",includeDir,nd,orderOfAccuracy,gridTypeName[gridType]); 
  # declareFileName := sprintf("%s/declare%ddOrder%d%s.h",includeDir,nd,orderOfAccuracy,gridTypeName[gridType]);  # for declarations of variables 


  file := open(myFileName, WRITE);

  scName := "sc";
  if gridType=rectangular then
    scName := "scr";
  end if:
  fprintf(file,"! Stencil: nd=%d, orderOfAccuracy=%d, gridType=%s\n",nd,orderOfAccuracy,gridTypeName[gridType]);
  fprintf(file,"un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\\\n");
  icoeff:=1; nonZeroCoeff:=0; 
  for i3 from -hw3 to hw3 do
  for i2 from -hw to hw do
    for i1 from -hw to hw do
     if sc[icoeff]=0 then
        # printf("** skip: sc[%d]=0\n",icoeff);
        if gridType=curvilinear then
          # print spaces to make formula look nice
          if nd=2 then
            fprintf(file,"                                  ",icoeff,i1,i2);
          else
            fprintf(file,"                                       ",icoeff,i1,i2,i3);
          end if: 
        else
          if nd=2 then
            fprintf(file,"                             ",icoeff,i1,i2);
          else
            fprintf(file,"                                ",icoeff,i1,i2,i3);
          end if:            
        end if;  

     else
       nonZeroCoeff:=nonZeroCoeff+1;
        if gridType=curvilinear then
          if nd=2 then
            fprintf(file," + sc(%3d,i1,i2)*u(i1%+d,i2%+d,i3,m)",icoeff,i1,i2);
          else
            fprintf(file," + sc(%3d,i1,i2,i3)*u(i1%+d,i2%+d,i3%+d,m)",icoeff,i1,i2,i3);
          end if:
        else
          if nd=2 then
            fprintf(file," + scr(%3d)*u(i1%+d,i2%+d,i3,m)",icoeff,i1,i2);
          else
            fprintf(file," + scr(%3d)*u(i1%+d,i2%+d,i3%+d,m)",icoeff,i1,i2,i3);
          end if:          
        end if;

    end if:
     icoeff:=icoeff+1;
    end do:
    fprintf(file,"\\\n"); # new line
  end do:
  end do:
  fprintf(file," FV(m)\n");

  printf("Done: numStencilCoeff=%d, nonZeroCoeff=%d, zeroCoeff=%d\n",numStencilCoeff,nonZeroCoeff,numStencilCoeff-nonZeroCoeff);


  fclose(file);

  printf("Wrote file=%s\n", myFileName);

  end if:

