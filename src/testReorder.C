//===============================================================================
//
//  Test reordering routines
//
//==============================================================================


#include "Overture.h"  
#include "display.h"



int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture


  int num = 7;
  RealArray a(num);
  IntegerArray ia(num);


  if( 1==1 )
  {
    a(0)=8.; a(1)=6.; a(2)=7.; a(3)=5.; a(4)=3.; a(5)=0.; a(6)=9.;
    ia(0)=3; ia(1)=6; ia(2)=2; ia(3)=4; ia(4)=0; ia(5)=1; ia(6)=5; 
  }
  else if( 1==0 )
  {
    // no changes needed:
    a(0)=8.; a(1)=6.; a(2)=7.; a(3)=5.; a(4)=3.; a(5)=0.; a(6)=9.;
    ia(0)=0; ia(1)=1; ia(2)=2; ia(3)=3; ia(4)=4; ia(5)=5; ia(6)=6;    
  }
  else if( 1==0 )
  {
    // one swap needed
    a(0)=8.; a(1)=6.; a(2)=7.; a(3)=5.; a(4)=3.; a(5)=0.; a(6)=9.;
    ia(0)=0; ia(1)=2; ia(2)=1; ia(3)=3; ia(4)=4; ia(5)=5; ia(6)=6;    
  }
  else 
  {
    // reverse order
    a(0)=8.; a(1)=6.; a(2)=7.; a(3)=5.; a(4)=3.; a(5)=0.; a(6)=9.;
    ia(0)=6; ia(1)=5; ia(2)=4; ia(3)=3; ia(4)=2; ia(5)=1; ia(6)=0;    
  }  


   ::display(a,"a at start");
   ::display(ia,"ia at start");

    RealArray b(num);
    b = a(ia);
    ::display(b," b = a(ia) : reordered");

   int numAssigns=0; 
   for( int i=0; i<num; i++ )
   {
      bool entrySaved=false; 
      Real x=-1.;
      // delay this copy:
      // Real x = a(i); numAssigns++;
      // printF("step: i=%d : save x= a[i=%d]\n",i,i);
      int j=i;
      for( int ii=0; ii<num; ii++ )
      {
        int k = ia(j);
        ia(j)=j;
        if( k==i )
          break;

        if( j!=k )
        {
          if( !entrySaved )
          {
            x = a(i); numAssigns++; entrySaved=true;
            printF("step: i=%d : save x= a[i=%d]\n",i,i);
          }
          a(j)=a(k); numAssigns++;
          printF("step: i=%d : set a[j=%d] = a[k=%d]\n",i,j,k);
        }
        j=k;
      }
      if( i!=j )
      {
        assert( entrySaved );
        a(j)=x; numAssigns++;
        printF("step: i=%d : set a[j=%d] = a[i=%d] = x\n",i,j,i);
      }
   }

   printF("numAssigns=%d\n",numAssigns);
   ::display(a,"a at end");
   ::display(ia,"ia at end");

   Real maxErr = max(fabs(a-b));
   printf("maxErr = %8.2e\n",maxErr);


   Overture::finish();          

   return(0);
}

