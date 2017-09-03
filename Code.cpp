// Advanced Numerical Analysis, Mech. Eng. Dept., Sharif University of Technology
// Developed by A. F. Forughi (Dec. 2010)";
// This program is able to read the Harwell-Boeing matrix [A][x]=[B] (non-zero values of [A] and [B]), saved as "input.txt";
// The BiCGSTAB algorithm solves the linear system in order to find [x];

# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string.h>
# include <conio.h>
using namespace std;
# include "hb_io.H"
# include "multi_hb.h"

void bicgstab(double [],int [],int [],double [],double [],int,int,double [],double []);
double vecmul(double [],double [],int);
void matmul(double [],int [],int [],double [],int,int,double []);

int main ( void )
{
  //Variables
  int *colptr = NULL;
  double *exact = NULL;
  double *guess = NULL;
//  int i;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
//  int j;
  char *key = NULL;
//  int khi;
//  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero; //No. of non-zeros
  int nrhs;   //rhs= right hand side [B]
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  int *rhsind = NULL;
  int *rhsptr = NULL;
  char *rhstyp = NULL;
  double *rhsval = NULL;
  double *rhsvec = NULL;
  int *rowind = NULL; //Row index
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  double *values = NULL;

  input.open ( "input.txt" ); // inputed file
  ofstream output;
  output.open("output.txt");  // output file


  if ( !input )
  {
    cout << "\n";
    cout << "  Error opening the file!\n";
    cout << " (input.txt not found)\n \n";
    getch();
    return 0;
  }

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd,
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = new int[ncol+1];

  if ( mxtype[2] == 'A' )
  {
    rowind = new int[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = new int[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    getch();
    exit ( 1 );
  }

// Matrix general information:
cout<<"Matrix General Informations: \n";
cout<<"Title: "<<title<<'\n';
cout<<"Key name="<<key<<'\n';
cout<<"Matrix type="<<mxtype<<'\n';
cout<<"Matrix is "<<nrow<<" * "<<nrow<<'\n';
cout<<"Number of non-zero elements= "<<nnzero<<"\n \n";


  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl,
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  if ( mxtype[2] == 'A' )
  {
    values = new double[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    values =  new double[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    getch();
    exit ( 1 );
  }

cout<<"Reading the Structure... Done! \n";

//   [A] values.

  hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );
  cout<<"Reading the Values... Done! \n";


//  Right hand sides.

if (0==rhscrd){
      cout<<"\n Error! \n Not a Standard Right Hand Side. \n";
      getch();
      return 0;
      }

  if ( 0 < rhscrd )
  {

    if ( rhstyp[0] == 'F' )
    {
      rhsval = new double[nrow*nrhs];
    }
    else if ( rhstyp[0] == 'M' )
    {
      // Not Standard!
      cout<<"\n Fatal Error! \n Not a Standard Right Hand Side. \n";
      getch();
      return 0;

      if ( mxtype[2] == 'A' )
      {
        rhsptr = new int[nrhs+1];
        rhsind = new int [nrhsix];
        rhsvec = new double[nrhsix];
      }
      else if ( mxtype[2] == 'E' )
      {
        rhsval = new double [nnzero*nrhs];
      }

      // Not standard!
       else
      {
      cout<<"\n Fatal Error! \n Not a Standard Right Hand Side. \n";
      getch();
      return 0;
      }
    }

    hb_rhs_read ( input, nrow, nnzero, nrhs, nrhsix,
      rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval,
      rhsind, rhsptr, rhsvec );

    cout << "Reading the right hand side... Done! \n \n";


//  Read the starting guesses if present.

    if ( rhstyp[1] == 'G' )
    {
      guess = new double[nrow*nrhs];

      hb_guess_read ( input, nrow, nrhs, rhscrd, rhsfmt, rhstyp, guess );

    }

//  Read the exact solutions if present.

    if ( rhstyp[2] == 'X' )
    {
      exact = new double[nrow*nrhs];

      hb_exact_read ( input, nrow, nrhs, rhscrd, rhsfmt, rhstyp, exact );

    }
  }

//************************ Reading is completed ***********************

//***************** Generating Column Position Vector *****************

int *colind;  // Column index
colind=new int [nnzero];

int dd=0;
int kk=0;
for (int ii=0;ii<=nrow-1;ii++){
    dd=colptr[ii+1]-colptr[ii];
    for (int jj=1;jj<=dd;jj++){
        colind[kk]=ii+1;
        kk=kk+1;
    }
}


//  Now; values[] , rowind[] & colind[] determined A matrix, and rhsval[] is b


//*********************** Calling BiCGSTAB Solver ******************************

cout<<"Matrix reading... Completed! \n \n"<<"BiCGSTAB Solver: \n \n";

double *tmpv;        //tmpv is a temp vector to carry [A].v (matmul answer)
tmpv=new double [nrow];
double *fguess;
fguess=new double [nrow];
double *bix;                   //x
bix=new double [nrow];
for (int iii=0;iii<nrow;iii++){     //First guess =0
    bix[iii]=fguess[iii]=0;
}

// Solving:
bicgstab(values,rowind,colind,rhsval,bix,nrow,nnzero,tmpv,fguess);


// Saving {x} and residuals
for (int iii=0;iii<nrow;iii++){
    output<<"x["<<iii+1<<"]= "<<bix[iii]<<"\n";
    }
cout<<"\n Residuals Writed to [residual.txt] \n";
cout<<"\n Answer Vector {x} Writed to [output.txt] \n";
  getch();     // pause

  return 0;
}


//*********************** The BiCGSTAB Solver ****************************

void bicgstab(double values[],int rowind[],int colind[],double rhsval[],double bix[],
              int nrow,int nnzero,double tmpv[],double fguess[]){

  ofstream residual;
  residual.open("residual.txt");         //Residual output
  residual<<"Iter  Residual\n";

// 0)    define variables :
double *bir0;                  //r0
bir0=new double [nrow];
double *bir;                   //r
bir=new double [nrow];
double *birr;                  //rr
birr=new double [nrow];
double *bip;                   //p
bip=new double [nrow];
double *bis;                   //s
bis=new double [nrow];
//double *bix;                   //x
for (int iii=0;iii<nrow;iii++){     // =0
    bir0[iii]=bir[iii]=birr[iii]=bip[iii]=bis[iii]=bix[iii]=0;
}
double bialpha=0;              //alpha
double biw=0;                  //w
double bibeta=0;               //beta

// 1)
for (int iii=0;iii<nrow;iii++){    // First Guess
    bix[iii]=fguess[iii];
}

matmul(values,rowind,colind,bix,nrow,nnzero,tmpv);
for (int iii=0;iii<nrow;iii++){
    bir0[iii]=rhsval[iii]-tmpv[iii];
}

// 2)

for (int iii=0;iii<nrow;iii++){
    bip[iii]=bir0[iii];
    bir[iii]=bir0[iii];
}

// 3)
int noiter=0;
double eps=0;
cout<<"Enter No. of Max Iterations (Integer)=";
cin>>noiter;
cout<<"\n"<<"Enter Epsilon (Real)="; // Minimum acceptable residual
cin>>eps;

for (int biter=0;biter<noiter;biter++){

    // 4)
    matmul(values,rowind,colind,bip,nrow,nnzero,tmpv);
    bialpha=(vecmul(bir,bir0,nrow))/(vecmul(tmpv,bir0,nrow));

    // 5)
    for (int iii=0;iii<nrow;iii++){
    bis[iii]=bir[iii]-bialpha*tmpv[iii];
    }

    // 6)
    matmul(values,rowind,colind,bis,nrow,nnzero,tmpv);
    biw=(vecmul(tmpv,bis,nrow))/(vecmul(tmpv,tmpv,nrow));

    // 7)
    for (int iii=0;iii<nrow;iii++){
    bix[iii]=bix[iii]+bialpha*bip[iii]+biw*bis[iii];
    }

    // 8)
    for (int iii=0;iii<nrow;iii++){
    birr[iii]=bis[iii]-biw*tmpv[iii];
    }

    // 9)
    bibeta=((vecmul(birr,bir0,nrow))/(vecmul(bir,bir0,nrow)))*(bialpha/biw);

    // 10)
    matmul(values,rowind,colind,bip,nrow,nnzero,tmpv);
    for (int iii=0;iii<nrow;iii++){
    bip[iii]=birr[iii]+bibeta*(bip[iii]-biw*tmpv[iii]);
    }

    // 11)
    for (int iii=0;iii<nrow;iii++){
    bir[iii]=birr[iii];
    }

    // residuals:
    double bires=sqrt((vecmul(bir,bir,nrow))/nrow);
    cout<<"iter="<<biter+1<<"  Res.="<<bires<<'\n';
    //if (biter%3==0){
    residual<<biter+1<<"  "<<bires<<"\n";
    //}
    if (bires<=eps){
    cout<<"\n Solution Converged!"<<'\n';
    break;
    }
}

return;
}

//*********** Multiply A matrix to V vector (efficiently!): *************

void matmul(double values[],int rowind[],int colind[],double vec[],int nrow,int nnzero,double tmpv[]){

 for (int i=0;i<nrow;i++){
    tmpv[i]=0;             // Initialize the output tmpv[]
 }

 for (int j=0;j<nnzero;j++){  //:)
    tmpv[rowind[j]-1]=tmpv[rowind[j]-1]+values[j]*vec[colind[j]-1];
 }

 return;
}


//*********************** Dot product: ***************************
double vecmul(double a[],double b[],int nrow){
 double mul=0;
      for (int i=0;i<nrow;i++){
          mul=mul+a[i]*b[i];
      }
 return mul;
}
