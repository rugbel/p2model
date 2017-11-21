#define TMB_LIB_INIT R_init_mypkg
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* data section */
  DATA_MATRIX(XnS);
  DATA_MATRIX(XnR);
  DATA_ARRAY(XvD);
  DATA_ARRAY(XvC);
  DATA_MATRIX(y); //network
  DATA_INTEGER(penflag); // 0 for no penalty 1 for penalized log-likelihood 2 for informative penalization
  DATA_VECTOR(penSigma); // constant matrix (three-dim vector) to use in informative penalization
  

  
  /* Parameter section */
  PARAMETER_VECTOR(gamma1);
  PARAMETER_VECTOR(gamma2);
  PARAMETER_VECTOR(delta1);
  PARAMETER_VECTOR(delta2); 
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(b);

  int kd = XvD.cols();
  int kc = XvC.cols();
  int g = y.cols();

  
  Type nll=0.0;     // Negative log likelihood function


  vector<Type> S(3);
  vector<Type> S1(3);
  S(0) = alpha(0) * alpha(0);   
  S(1) = alpha(0) * alpha(1);
  S(2) = alpha(1) * alpha(1) + alpha(2) * alpha(2); 
  Type DS =  S(0) * S(2) - S(1) * S(1);
  S1(0) = S(2) / DS;
  S1(1) = - S(1) / DS;
  S1(2) = S(0) / DS;
  Type quadu = ((a*a).sum()) * S1(0) + ((a*b).sum()) * S1(1) * 2.0 + ((b*b).sum()) * S1(2);
  nll -= -0.5 * quadu - g * 0.5 * log(DS);

  
															   
  ADREPORT(S);
  
  // nll from y
  vector<Type> XSg1 = XnS * gamma1; //alpha_i in the model
  vector<Type> XRg2 = XnR * gamma2; //beta_i 
 
 
  for(int i=0;i<(g-1);i++)
    for(int j=i+1;j<g;j++)
    {
     Type y1 = y(i,j);
     Type y2 = y(j,i);
     Type alphai = XSg1(i) + a(i);
     Type alphaj = XSg1(j) + a(j);
     Type betai = XRg2(i) + b(i);
     Type betaj = XRg2(j) + b(j);
     Type muij = 0.0;
     Type muji = 0.0;
     for(int k=0;k<kd;k++)
       {
        muij += XvD(i,j,k) * delta1(k);
	muji += XvD(j,i,k) * delta1(k);
       } 
     Type rhoij = 0.0;
     for(int k=0;k<kc;k++)
	  rhoij += XvC(i,j,k) * delta2(k);
     Type xi1 =  muij + alphai + betaj;
     Type xi2 =  muji + alphaj + betai;
     Type xi3 =  muij + muji + alphai + betaj + alphaj + betai + rhoij;
     Type den = 1.0 + exp(xi1) + exp(xi2) + exp(xi3);
     nll -=  y1 * xi1  + y2 * xi2 + rhoij * y1 * y2 - log(den);
  }

  Type pen = 0.0;
  Type rho = S(1) / sqrt(S(0) * S(2));
  
   if(penflag==1)
     pen = 0.5 * log(DS);

   if(penflag==2)
    {
      Type rhoinfo = penSigma(1) / sqrt(penSigma(0) * penSigma(2));
      pen = 0.5 * log(DS) + Type(2) * log(S(0)) + Type(2) * log(S(2))  - (Type(2) / sqrt(penSigma(0))) * sqrt(S(0)) -  (Type(2) / sqrt(penSigma(2))) * sqrt(S(2))  - pow(rho-rhoinfo, Type(2)) / 0.125; 
    }

   nll -= pen;

   ADREPORT(rho);
   REPORT(pen);
   
  return nll;
}
