#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
typedef std::vector<double> stdvec;
// [[Rcpp::export]]
SEXP  slidingWindow(NumericVector slidestarts, DataFrame& df,double window, double binx, double biny,Rcpp::Function func, SEXP obj,arma::sp_mat m_tile,StringVector bc, StringVector gene, int tile,int groupid) {
  int n = slidestarts.length();
  for (int j=0; j<n;j++)
  {
    //std::cout<<j<<std::endl;
    for (int t=0;t<n;t++) {
    //std::cout<<t<<std::endl;

      //double miny = min(df['y_miseq']);
      obj=func(df,binx,biny,window,j,t,obj,m_tile,bc,gene,tile,groupid);
    }
  }
  return obj;

 }

