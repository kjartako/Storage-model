#include <RcppEigen.h>
#include "Eigen/Dense"
#include <math.h>
#include <cmath>
using namespace Eigen;
using namespace std;
using namespace Rcpp;

const double unitol=pow(std::numeric_limits<double>::epsilon(),0.25);
const double c_s0 = 0.01;
const double c_df = 10;
const double c_alpha = 2;
const double c_beta = 20;

using mat = Eigen::Matrix<double,Dynamic,Dynamic>;
using vec = Eigen::Matrix<double,Dynamic,1>;

/////////////////////////////
double SIGN(const double &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
/////////////////////////////
// [[Rcpp::export]]
Eigen::VectorXd systematic_resampling(Eigen::VectorXd x, int n, Eigen::VectorXd w)
{
  //Cumulative sum
  for(int f=1; f<n; f++)
  {
    w(f) += w(f-1);
  }
  
  vec U(n);
  U(0) = R::runif(0,1.0/n);
  
  for(int j=1;j<n;j++)
  {
    U(j)=U(0)+((double)j/(double)n);
  }
  
  vec Ncount(n);
  
  int h = 0;
  
  for(int i=0;i<n;i++)
  {
    int j = h;
    while (U(h) < w(i) && h < n)
    {
      h = h + 1; 
    }
    Ncount(i) = h - j;
  }
  
  vec res(n);
  
  //rep(x,Ncount)
  int ind=0;
  for (int k=0; k < n; k++)
  {
    int p = Ncount(k);
    res.segment(ind,p).fill(x(k));
    ind += p;
  }

  return(res);
}
/////////////////////////////
// [[Rcpp::export]]
Eigen::VectorXd sample_replace(int n, int size, Eigen::VectorXd p)
{
  // Normalize weights
  
  double sum = p.sum();
  p = p/sum;
  
  // Sample using normalized weights
  
  vec a(n);
  vec ans(size);
  
  int i, j, k;
  std::vector<double> q(n);
  double rU;
  
  std::vector<int> HL(n);
  std::vector<int>::iterator H, L;
  
  H = HL.begin() - 1; L = HL.begin() + n;
  for (i = 0; i < n; i++) {
    q[i] = p(i) * n;
    if (q[i] < 1.0) {
      *++H = i;
    } else {
      *--L = i;
    }
  }
  
  if (H >= HL.begin() && L < HL.begin() + n) {
    for (k = 0; k < n - 1; k++) {
      i = HL[k];
      j = *L;
      a(i) = j;
      q[j] += q[i] - 1;
      
      L += (q[j] < 1.0);
      
      if (L >= HL.begin() + n) {
        break;
      }
    }
  }
  
  for (i = 0; i < n; i++) {
    q[i] += i;
  }
  
  for (i = 0; i < size; i++) {
    rU = unif_rand() * n;
    k = static_cast<int>(rU);
    ans(i) = (rU < q[k]) ? k : a(k);
  }
  
  return(ans);
  
}
/////////////////////////////
template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol, double xleft, double xright, vec grid, vec vals, double x, double apar, double bpar, double Cpar, double delta, double beta, vec intgrid, vec intwts, int nint)
{
  const int ITMAX=100;
  const double EPS=std::numeric_limits<double>::epsilon(); 
  double a=x1,b=x2,c=x2,d,e,fa=func(a,xleft,xright,grid,vals, x, apar, bpar, Cpar, delta, beta, intgrid, intwts, nint),fb=func(b,xleft,xright,grid,vals,x, apar, bpar, Cpar, delta, beta, intgrid, intwts, nint),fc,p,q,r,s,tol1,xm;
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    throw("Root must be bracketed in zbrent");
  fc=fb;
  for (int iter=0;iter<ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (abs(fc) < abs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*abs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (abs(xm) <= tol1 || fb == 0.0) return b;
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=abs(p);
      double min1=3.0*xm*q-abs(tol1*q);
      double min2=abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (abs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=func(b,xleft,xright,grid,vals, x, apar, bpar, Cpar, delta, beta, intgrid, intwts, nint);
  }
  throw("Maximum number of iterations exceeded in zbrent");
}

/////////////////////////////
template <class T>
double zbrent2(T &func, const double x1, const double x2, const double tol, double Dval, double xleft, double xright, vec grid, vec vals, double apar, double bpar, double Cpar)
{
  const int ITMAX=100;
  const double EPS=std::numeric_limits<double>::epsilon(); 
  double a=x1,b=x2,c=x2,d,e,fa=func(a,Dval,xleft,xright,grid,vals, apar, bpar, Cpar),fb=func(b,Dval,xleft,xright,grid,vals, apar, bpar, Cpar),fc,p,q,r,s,tol1,xm;
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    throw("Root must be bracketed in zbrent");
  fc=fb;
  for (int iter=0;iter<ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (abs(fc) < abs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*abs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (abs(xm) <= tol1 || fb == 0.0) return b;
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=abs(p);
      double min1=3.0*xm*q-abs(tol1*q);
      double min2=abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (abs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=func(b,Dval,xleft,xright,grid,vals, apar, bpar, Cpar);
  }
  throw("Maximum number of iterations exceeded in zbrent");
}

/////////////////////////////

vec seq(double from, double to, int n)
{
  vec res(n);
  res(0)=from;
  res(n-1)=to;
  
  double diff = (to-from)/(n-1);
  
  for(int i=1;i<(n-1);i++)
  {
    res(i)= res(i-1)+diff;
  }
  
  return(res);
}

/////////////////////////////

double trapes(vec w, vec f, int li)
{
  double tmp = 2*w.dot(f);
  return(tmp-(w(0)*f(0))-(w(li)*f(li)));
}

/////////////////////////////

double LinInter(double x, vec xgrid, vec ygrid)
{
  double diff = xgrid(1)-xgrid(0);
  int j = floor((x-xgrid(0))/diff);
  
  return(ygrid[j]+((x-xgrid[j])/(xgrid[j+1]-xgrid[j]))*(ygrid[j+1]-ygrid[j]));
}

/////////////////////////////

double P(double x, double apar, double bpar)
{
  return(exp(apar-bpar*x));
}
/////////////////////////////

double P_dx(double x, double apar, double bpar)
{
  return(-bpar*exp(apar-bpar*x));
}
/////////////////////////////

double logP(double x, double apar, double bpar)
{
  return(apar-bpar*x);
}

/////////////////////////////

double D(double p, double apar, double bpar)
{
  return((apar-log(p))/bpar);
}

/////////////////////////////
// [[Rcpp::export]]
double evalSigma(double x, double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double C)
{
  double ret;

  if(x < xleft)
  {
    ret = 0.0;
  }
  else if(x > xright)
  {
    ret = C;
  }
  else
  {
    ret = LinInter(x,grid,vals);
  }
  return(ret);
}
/////////////////////////////
// [[Rcpp::export]]
double evalf(double x, double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double apar, double bpar,double Cpar)
{
  double sig;
  
  if(x < xleft)
  {
    sig = 0.0;
  }
  else if(x > xright)
  {
    sig = Cpar;
  }
  else
  {
    sig = LinInter(x,grid,vals);
  }
  
  double val = x - sig;
  
  return(P(val,apar,bpar));
}

/////////////////////////////
// [[Rcpp::export]]
double evalf_dx(double x, double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double apar, double bpar,double Cpar)
{
  double res;
  
  if(x < xleft)
  {
    res = P_dx(x,apar,bpar);
  }
  else if(x > xright)
  {
    res = P_dx(x-Cpar,apar,bpar);
  }
  else
  {
	  double diff = grid(1)-grid(0);
	  int j = floor((x-grid(0))/diff);
	  double slope = (vals[j+1]-vals[j])/(grid[j+1]-grid[j]);
	  double sigma = vals[j]+(x-grid[j])*slope;
	  
    res = P_dx(x-sigma,apar,bpar)*(1-slope);
  }
  
  return(res);
}
/////////////////////////////
double xEq(double x, double Dval,double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double apar, double bpar,double Cpar)
{
  return(x-evalSigma(x,xleft, xright, grid, vals,Cpar)-Dval);
}
/////////////////////////////
// [[Rcpp::export]]
double finv(double fval,double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double apar, double bpar,double Cpar)
{
  double Dval = D(fval,apar,bpar);
  
  double res = zbrent2(xEq, -200, 200, unitol, Dval, xleft, xright, grid, vals, apar, bpar, Cpar);
  
  return(res);
  
}
/////////////////////////////
// [[Rcpp::export]]
double evalf_log(double x, double sigmaeval, double apar, double bpar)
{
  double arg = x - sigmaeval;
  
  return(logP(arg,apar,bpar));
}

/////////////////////////////
double sigmaEq(double sigma, double xleft, double xright, vec grid, vec vals, double x, double apar, double bpar, double Cpar, double delta, double beta, vec intgrid, vec intwts, int nint)
{
  vec tmp(nint);

  for (int i=0;i<nint;i++)
  {
    tmp(i) = evalf((1.0-delta)*sigma+intgrid(i),xleft,xright,grid,vals,apar,bpar,Cpar);
  }
  
  double eq = x - D(beta*trapes(intwts,tmp,nint-1),apar,bpar) - sigma;
    
  return(eq);  
}
/////////////////////////////
double qt(double mean, double v_x)
{
  return(R::rnorm(mean,v_x));
}
/////////////////////////////
double fx_c(double x, double xprev, double xleft, double xright, vec grid, vec vals, double C, double delta, double v_x)
{
  double mean = (1.0-delta)*evalSigma(xprev,xleft,xright,grid,vals,C);
  double res = R::dnorm(x,mean,v_x,0);
  return(res);
}
/////////////////////////////

// w_f = g - px = p(p_t,x_t|p_{t-1},x_{t-1}) - p(x_t|x_{t-1})
double w_f(double p, double p_prev,double x, double xprev, double a, double b, double v, double v_x, double sig, double sigprev,double mean)
{
  //g-px 
  double res = (-0.5*log(2*PI)-log(v))+(-1.0/(2.0*pow(v,2)))*pow((p-p_prev-evalf_log(x,sig,a,b)+evalf_log(xprev,sigprev,a,b)),2);//-(0.5*pow(x-mean,2))+(0.5*pow(x-mean,2));
  
  return(res);
}

/////////////////////////////

double w_rho(double rho, double p, double p_prev,double x, double xprev, double a,double b,double v,double v_x,double sig, double sigprev, double mean)
{
	double tmpval = pow(v,2)*(1-pow(rho,2));
	double A = x-mean;
	double B = p-p_prev-evalf_log(x,sig,a,b)+evalf_log(xprev,sigprev,a,b);
	
	double g = -0.5*log(2*PI*tmpval) - 0.5*((pow(A*v,2)+pow(B,2)-2*rho*v*A*B)/tmpval);
	
	double px = R::dnorm(x,mean,v_x,1);
	
	return(exp(g-px));
}

/////////////////////////////
vec vecexp(vec w, int npar)
{
  vec res(npar);
  
  for(int i=0;i<npar;i++)
  {
    res(i)=exp(w(i));
  }
  
  return(res);
}

/////////////////////////////
double ESScalc(vec wnorm, int npar)
{
  double Wsumsq = 0;
  
  for(int h=0;h<npar;h++)
  {
    Wsumsq = Wsumsq + pow(wnorm(h),2);
  }
  double ESSval = 1.0/Wsumsq;
  
  return(ESSval);
}

/////////////////////////////
// [[Rcpp::export]]
Rcpp::List Storagefunc(double xleft,double xright,Eigen::VectorXd grid, Eigen::VectorXd vals, double apar, double bpar, double Cpar, double delta, double beta, Eigen::VectorXd intgrid, Eigen::VectorXd intwts, int nint, int ngrid, int niter)
{
  for (int i=0;i<niter;i++)
  {
    // locate left endpoint
    vec tmp(nint);

    for (int j=0;j<nint;j++)
    {
      tmp(j) = evalf(intgrid(j),xleft,xright,grid,vals,apar,bpar,Cpar);
    }

    double xleft_new = D(beta*trapes(intwts,tmp,nint-1),apar,bpar);

    // locate right endpoint
    for (int k=0;k<nint;k++)
    {
      tmp(k) = evalf((1.0-delta)*Cpar+intgrid(k),xleft,xright,grid,vals,apar,bpar,Cpar);
    }
    
    double xright_new = D(beta*trapes(intwts,tmp,nint-1),apar,bpar)+Cpar;

    vec grid_new = seq(xleft_new,xright_new,ngrid);

    vec vals_new(ngrid);
    vals_new(0) = 0.0;
    vals_new(ngrid-1) = Cpar;

    for (int m=1;m<(ngrid-1);m++)
    {
      vals_new(m) = zbrent(sigmaEq, 0, Cpar, unitol, xleft, xright, grid, vals, grid_new(m), apar, bpar, Cpar, delta, beta, intgrid, intwts, nint);
    }

    xleft=xleft_new;
    xright=xright_new;
    grid=grid_new;
    vals=vals_new;
  }
  
  Rcpp::List rep;
  rep["xleft"]=xleft;
  rep["xright"]=xright;
  rep["grid"]=grid;
  rep["vals"]=vals;
  
  return(rep);
}

/////////////////////////////

// [[Rcpp::export]]
Rcpp::List SMC_sample(Eigen::VectorXd theta, Eigen::VectorXd p, int npar, int N, double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double a, double C, double r, double v_x, Eigen::VectorXd intgrid, Eigen::VectorXd intwts, int nint, int ngrid, int niter, double resampthresh)
{
  double loglike = 0;
  vec logw(npar);

  double v = exp(theta(0));
  double delta = 0.5+0.5*tanh(theta(1));
  double b = exp(theta(2));
  double beta = (1.0-delta)/(1.0+r);

  Rcpp:List rep = Storagefunc(xleft,xright,grid,vals,a,b,C,delta,beta,intgrid,intwts,nint,ngrid,niter);
  
  double xleft = rep["xleft"];
  double xright = rep["xright"];
  vec grid = rep["grid"];
  vec vals = rep["vals"];

  mat wmat(npar,N-1);
  
  // Time t=1
  mat Xp(npar,N);

  for (int i=0;i<npar;i++)
  {
    Xp(i,0)= R::runif(-2.0,2.0+C);
  }

  logw.setZero();
  
  // vec w = vecexp(logw,npar);
  // 
  // loglike = loglike + log(w.mean());

  // Time t = 2,3,..,N

  for(int t=1;t<N;t++)
  {
  	for(int k=0;k<npar;k++)
    {
        // Sample new particles
        double sigprev = evalSigma(Xp(k,t-1),xleft,xright,grid,vals,C);
        double mean = (1.0-delta)*sigprev;
        Xp(k,t)=qt(mean,v_x);
        double sig = evalSigma(Xp(k,t),xleft,xright,grid,vals,C);
        // Compute weights, and add to total logweight since last resampling
        logw(k)=logw(k)+w_f(p(t),p(t-1),Xp(k,t),Xp(k,t-1),a,b,v,v_x,sig,sigprev,mean);
    }
	
	  wmat.col(t-1)=logw;
    
    // calculate normalized weights
    double maxtmp = logw.maxCoeff();
    vec wprop(npar);
    
    for (int u=0;u<npar;u++)
    {
      wprop(u) = exp(logw(u)-maxtmp);
    }
    
    double wprop_sum = wprop.sum();
    vec wnorm = wprop/wprop_sum;
    
    double ESSval = ESScalc(wnorm, npar);
    
    if(ESSval < (resampthresh*npar))
    {
      // Resample particles
      
      Xp.col(t)=systematic_resampling(Xp.col(t),npar,wnorm);
  	  	
	  	// Update loglike, using log-sum-exp

	  	loglike = loglike + (maxtmp +log(wprop_sum/npar));	
  	  
  	  logw.setZero();
    }
  }
  
  // Update loglike, using log-sum-exp
  double logwmax = logw.maxCoeff();
  double expsum = 0;
  
  for (int h=0;h<npar;h++)
  {
    expsum = expsum + exp(logw(h) - logwmax);
  }
  
  loglike = loglike + (logwmax +log(expsum/npar));
  
  double v_prior = -c_df*log(v)+0.5*c_df*log(0.5*c_s0*c_df)-log(tgamma(0.5*c_df))+log(2.0)-(0.5*c_df*c_s0)/pow(v,2);
  double delta_prior = -log(tgamma(c_alpha))-log(tgamma(c_beta))+log(tgamma(c_alpha+c_beta))+log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*pow(2.0*delta-1.0,2)));
  double b_prior = -0.5*log(2.0*PI)-0.5*pow(theta(2),2);

  Rcpp::List res;
  res["loglike_noprior"]=loglike;
  res["loglike"]=loglike+v_prior+delta_prior+b_prior;
  //res["Xp"]=Xp;
  //res["logw"]=wmat;
  //res["rep"]=rep;

  return(res);
  
}

/////////////////////////////

// [[Rcpp::export]]
Rcpp::List SMC_sample_smooth(Eigen::VectorXd theta, Eigen::VectorXd p, int npar, int nrep, int N, double xleft, double xright, Eigen::VectorXd grid, Eigen::VectorXd vals, double a, double C, double r, double v_x, Eigen::VectorXd intgrid, Eigen::VectorXd intwts, int nint, int ngrid, int niter, double resampthresh)
{
  double loglike = 0;
  vec logw(npar);
  
  double v = exp(theta(0));
  double delta = delta = 0.5+0.5*tanh(theta(1));
  double b = exp(theta(2));
  double beta = (1.0-delta)/(1.0+r);
  
  Rcpp:List rep = Storagefunc(xleft,xright,grid,vals,a,b,C,delta,beta,intgrid,intwts,nint,ngrid,niter);
  
  double xleft = rep["xleft"];
  double xright = rep["xright"];
  vec grid = rep["grid"];
  vec vals = rep["vals"];
  
  mat wmat(npar,N);
  
  // Time t=1
  mat Xp(npar,N);
  
  for (int i=0;i<npar;i++)
  {
    Xp(i,0)= R::runif(-2.0,2.0+C);
  }
  
  logw.setZero();
  
  wmat.col(0)=vecexp(logw,npar);
  
  // vec w = vecexp(logw,npar);
  // 
  // loglike = loglike + log(w.mean());
  
  // Time t = 2,3,..,N
  
  for(int t=1;t<N;t++)
  {
    for(int k=0;k<npar;k++)
    {
      // Sample new particles
      double sigprev = evalSigma(Xp(k,t-1),xleft,xright,grid,vals,C);
      double mean = (1.0-delta)*sigprev;
      Xp(k,t)=qt(mean,v_x);
      double sig = evalSigma(Xp(k,t),xleft,xright,grid,vals,C);
      // Compute weights, and add to total logweight since last resampling
      logw(k)=logw(k)+w_f(p(t),p(t-1),Xp(k,t),Xp(k,t-1),a,b,v,v_x,sig,sigprev,mean);
    }
    
    // calculate normalized weights
    double maxtmp = logw.maxCoeff();
    vec wprop(npar);
    
    for (int u=0;u<npar;u++)
    {
      wprop(u) = exp(logw(u)-maxtmp);
    }
    
    double wprop_sum = wprop.sum();
    vec wnorm = wprop/wprop_sum;
    
    wmat.col(t)=wnorm;
    
    double ESSval = ESScalc(wnorm, npar);
    
    if(ESSval < (resampthresh*npar))
    {
      // Resample particles
      
      Xp.col(t)=systematic_resampling(Xp.col(t),npar,wnorm);
      
      // Update loglike, using log-sum-exp
      
      loglike = loglike + (maxtmp +log(wprop_sum/npar));	
      
      logw.setZero();
    }
  }
  
  // Update loglike, using log-sum-exp
  double logwmax = logw.maxCoeff();
  double expsum = 0;
  
  for (int h=0;h<npar;h++)
  {
    expsum = expsum + exp(logw(h) - logwmax);
  }
  
  loglike = loglike + (logwmax +log(expsum/npar));
  
  double v_prior = ((-0.5*c_df*c_s0)/pow(v,2))-c_df*log(v);
  double delta_prior = log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*pow(2.0*delta-1.0,2)));
  //double b_prior = -0.5*pow(b,2)+log(b);
  double b_prior = -0.5*pow(theta(2),2);
  
  
  //////////////////////////////////

  // smoothing  
  mat Xsmooth(nrep,N);
  vec wsmooth(npar);
  vec tmpT(1);
  vec tmpj(1);

  for(int h=0;h<nrep;h++)
  {
    tmpT = sample_replace(npar,1, wmat.col(N-1));

    Xsmooth(h,N-1)=Xp(tmpT(0),N-1);

    for(int j=(N-2);j>=0;j--)
    {
      for(int i=0;i<npar;i++)
      {
        wsmooth(i) = wmat(i,j)*fx_c(Xsmooth(h,j+1), Xp(i,j), xleft, xright,grid,vals, C, delta,v_x);
      }

      tmpj=sample_replace(npar,1, wsmooth);
      Xsmooth(h,j)=Xp(tmpj(0),j);
    }

  }
  
  Rcpp::List res;
  res["loglike_noprior"]=loglike;
  res["loglike"]=loglike+v_prior+delta_prior+b_prior;
  res["Xp"]=Xp;
  res["logw"]=wmat;
  res["Xsmooth"]=Xsmooth;
  res["rep"]=rep;
  
  return(res);
  
}