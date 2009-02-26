// -------------------------------------------
// fitting functions
// -------------------------------------------

// crystal ball fit
double crystalball(double *x, double *par) {
  // par[0]:  mean
  // par[1]:  sigma
  // par[2]:  alpha, crossover point
  // par[3]:  n, length of tail
  // par[4]:  N, normalization
                                
  double cb = 0.0;
  double exponent = 0.0;
  
  if (x[0] > par[0] - par[2]*par[1]) {
    exponent = (x[0] - par[0])/par[1];
    cb = exp(-exponent*exponent/2.);
  } else {
    double nenner  = pow(par[3]/par[2], par[3])*exp(-par[2]*par[2]/2.);
    double zaehler = (par[0] - x[0])/par[1] + par[3]/par[2] - par[2];
    zaehler = pow(zaehler, par[3]);
    cb = nenner/zaehler;
  }
  
  if (par[4] > 0.) {
    cb *= par[4];
  }
  return cb;
}

// crystal ball fit
double cruijffFunction(double *x, double *par) {

  double dx    = 0.0; 
  double sigma = 0.0; 
  double alpha = 0.0; 
  double f     = 0.0; 
  double func  = 0.0;

  dx = (x[0] - par[0]);  
  sigma = dx<0 ? par[1]: par[1] ;  
  alpha = dx<0 ? par[2]: par[3] ;  
  f = 2*sigma*sigma + alpha*dx*dx ;  
  func = exp(-dx*dx/f); 

  if (par[4] > 0.) {
    func *= par[4];
  }

  return func;
}


