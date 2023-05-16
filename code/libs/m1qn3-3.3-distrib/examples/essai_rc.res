 -------------------------------------------------------------------------------
 M1QN3 Copyright (C) 2008, J. Ch. Gilbert, Cl. Lemarechal.
 -------------------------------------------------------------------------------
 This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you
 are welcome to redistribute it under certain conditions. See the file COPYING 
 in the root directory of the M1QN3 distribution for details.
 -------------------------------------------------------------------------------

 M1QN3 (Version 3.3, October 2009): entry point
     dimension of the problem (n):             2
     absolute precision on x (dxmin):          1.00D-10
     expected decrease for f (df1):            6.28D+03
     relative precision on g (epsg):           1.00D-05 (dfn-norm)
     maximal number of iterations (niter):   200
     maximal number of simulations (nsim):   200
     printing level (impres):                  5
     reverse communication

 m1qn3: Diagonal Initial Scaling mode

     allocated memory (ndz) :    20000
     used memory :               19998
     number of updates :          3998
     (y,s) pairs are stored in core memory

 m1qn3: cold start

     f             =  6.27500000D+03
     dfn-norm of g =  5.00001000D+03

 m1qn3a: descent direction -g: precon =  0.502D-03

 -------------------------------------------------------------------------------

 m1qn3: iter 1, simul 1, f= 6.27500000D+03, h'(0)=-1.25500D+04

 m1qn3: line search

     mlis3       fpn=-1.255D+04 d2= 6.30D+00  tmin= 3.98D-11 tmax= 1.00D+20
     mlis3     1.000D+00 -5.866D+03 -1.550D+03

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  5.727D-04
     updated diagonal: average value =  5.727D-04, sqrt(variance) =  2.291D-09

 m1qn3: stopping criterion on g:  1.23523D-01

 m1qn3: descent direction d: angle(-g,d) =   0.1 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 2, simul 2, f= 4.09368426D+02, h'(0)=-2.18485D+02

 m1qn3: line search

     mlis3       fpn=-2.185D+02 d2= 1.25D-01  tmin= 2.83D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -1.762D+02 -1.380D+02

 m1qn3: convergence rate, s(k)/s(k-1) =  1.40939D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  2.714D+00
     updated diagonal: average value =  1.554D-03, sqrt(variance) =  5.195D-07

 m1qn3: stopping criterion on g:  7.80236D-02

 m1qn3: descent direction d: angle(-g,d) =   1.0 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 3, simul 3, f= 2.33171494D+02, h'(0)=-2.36787D+02

 m1qn3: line search

     mlis3       fpn=-2.368D+02 d2= 3.69D-01  tmin= 1.65D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -1.538D+02 -8.712D+01

 m1qn3: convergence rate, s(k)/s(k-1) =  1.71604D+00

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.582D+00
     updated diagonal: average value =  2.458D-03, sqrt(variance) =  5.504D-06

 m1qn3: stopping criterion on g:  2.87118D-02

 m1qn3: descent direction d: angle(-g,d) =   2.5 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 4, simul 4, f= 7.93960556D+01, h'(0)=-5.09545D+01

 m1qn3: line search

     mlis3       fpn=-5.095D+01 d2= 1.26D-01  tmin= 2.83D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -3.599D+01 -2.341D+01

 m1qn3: convergence rate, s(k)/s(k-1) =  5.85225D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.845D+00
     updated diagonal: average value =  4.536D-03, sqrt(variance) =  6.687D-05

 m1qn3: stopping criterion on g:  1.31855D-02

 m1qn3: descent direction d: angle(-g,d) =   6.2 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 5, simul 5, f= 4.34069529D+01, h'(0)=-2.03263D+01

 m1qn3: line search

     mlis3       fpn=-2.033D+01 d2= 9.62D-02  tmin= 3.33D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -1.404D+01 -8.854D+00

 m1qn3: convergence rate, s(k)/s(k-1) =  8.72873D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.758D+00
     updated diagonal: average value =  8.012D-03, sqrt(variance) =  6.226D-04

 m1qn3: stopping criterion on g:  5.73203D-03

 m1qn3: descent direction d: angle(-g,d) =  12.3 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 6, simul 6, f= 2.93708145D+01, h'(0)=-7.63349D+00

 m1qn3: line search

     mlis3       fpn=-7.633D+00 d2= 7.43D-02  tmin= 4.33D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -5.553D+00 -3.848D+00

 m1qn3: convergence rate, s(k)/s(k-1) =  8.78988D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.951D+00
     updated diagonal: average value =  1.701D-02, sqrt(variance) =  5.475D-03

 m1qn3: stopping criterion on g:  2.85948D-03

 m1qn3: descent direction d: angle(-g,d) =  16.2 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 7, simul 7, f= 2.38177322D+01, h'(0)=-5.87723D+00

 m1qn3: line search

     mlis3       fpn=-5.877D+00 d2= 1.83D-01  tmin= 2.77D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -4.701D+00 -3.786D+00

 m1qn3: convergence rate, s(k)/s(k-1) =  1.57068D+00

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  2.843D+00
     updated diagonal: average value =  7.998D-02, sqrt(variance) =  5.868D-02

 m1qn3: stopping criterion on g:  1.82670D-03

 m1qn3: descent direction d: angle(-g,d) =   7.1 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 8, simul 8, f= 1.91170296D+01, h'(0)=-1.53986D+01

 m1qn3: line search

     mlis3       fpn=-1.540D+01 d2= 2.89D+00  tmin= 6.00D-11 tmax= 1.00D+20
     mlis3     1.000D+00 -1.198D+01 -8.916D+00

 m1qn3: convergence rate, s(k)/s(k-1) =  3.96859D+00

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  3.786D+00
     updated diagonal: average value =  2.834D-01, sqrt(variance) =  1.938D-01

 m1qn3: stopping criterion on g:  1.06837D-03

 m1qn3: descent direction d: angle(-g,d) =   6.8 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 9, simul 9, f= 7.13424650D+00, h'(0)=-1.29243D+01

 m1qn3: line search

     mlis3       fpn=-1.292D+01 d2= 5.94D+00  tmin= 4.14D-11 tmax= 1.00D+20
     mlis3     1.000D+00 -7.050D+00 -1.106D+00

 m1qn3: convergence rate, s(k)/s(k-1) =  1.43412D+00

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.058D+00
     updated diagonal: average value =  2.859D-01, sqrt(variance) =  1.834D-01

 m1qn3: stopping criterion on g:  1.27599D-04

 m1qn3: descent direction d: angle(-g,d) =  40.2 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 10, simul 10, f= 8.47318698D-02, h'(0)=-1.11118D-01

 m1qn3: line search

     mlis3       fpn=-1.111D-01 d2= 5.20D-02  tmin= 4.39D-10 tmax= 1.00D+20
     mlis3     1.000D+00 -5.890D-02 -6.680D-03

 m1qn3: convergence rate, s(k)/s(k-1) =  9.35907D-02

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.068D+00
     updated diagonal: average value =  3.041D-01, sqrt(variance) =  1.940D-01

 m1qn3: stopping criterion on g:  9.05234D-05

 m1qn3: descent direction d: angle(-g,d) =  34.6 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 11, simul 11, f= 2.58275399D-02, h'(0)=-9.84247D-03

 m1qn3: line search

     mlis3       fpn=-9.842D-03 d2= 6.98D-04  tmin= 5.01D-09 tmax= 1.00D+20
     mlis3     1.000D+00 -8.414D-03 -7.053D-03

 m1qn3: convergence rate, s(k)/s(k-1) =  1.15870D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.449D+00
     updated diagonal: average value =  4.049D-01, sqrt(variance) =  2.131D-01

 m1qn3: stopping criterion on g:  6.80873D-05

 m1qn3: descent direction d: angle(-g,d) =   8.4 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 12, simul 12, f= 1.74136823D-02, h'(0)=-2.07702D-02

 m1qn3: line search

     mlis3       fpn=-2.077D-02 d2= 3.80D-03  tmin= 1.65D-09 tmax= 1.00D+20
     mlis3     1.000D+00 -1.315D-02 -7.070D-03

 m1qn3: convergence rate, s(k)/s(k-1) =  2.33420D+00

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.402D+00
     updated diagonal: average value =  5.608D-01, sqrt(variance) =  2.854D-01

 m1qn3: stopping criterion on g:  2.37288D-05

 m1qn3: descent direction d: angle(-g,d) =   1.6 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 13, simul 13, f= 4.26478840D-03, h'(0)=-3.83800D-03

 m1qn3: line search

     mlis3       fpn=-3.838D-03 d2= 1.05D-03  tmin= 3.10D-09 tmax= 1.00D+20
     mlis3     1.000D+00 -2.726D-03 -1.785D-03

 m1qn3: convergence rate, s(k)/s(k-1) =  5.24661D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.817D+00
     updated diagonal: average value =  1.007D+00, sqrt(variance) =  4.933D-01

 m1qn3: stopping criterion on g:  1.10527D-05

 m1qn3: descent direction d: angle(-g,d) =   0.2 degrees

 -------------------------------------------------------------------------------

 m1qn3: iter 14, simul 14, f= 1.53846171D-03, h'(0)=-1.55780D-03

 m1qn3: line search

     mlis3       fpn=-1.558D-03 d2= 7.95D-04  tmin= 3.55D-09 tmax= 1.00D+20
     mlis3     1.000D+00 -1.060D-03 -6.483D-04

 m1qn3: convergence rate, s(k)/s(k-1) =  8.71075D-01

 m1qn3: matrix update:
     fitting the ellipsoid: factor =  1.686D+00
     updated diagonal: average value =  1.688D+00, sqrt(variance) =  8.112D-01

 m1qn3: stopping criterion on g:  4.60536D-06

 -------------------------------------------------------------------------------

 m1qn3: output mode is  1
     number of iterations:             14
     number of simulations:            15
     realized relative precision on g:  4.61D-06
     f             =  4.78867704D-04
     dfn-norm of g =  2.30268456D-02
