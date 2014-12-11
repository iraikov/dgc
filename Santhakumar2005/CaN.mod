TITLE CaN
 
UNITS {
    (mA) =(milliamp)
    (mV) =(millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    FARADAY = 96520 (coul)
    R = 8.3134	(joule/degC)
}
 
NEURON { 
    SUFFIX CaN
    USEION nca READ enca WRITE inca VALENCE 2 
    RANGE  gnca
    RANGE gncabar
    RANGE cinf, ctau, dinf, dtau, inca
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
        v (mV) 
        celsius (degC)
	gncabar (mho/cm2)
}
 
STATE {
	c d
}
 
ASSIGNED {
    gnca (mho/cm2)
    inca (mA/cm2)
    enca (mV)
    
    cinf dinf
    ctau (ms) 
    dtau (ms) 
    
    q10
} 

BREAKPOINT {
    SOLVE states METHOD cnexp
    gnca = gncabar*c*c*d
    inca = gnca*(v-enca)
}
 
UNITSOFF
 
INITIAL {
    q10 = 3^((celsius - 6.3)/10)
    rates(v)
    c = cinf
    d = dinf
}

DERIVATIVE states {	:Computes state variables m, h, and n 
    rates(v)	:      at the current v and dt.
    c' = (cinf-c) / ctau
    d' = (dinf-d) / dtau
}
 

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
                      
      LOCAL  alpha, beta, sum
      :"c" NCa activation system
      alpha = -0.19*vtrap(v-19.88,-10)
      beta = 0.046*exp(-v/20.73)
      sum = alpha+beta        
      ctau = 1/sum * q10
      cinf = alpha/sum
      :"d" NCa inactivation system
      alpha = 0.00016/exp(-v/48.4)
      beta = 1/(exp((-v+39)/10)+1)
      sum = alpha+beta        
      dtau = 1/sum * q10
      dinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

