: Fast delayed rectifier

NEURON {
	SUFFIX fKDR
	USEION k READ ek WRITE ik
	RANGE gbar, g
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
} 

PARAMETER {
	v0 = 23       (mV)
	taumult = 1
	gbar = 0   		(mho/cm2)
        Q10 = 3 (1)

} 


ASSIGNED {
	v 			(mV)
	ek 			(mV)
	ik			(mA/cm2)
	g				(mho/cm2)
	malpha	(/ms)		
	mbeta		(/ms)
	minf
	mtau 		(ms)
        celsius (degC)

}

STATE { m }

INITIAL { 
	rates(v)
	m = minf
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = gbar*m^4
  ik = g*(v - ek)
} 

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
}

PROCEDURE rates(v (mV)) { LOCAL tcorr
    
  tcorr = Q10^((celsius - 25(degC))/10 (degC))

  malpha = tcorr*0.07*(v+v0)/(1-exp(-0.166*(v+v0)))
  mbeta = tcorr*0.264*exp(-0.025*(v+48))
  mtau = taumult/(malpha + mbeta)
  minf = malpha/(malpha + mbeta)
}