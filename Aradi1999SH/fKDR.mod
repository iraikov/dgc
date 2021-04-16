: fKDR channel

NEURON {
	SUFFIX fKDR
	USEION k READ ek WRITE ik
	RANGE gmax
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	gmax = .001	(S/cm2)
}

ASSIGNED {
    v	(mV)
    ek (mV)
    ik	(mA/cm2)
    gk (S/cm2) an (/ms) bn (/ms)
}

STATE { n }


BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gmax * n^4
    ik = gk * (v - ek)
}


DERIVATIVE states {
    rates(v)
    
    n' = an*(1 - n) - bn*n
}


INITIAL {
    rates(v)
    n = an/(an+bn)
}


PROCEDURE rates (v) {
    
    LOCAL anx
    
    anx = 0.1667*(v + 23) 
    an = 0.42*anx/(1 - exp(-anx)) 
    bn = 0.264*exp(-0.025*(v + 48))
}
