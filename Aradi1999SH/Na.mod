: Na channel

NEURON {
	SUFFIX Na
	USEION na READ ena WRITE ina
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
	TEMP = 25 (degC)
}

PARAMETER {
	gmax = .001	(S/cm2)
}

ASSIGNED {
    v	(mV)
    ena (mV)
    ina		(mA/cm2)
    gna (S/cm2) am (/ms) bm (/ms) ah (/ms) bh (/ms)
}

STATE { m h }


BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gmax * m^3 * h
    ina = gna * (v - ena)
}


DERIVATIVE states {
    rates(v)
    
    m' = am*(1 - m) - bm*m
    h' = ah*(1 - h) - bh*h
}


INITIAL {
    rates(v)
    m = am/(am+bm)
    h = ah/(ah+bh)
}


PROCEDURE rates (v) {
    
    LOCAL amx, bmx
    
    amx = 0.2*(v + 45)
    am = 2.5*amx/(1 - exp(-amx))
    
    bmx = -0.2*(v + 17) 
    bm = 1.5*bmx/(1 - exp(-bmx)) 
    
    ah = 0.23*exp(-0.05*(v + 67))
    bh = 3.33/(1 + exp(0.1*(-14.5 - v)))
    

}
