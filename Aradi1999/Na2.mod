: Na channel

NEURON {
	SUFFIX Na2
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
    printf ("t = %g gmax = %g gna = %g ena = %g m = %g h = %g\n", t, gmax, gna, ena, m, h)
    printf ("t = %g am = %g bm = %g\n", t, am, bm)
}


PROCEDURE rates (v) {
    
    am = -0.3 * (v - 25) / (exp ((v - 25) / -5) - 1)
    bm = 0.3 * (v - 53) / (exp ((v - 53) / 5) - 1)
    
    ah = 0.23 / exp ((v - 3) / 20)
    bh = 3.33 / (exp ((v - 55.5) / -10) + 1)
    

}
