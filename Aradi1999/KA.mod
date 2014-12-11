: KA channel

NEURON {
	SUFFIX KA
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
    gk (S/cm2) ak (/ms) bk (/ms) al (/ms) bl (/ms)
}

STATE { k l }


BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gmax * k * l
    ik = gk * (v - ek)
}


DERIVATIVE states {
    rates(v)
    
    k' = ak*(1 - k) - bk*k
    l' = al*(1 - l) - bl*l
}


INITIAL {
    rates(v)
    k = ak/(ak+bk)
    k = al/(al+bl)
}


PROCEDURE rates (v) {
    
    LOCAL akx, bkx
    
    akx = 0.06667*(v + 25) 
    ak = 0.75*akx/(1 - exp(-akx)) 
    
    bkx = -0.125*(v + 15) 
    bk = 0.8*bkx/(1 - exp(-bkx))
    
    al = 0.00015*exp(-0.06667*(v + 13))            
    bl = 0.06/(1 + exp(0.08333*(-68 - v)))

}
