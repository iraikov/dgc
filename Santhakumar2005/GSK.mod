TITLE Calcium-activated potassium channel (non-voltage-dependent)


UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
}

NEURON {
	SUFFIX GSK
	USEION sk READ esk WRITE isk VALENCE 1
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	RANGE gskbar, qinf, qtau, isk
}

PARAMETER {
	celsius (degC)
	v	(mV)
	gskbar  (mho/cm2)
	esk	(mV)
}

STATE { q }

ASSIGNED {
	isk (mA/cm2) gsk (mho/cm2) qinf qtau (ms) q10
	ncai (mM) lcai (mM) tcai (mM)
    }
    
BREAKPOINT {          :Computes i=g*q^2*(v-esk)
	SOLVE state METHOD cnexp
        gsk = gskbar * q * q
	isk = gsk * (v-esk)
}

UNITSOFF

INITIAL {
	q10 = 3^((celsius - 6.3)/10)
	rate(ncai + lcai + tcai	)
	q=qinf
}


DERIVATIVE state {  :Computes state variable q at current v and dt.
	rate(ncai + lcai + tcai)
	q' = (qinf-q)/qtau
}

PROCEDURE rate(cac) {  :Computes rate and other constants at current v.
    
    LOCAL alpha, beta
        
	:"q" activation system
        alpha = 1.25e1 * cac * cac
        beta = 0.00025 

	qtau = 1 / (alpha + beta) * q10
	qinf = alpha * qtau
}

UNITSON
