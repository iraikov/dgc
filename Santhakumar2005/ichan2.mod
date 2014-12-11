TITLE ichan2.mod  
 
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
    SUFFIX ichan2 
    USEION nat READ enat WRITE inat VALENCE 1
    USEION kf READ ekf WRITE ikf  VALENCE 1
    USEION ks READ eks WRITE iks  VALENCE 1
    NONSPECIFIC_CURRENT il 
    RANGE  gnat, gkf, gks
    RANGE gnatbar, gkfbar, gksbar
    RANGE gl, el
    RANGE minf, mtau, hinf, htau, nfinf, nftau, inat, ikf, nsinf, nstau, iks
}
 
PARAMETER {
        v (mV) 
        celsius (degC)
        enat  (mV)
	gnatbar (mho/cm2)   
        ekf  (mV)
	gkfbar (mho/cm2)
        eks  (mV)
	gksbar (mho/cm2)
	gl (mho/cm2)    
 	el (mV)
}
 
STATE {
	m h nf ns
}
 
ASSIGNED {
         
        gnat (mho/cm2) 
        gkf (mho/cm2)
        gks (mho/cm2)

        inat (mA/cm2)
        ikf (mA/cm2)
        iks (mA/cm2)

	il (mA/cm2)

	minf hinf nfinf nsinf
 	mtau (ms) htau (ms) nftau (ms) nstau (ms)
        
        q10
} 

BREAKPOINT {
	SOLVE states METHOD cnexp
        gnat = gnatbar*m*m*m*h  
        inat = gnat*(v - enat)
        gkf = gkfbar*nf*nf*nf*nf
        ikf = gkf*(v-ekf)
        gks = gksbar*ns*ns*ns*ns
        iks = gks*(v-eks)

	il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
        q10 = 3^((celsius - 6.3)/10)
        
	rates(v)
	
	m = minf
	h = hinf
        nf = nfinf
        ns = nsinf
}

DERIVATIVE states {	:Computes state variables m, h, and n 
        rates(v)	:      at the current v and dt.
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        nf' = (nfinf-nf)/nftau
        ns' = (nsinf-ns)/nstau
}
 

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
                :"m" sodium activation system - act and inact cross at -40
	alpha = -0.3*vtrap((v+60-17),-5)
	beta = 0.3*vtrap((v+60-45),5)
	sum = alpha+beta        
	mtau = 1/sum * q10
        minf = alpha/sum
                :"h" sodium inactivation system
	alpha = 0.23/exp((v+60+5)/20)
	beta = 3.33/(1+exp((v+60-47.5)/-10))
	sum = alpha+beta
	htau = 1/sum * q10
        hinf = alpha/sum 
             :"ns" sKDR activation system
        alpha = -0.028*vtrap((v+65-35),-6)
	beta = 0.1056/exp((v+65-10)/40)
	sum = alpha+beta        
	nstau = 1/sum * q10
        nsinf = alpha/sum
            :"nf" fKDR activation system
        alpha = -0.07*vtrap((v+65-47),-6)
	beta = 0.264/exp((v+65-22)/40)
	sum = alpha+beta        
	nftau = 1/sum * q10
        nfinf = alpha/sum
	
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

