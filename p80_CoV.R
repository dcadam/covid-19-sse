#proportion resposible formula
propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}

#Empirical Estiamtes P80% Total
propresponsible(0.58,0.43,0.8)
propresponsible(0.45,0.29,0.8)
propresponsible(0.72,0.67,0.8)

#Branching Process p80% Total
propresponsible(0.74,0.33,0.8)
propresponsible(0.58,0.14,0.8)
propresponsible(0.97,0.98,0.8)

#Hypothetical Scenario One P80% Total
propresponsible(0.62,0.35,0.8)
propresponsible(0.49,0.25,0.8)
propresponsible(0.70,0.56,0.8)

#Hypothetical Scenario Two P80% Total
propresponsible(0.72,0.19,0.8)
propresponsible(0.53,0.13,0.8)
propresponsible(0.94,0.26,0.8)

#Empirical Estiamtes P80% Wave 1
propresponsible(0.61,0.40,0.8)
propresponsible(0.38,0.21,0.8)
propresponsible(0.88,0.86,0.8)

#Empirical Estiamtes P80% Wave 2
propresponsible(0.57,0.44,0.8)
propresponsible(0.42,0.27,0.8)
propresponsible(0.73,0.82,0.8)

#Branching Process p80% Wave 1
propresponsible(0.66,2.31,0.8)
propresponsible(0.44,0.22,0.8)
propresponsible(0.99,1*10^100,0.8)

#Branching Process p80% Wave2
propresponsible(0.77,0.20,0.8)
propresponsible(0.54,0.08,0.8)
propresponsible(1.13,0.63,0.8)


#Coeffcient of variation formula
nbinomcoefficient=function(r,k){
  var=r*(1+(r/k))
  sd=sqrt(var)
  sd/r
}

#Empirical Estiamtes CV
nbinomcoefficient(0.58,0.43)
nbinomcoefficient(0.45,0.29)
nbinomcoefficient(0.72,0.67)

#Hypothetical Scenario One CV
nbinomcoefficient(0.62,0.35)
nbinomcoefficient(0.49,0.25)
nbinomcoefficient(0.70,0.56)

#Hypothetical Scenario Two CV
nbinomcoefficient(0.72,0.19)
nbinomcoefficient(0.53,0.13)
nbinomcoefficient(0.94,0.26)

#Empirical Estiamtes CV Wave 1
nbinomcoefficient(0.61,0.40)
nbinomcoefficient(0.38,0.21)
nbinomcoefficient(0.88,0.86)

#Empirical Estiamtes CV Wave 2
nbinomcoefficient(0.57,0.44)
nbinomcoefficient(0.42,0.27)
nbinomcoefficient(0.73,0.82)

#Branching Process CV
nbinomcoefficient(0.74,0.33)
nbinomcoefficient(0.58,0.14)
nbinomcoefficient(0.97,0.98)

#Branching Process CV Wave 1
nbinomcoefficient(0.66,2.31)
nbinomcoefficient(0.44,0.22)
nbinomcoefficient(0.99,Inf)

#Branching Process CV Wave2
nbinomcoefficient(0.77,0.21)
nbinomcoefficient(0.54,0.08)
nbinomcoefficient(1.13,0.63)

