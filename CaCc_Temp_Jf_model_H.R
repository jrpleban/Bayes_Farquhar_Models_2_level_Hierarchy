CaCc_Temp_Jf <- "model {
for (i in 1:N)
{
An[i] ~ dnorm( mu.A[i] , tau )
mu.A[i] <- min(Ac[i], Aj[i])

### Temp Depedencies on Pars (Arrhenius Temp Functions)###
#  Arrhenius Temp function 
#                                          ( Ee (Tobs - 298))
##     Y = f(Y25, Ee, Tobs) = Y25 * exp (-----------------------)
#                                          ( 298 * R * Tobs)

gammaS[i] <- gammaS25[ID[i]]*exp((T[i]-Tref)*EgammaS/(Kref*R*K[i]))
Rd[i] <- Rd25[ID[i]]* exp((T[i]-Tref)*ERd/(Kref*R*K[i]))
Kc[i] <- Kc25*exp((T[i]-Tref)*EKc/(Kref*R*K[i]))
Ko[i] <- Ko25*exp((T[i]-Tref)*EKo/(Kref*R*K[i]))
Vcmax[i] <- Vcmax25[ID[i]]*exp((T[i]-Tref)*EVcmax/(Kref*R*K[i]))


#quadratic solution for net A if limited by Rubisco
a1[i]<-(-1/g[i])
b1[i]<-((Vcmax[i]-Rd[i])/g[i])+(Ca[i]+(Kc[i]*((1+O[i])/Ko[i] )))
c1[i]<-Rd[i]*(Ca[i]+(Kc[i]*((1+O[i])/Ko[i] )))-Vcmax[i]*(Ca[i]-gammaS[i])
bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])

# quadratic solution for net A if limited by light (RuBP regeneration)
a2[i]<-(-4.0/g[i])
b2[i]<-(4.0*Ca[i]) + (8.0*gammaS[i]) + (Jf[i]/g[i]) - (4.0*Rd[i] /g[i])
c2[i]<-Jf[i] *(gammaS[i] - Ca[i] ) + (Rd[i] *(4.0*Ca[i] + 8.0*gammaS[i]))
bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])


}
for (k in 1:NID){
    ### nromal prior from Agricultural Species WULLSCHLEGER
    Vcmax25[k] ~ dnorm(mu.Vcmax25,tau.Vcmax25)T(0,)
    ### normal prior Dark Respiration from LR data (umol m-2 s-1)
    Rd25[k] ~ dnorm(mu.Rd25, tau.Rd25)
    ### prior on gamma star (Pa)
    gammaS25[k] ~ dnorm(mu.gammaS25,tau.gammaS25)T(0,)
    #Kc[k] ~ dnorm(mu.Kc, tau.Kc)
    #Ko[k] ~ dnorm(mu.Ko, tau.Ko)

}
### Activation Energy Prior(s) for Arrhenius Temp Functions
#### dnorm ( mean, precision)  precision = 1/(sd^2)

EgammaS ~ dnorm(26.8, 0.5) # weakly informed prior
ERd ~ dnorm(63.9, 0.1) # weakly informed prior
EKc ~ dnorm(70.4,0.5)  # narrow prior (evolutionary contraint)
EKo ~ dnorm(36.0,0.5)  # narrow prior  (evolutionary contraint)
EVcmax ~ dnorm(65.4, 0.01) # weakly informed prior (greater variance than other E's)

### Maximum rate Carboxylation -- nromal prior from Agricultural Species WULLSCHLEGER
mu.Vcmax25 ~ dnorm(90, 0.000625) ## SD of 40
tau.Vcmax25 ~ dgamma(5,5000)  # as variance distribution has mean = 34 ; SD = 9

### normal prior Dark Respiration from LR data (umol m-2 s-1)
mu.Rd25 ~ dnorm(3.97, 1)
tau.Rd25 ~ dgamma(9, 16)  # as variance distribution has mean = 1.4 ; SD = .25

### prior on gamma star (Pa) -- CO2 compensation point with photorespiration
mu.gammaS25 ~ dnorm(3.86,10)  ### provides narrow prior based as thought highly conserver
tau.gammaS25 ~ dgamma(50,5 )  ### mean is centered around 10 (or .3 as variance)

## Michaelis-Menten constant for CO2
Kc25 ~ dnorm(32.675,0.5)  # narrow prior (SD of 1.4) (RUBSICO evolutionary contraints)
#Michaelis-Menten constant for O2
Ko25 ~ dnorm(28612.4282, 1e-05)  # narrow prior (SD or 316) (RUBSICO evolutionary contraints)

### truncated T distribution with varaince set at 2.5 to reflect reasonalbe measurement error
tau ~ dt(0 , pow(2.5,-2), 1)T(0,)
}
"

