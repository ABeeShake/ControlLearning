using LinearAlgebra, Statistics, Distributions, Plots

##system states
Isc1 = 0
Il = 4.57
Ip = 25.56
X = 0
I_prime = 25.56	
XL = 25.56	
Isc2 = 0

y = [Isc1;Il;Ip;X;I_prime;XL;Isc2]

##hyperparameters
###non-diabetic values
m1 = 0.19
m2 = 0.484
m4 = 0.194
m5 = 3.04e-2
m6 = 0.6471
HEb = 0.6
m3 = (HEb * m1)/(1-HEb)
ka1 = 1.8e-3
ka2 = 1.82e-2
ki = 7.9e-3
Vi = 5e-2
p2U = 3.31e-2

## Non-Noisy Dynamics

#Renal Extraction
function plasma_insulin!(y,∂y)
    
    ∂y[2] = -(m1 + m3)*y[2] + m2*y[3]
    ∂y[3] = m1*y[2] - (m2 + m4)*y[3] + ka1*y[1] + ka2*y[7]

end

#Insulin Action

function insulin_action!(y,∂y)
    
    ∂y[4] = -p2U*y[4] + p2U*(y[3]/Vi-25.56)

end

#Production

function production!(y,∂y)
    
    ∂y[5] = ki*(y[3]/Vi-y[5])
    ∂y[6] = ki*(y[5]-y[6])

end

#Generate Sequence Data

function step!(y,∂y,Δt)

    plasma_insulin!(y,∂y)
    insulin_action!(y,∂y)
    production!(y,∂y)

    y += Δt*∂y
    
end

