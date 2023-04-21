#Packages
using Plots

#Initalizations

##system states
Qsto1 = 10
Qsto2 = 20
Qgut = 30
Gp = 75.18
Gt = 106.16
Gs = 75.18
Isc1 = 0
Il = 4.57
Ip = 25.56
X = 0
I_prime = 25.56	
XL = 25.56	
Isc2 = 0

dose = 1
conta = 0
dosekempt = 50

##hyperparameters
###non-diabetic values
ke1 = 5e-4
ke2 = 339
kp1 = 2.7
kp2 = 2.1e-3
kp3 = 9e-3
kp4 = 6.18e-2
ki = 7.9e-3
Fsnc = 1
Vm0 = 2.5
Vmx = 4.7e-2
km0 = 225.59
p2U = 3.31e-2
Vg = 1.88
k1 = 6.5e-2
k2 = 7.9e-2
Vi = 5e-2
m1 = 0.19
m2 = 0.484
m4 = 0.194
m5 = 3.04e-2
m6 = 0.6471
HEb = 0.6
m3 = (HEb * m1)/(1-HEb)
kmax = 5.58e-2
kmin = 8e-3
kabs = 5.7e-2
f = 0.9
a = 1.3e-4
b = 0.82
d = 2.36e-3
ka1 = 1.8e-3
ka2 = 1.82e-2
kd = 1.64e-2
kgut = 3.78e-2
BW = 150

#PID parameters
I_dir = 0.4 #units/kg/day Daily Insulin Requirement
Kp = I_dir/135
TI = 450 #150 (nighttime)
TD = 90 #60 (nighttime)
G_target = 90 #110 (nighttime)
#K1 = 1.966308 #unused in reference code
#K2 = 0.966584 #unused in reference code


#Renal Extraction
function renal_extraction(y)
    
    if y[4] > ke2

        escr = ke1 * (y[4] - k2)
    
    else
    
        escr = 0
    
    end

    return escr

end

#Plasma Glucose

function plasma_glucose(y,∂y,escr)

    ∂y[4] = max((kp1 - kp2*y[4]-kp3*y[12]),0) + max((kabs*y[3]*f/BW),0) - Fsnc - k1*y[4] + k2*y[5] - escr
    ∂y[5] = k1*y[4] - k2*y[5] - (Vm0 + Vmx*y[10])/(km0 + y[5])*y[5]

    return ∂y

end

#Absorption

function absorption(y,∂y)
    
    ∂y[1] = -kmax*y[1]
    Qsto = y[1] + y[2]

    if (dose == zeros(5,1)) || conta == 1

        ∂y[2] = 0
        ∂y[3] = 0

    else

        aa = 5/2/(1-b)/dosekempt
        cc = 5/2/d/dosekempt
        kgut = kmin + (kmax - kmin)/2*(tanh(aa*(Qsto-b*dosekempt)) - tanh(cc*(Qsto - d*dosekempt))+2)
        ∂y[2] = kmax*y[1] - y[2]*kgut
        ∂y[3] = kgut*y[2] - kabs*y[3]
    
    end

    return ∂y

end

#Subcutaneous glucose kinetics

function subcutaneous_glucose(y,∂y)
    
    ∂y[6] = -0.1*y[6] + 0.1*y[4]

    return ∂y

end

#Plasma Insulin and PID

function plasma_insulinPID(y,∂y)
    P = Kp*(y[6]/Vg - G_target)
    D = Kp*TD*∂y[6]/Vg
    ∂y[13] = P/TI
    PID = (P + D + y[13])/BW

    if PID > 0

        ∂y[7] = PID - (kd + ka1)*y[7]

    else

        ∂y[7] = -(kd + ka1)*y[7]

    end

    ∂y[14] = kd*y[7] - ka2*y[14]

    return ∂y

end

#Plasma Insulin

function plasma_insulin(y,∂y)
    
    ∂y[8] = -(m1 + m3)*y[8] + m2*y[9]
    ∂y[9] = m1*y[8] - (m2 + m4)*y[9] + ka1*y[7] + ka2*y[14]

    return ∂y

end

#Insulin Action

function insulin_action(y,∂y)
    
    ∂y[10] = -p2U*y[10] + p2U*(y[9]/Vi-25.56)

    return ∂y

end

#Production

function production(y,∂y)
    
    ∂y[11] = ki*(y[9]/Vi-y[11])
    ∂y[12] = ki*(y[11]-y[12])

    return ∂y

end

#Run Simulation

function run_sim(y,∂y, trials, Δt; control = true)

    Gp_list = zeros(trials)

    Ip_list = zeros(trials)

    for i in 1:trials

        Gp_list[i] = y[4]
        Ip_list[i] = y[9]

        escr = renal_extraction(y)
        ∂y = plasma_glucose(y,∂y,escr)
        ∂y = absorption(y,∂y)
        ∂y = subcutaneous_glucose(y,∂y)
        
        if control == true
        
            ∂y = plasma_insulinPID(y,∂y)
        
        end
        
        ∂y = plasma_insulin(y,∂y)
        ∂y = insulin_action(y,∂y)
        ∂y = production(y,∂y)

        y = Δt*∂y + y
        
    end

    title_text = control == 1 ? "Controlled" : "Uncontrolled"

    plt1 = plot(Gp_list)
    title!("$(title_text) Plasma Glucose Dynamics")
    xlabel!("k")
    ylabel!("Glucose (mg/kg)")
    display(plt1)

    plt2 = plot(Ip_list)
    title!("$(title_text) Plasma Insulin Dynamics")
    xlabel!("k")
    ylabel!("Insulin (pmol/kg)")
    display(plt2)
    
    return y, ∂y

end

function main()

    y = [Qsto1;Qsto2;Qgut;Gp;Gt;Gs;Isc1;Il;Ip;X;I_prime;XL;0;Isc2]
    ∂y = zeros(size(y))

    y_end, ∂y_end = run_sim(y,∂y,100000,1, control = false)

    println(y_end)
    println(∂y_end)

end

main()