using .ControlProj1, LinearAlgebra,Statistics,Distributions

function trial(;n,p,Δt,steps)

    A = [1 3 3;1 0 0;0 1 0]

    B = [1;0;0]

    x0 = ones(n)

    u0 = rand(Uniform(-3,3))

    x_steps,u_steps = get_trajectory(A,B,x0,u0,Δt, steps)

    A_hat,B_hat = DMD_SVD(x_steps,u_steps)

    err = error_rate(x_steps=x_steps,u_steps=u_steps,A_hat = A_hat, B_hat = B_hat)

    return mean(err)

end

function main(;trials)

    p = 3
    n = 3
    Δt = 0.2
    steps = 100

    err_list = zeros(trials)

    for i in 1:trials

        err_list[i] = trial(p=p,n=n,Δt=Δt,steps=steps)

    end

    println("Average error: $(mean(err_list))")

end

main(trials = 100)