using .ControlProj1, LinearAlgebra,Statistics,Distributions, Plots

function trial(;Δt,steps)

    A = [1 3 3;1 0 0;0 1 0]

    B = [1;0;0]

    x0 = ones(3)

    u0 = rand(Uniform(-3,3))

    x_steps,u_steps = get_trajectory(A,B,x0,u0,Δt, steps)

    A_hat,B_hat = DMD_SVD(x_steps=x_steps,u_steps=u_steps)

    err = error_rate(x_steps=x_steps,u_steps=u_steps,A_hat = A_hat, B_hat = B_hat)

    return err

end

function main(;trials)

    Δt = 0.2
    steps = 100

    err_list = zeros(trials)

    e = zeros(steps)

    for i in 1:trials

        e = trial(Δt=Δt,steps=steps)

        err_list[i] = norm(e)

    end

    println("Average error: $(mean(err_list))")
    
    plot(norm.(eachcol(e)))

end

main(trials = 100)