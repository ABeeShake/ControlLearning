using .ControlProj1, LinearAlgebra,Statistics

function trial(;n,p,Δt,steps)

    A,B = gen_params(n,p)

    x0 = ones(n)

    u0 = copy(x0)

    x_steps,u_steps = get_trajectory(A,B,x0,u0,Δt, steps)

    A_hat,B_hat = DMD_SVD(x_steps,u_steps)

    A_err = norm(A - A_hat)

    B_err = norm(B - B_hat)

    return A_err, B_err

end

function main(;trials)

    p = 3
    n = 3
    Δt = 0.2
    steps = 100

    A_list = zeros(trials)
    B_list = zeros(trials)

    for i in 1:trials

        A_list[i],
        B_list[i] = trial(p=p,n=n,Δt=Δt,steps=steps)

    end

    println("Average A error: $(mean(A_list))")
    println("Average B error: $(mean(B_list))")

end

main(trials = 100)