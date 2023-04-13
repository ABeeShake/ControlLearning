using .ControlProj1, LinearAlgebra,Statistics,Distributions, Plots

function trial(;A,B,Δt,steps)

    x0 = ones(3)

    u0 = rand(Uniform(-3,3))

    x_steps,u_steps = get_trajectory(A,B,x0,u0,Δt, 2*steps)

    x_train = x_steps[:,1:steps]
    u_train = u_steps[1:steps]

    x_test = x_steps[:,steps+1:end]
    u_test = u_steps[steps+1:end]

    A_hat,B_hat = DMD_SVD(x_steps=x_train,u_steps = u_train)

    err = error_rate(x_steps=x_train,u_steps=u_train,A_hat = A_hat, B_hat = B_hat)

    x_pred = trajectory_prediction(A_hat=A_hat,B_hat=B_hat, x_test=x_test, u_test=u_test, Δt=Δt, steps=steps)

    return err, x_test, x_pred

end

function main(;trials)

    A = [1 3 3;1 0 0;0 1 0]

    B = [1;0;0]

    Δt = 0.2
    steps = 100

    err_list = zeros(trials)

    e = zeros(steps)

    x_test = zeros(3,steps)
    x_pred = zeros(3,steps)

    for i in 1:trials

        e,x_test, x_pred = trial(A=A,B=B,Δt=Δt,steps=steps)

        err_list[i] = norm(e)

    end

    println("Average Training Error: $(mean(err_list))")

    plt1 = plot(norm.(eachcol(e)))

    display(plt1)

    plt2 = path3d(x_test[1,:],x_test[2,:],x_test[3,:],label = "true path")
    path3d!(x_pred[1,:],x_pred[2,:],x_pred[3,:], label = "prediction")

    display(plt2)

end

main(trials = 100)