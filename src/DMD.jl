using .ControlProj1, LinearAlgebra, TSVD

function get_trajectory(A,B,x0,u0,Δt, steps)
    
    x_steps = zeros(length(x0),steps)
    x_steps[:,1] = x0

    u_steps = zeros(steps)
    u_steps[1] = u0

    for i in 2:steps

        x_steps[:,i], 
        u_steps[i] = linsys(A = A,B = B,x = x_steps[:,i-1],u = u_steps[i-1],Δt=Δt,noise=true)

    end

    return x_steps, u_steps

end

function DMD_SVD(;x_steps, u_steps, low_rank = false)

    X_k = x_steps[:,1:end-1]
    X_k1 = x_steps[:,2:end]

    n,_ = size(X_k)

    Υ_k = u_steps[1:end-1]

    Ω = vcat(X_k, Υ_k')

    U,D,V = tsvd(Ω)

    D = Diagonal(D)

    G_hat = X_k1 * V * inv(D) * U'

    A_hat = G_hat[:,1:n]
    B_hat = G_hat[:,n+1:end]

    if low_rank == true

        U_k1, D_k1, V_k1 = tsvd(X_k1)

        A_hat = U_k1' * A_hat * U_k1
        B_hat = U_k1' * B_hat

        return A_hat, B_hat, U_k1

    end

    return A_hat, B_hat
    
end

function error_rate(;x_steps, u_steps, A_hat, B_hat, U = I)

    X_k = x_steps[:,1:end-1]
    X_k1 = x_steps[:,2:end]

    Υ_k = u_steps[1:end-1]

    X_hat = U * (A_hat * X_k + reduce(hcat,[B_hat .* u for u in Υ_k]))

    err = X_k1 - X_hat

    return err
    
end