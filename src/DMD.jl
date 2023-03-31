using .ControlProj1, LinearAlgebra

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

function DMD_SVD(x_steps, u_steps)

    X_k = x_steps[:,1:end-1]
    X_k1 = x_steps[:,2:end]

    n,_ = size(X_k)

    Υ_k = u_steps[1:end-1]

    Ω = vcat(X_k, Υ_k')

    U,D,V = svd(Ω)

    D = Diagonal(D)

    U1 = U[:,1:n]
    U2 = U[:,n+1:end]

    U_k1, D_k1, V_k1 = svd(X_k1)

    A_hat = U_k1' * X_k1 * V * inv(D) * U1' * U_k1
    B_hat = U_k1' * X_k1 * V * inv(D) * U2'

    return A_hat, B_hat
    
end

function error_rate(;x_steps, u_steps, A_hat, B_hat)

    X_k = x_steps[:,1:end-1]
    X_k1 = x_steps[:,2:end]

    Υ_k = u_steps[1:end-1]

    X_hat = A_hat * X_k + B_hat * Υ_k

    err = X_k1 - X_hat

    return err
    
end