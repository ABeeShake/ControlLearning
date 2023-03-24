using .ControlProj1, LinearAlgebra

function get_trajectory(A,B,x0,u0,Δt, steps)
    
    x_steps = zeros(length(x0),steps)
    x_steps[:,1] = x0

    u_steps = zeros(length(u0),steps)
    u_steps[:,1] = u0

    for i in 2:steps

        x_steps[:,i], 
        u_steps[:,i] = linsys(A,B,x_steps[:,i-1],u_steps[:,i-1],Δt)

    end

    return x_steps, u_steps

end

function DMD_SVD(x_steps, u_steps)

    X_k = x_steps[:,1:end-1]
    X_k1 = x_steps[:,2:end]

    n,_ = size(X_k)

    Υ_k = u_steps[:,1:end-1]

    Ω = vcat(X_k, Υ_k)

    U,D,V = svd(Ω)

    D = Diagonal(D)

    G_hat = X_k1 * V * inv(D) * U'

    A_hat = G_hat[:,1:n]
    B_hat = G_hat[:,n+1:end]

    return A_hat, B_hat
    
end