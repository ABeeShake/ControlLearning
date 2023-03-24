using Distributions, LinearAlgebra

function linsys(A,B,x,u,Δt)

    x_next = (Δt .* A + I)*x + Δt .* B * u

    u_next = copy(x_next)

    return x_next, u_next
    
end

function gen_params(n,p)

    A = randn(n,n)
    B = randn(n,p)

    return A,B

end