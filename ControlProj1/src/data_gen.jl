using Distributions, LinearAlgebra

function linsys(;A,B,x,u,Δt, noise=true)

    if noise == true

        x_next = (Δt .* A + I)*x + Δt .* B .* u + randn(size(x))

    else

        x_next = (Δt .* A + I)*x + Δt .* B .* u

    end
    
    u_next = rand(Uniform(-3,3))

    return x_next, u_next
    
end

function gen_params(n,p)

    A = randn(n,n)
    B = randn(n,p)

    return A,B

end