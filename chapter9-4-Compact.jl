using Plots, LinearAlgebra
gr()

function init(xs, dx, jmax):
    x = range(xs, stop = xs + dx * (jmax - 1), length=jmax) 
    q = sin.(4.0 * pi * x)
    return (x, q)
end

function filter_LHS_4th()
    M = diagm([1.0 for j =1:jmax])
    for j = 1:jmax-1
       M[j,j+1] = alpha_f
    end
    for j = 2:jmax
       M[j,j-1] = alpha_f
    end
    return M    
end

function filter_RHS_4th(rhs, f, alpha):        
    N = 2
    a = [5./8 + 3./4 * alpha, 1./2 + alpha, -1./8 + alpha/4]
    for j = -N:jmax-N
        rhs[j] = 0.0
        for n = 1:3
           rhs[j] += 0.5 * a[n] * (f[j+n]+f[j-n]) # Eq.(9.3)
        end
    end
end

function cmpt_LHS_4th():
    M = diagm([1.0 for j =1:jmax])
    for j = 1:jmax-1
       M[j,j+1] = 1/4
    end
    for j = 2:jmax
       M[j,j-1] = 1/4
    end
    return M    
end

function cmpt_RHS_4th(rhs, f):
    a = 3 / 2
    
    rhs[0]  = a / (2 * dx) * (f[1] - f[-1])
    rhs[-1] = a / (2 * dx) * (f[0] - f[-2])
    
    for j = 1: jmax-1
        rhs[j] = a / (2 * dx) * (f[j + 1] - f[j - 1])
    end
end

function cmpt_4th(x, q, c, dt, dx, nmax, enable_filter = False, interval = 2):
    
    # init 
    plt = plot(push!(x, x[jmax]+dx), push!(q, q[0]), marker=:circle, lw=2, label='n=0') 
    
    L = cmpt_LHS_4th()
    Linv = inv(L)
    
    if enable_filter
        FL = filter_LHS_4th()
        FLinv = inv(FL)
    end
    
    R = zeros(jmax)
    for n = 1:nmax
        
        qold = copy(q)

        coefs = [0.5, 1.0]
        # coefs = [1/4, 1/3, 1/2, 1]
        for coef in coefs
            f = c * qold
            cmpt_RHS_4th(R, f)
            dfdx = Linv * R
            qold = q - coef * dt * dfdx
        end

        if enable_filter
            f = c * qold
            filter_RHS_4th(R, f, alpha_f)
            qold = FLinv * R
        end

        q = qold
            
        # visualization
        if n % interval == 0:
            plt = plot!(push!(x, x[jmax]+dx), push!(q, q[0]), marker=:circle, lw=2, label="n=$n")
        end
    end
    display(plt)
    readline()
end

# without filter
xs, xe = -1, 1
jmax = 20 
c = 1
dx = 0.1
dt = 0.01
println("CFL= ", c * dt / dx )

nmax = 100
x, q = init(xs, dx, jmax)
cmpt_4th(x, q, c, dt, dx, nmax, enable_filter = False, interval = 50) 

# with filter
xs, xe = -1, 1
jmax = 20 
c = 1
dx = 0.1
dt = 0.01
println("CFL= ", c * dt / dx )

alpha_f = 0.49

nmax = 100
x, q = init(xs, dx, jmax)
cmpt_4th(x, q, c, dt, dx, nmax, enable_filter = True, interval = 50) 

)
