using Plots, LinearAlgebra
gr()

function init(xs, dx, jmax)
    x = range(xs, stop=xs+dx*(jmax-1), length=jmax) 
    q = sin.(4.0 * pi * x)
    return (x, q)
end

function filter_LHS_4th()
    M = diagm([1.0 for j =1:jmax])
    M[1,jmax] = alpha_f
    M[jmax,1] = alpha_f
    for j = 1:jmax-1
       M[j,j+1] = alpha_f
    end
    for j = 2:jmax
       M[j,j-1] = alpha_f
    end
    return M    
end

function filter_RHS_4th(rhs, f, alpha)       
    a = [5/8 + 3/4 * alpha, 1/2 + alpha, -1/8 + alpha/4]
    
    rhs[1] = 0.5 * a[1] * (f[1]+f[1]) + # Eq.(9.3)
             0.5 * a[2] * (f[2]+f[jmax]) +
             0.5 * a[3] * (f[3]+f[jmax-1])
 
    rhs[2] = 0.5 * a[1] * (f[2]+f[2]) + # Eq.(9.3)
             0.5 * a[2] * (f[3]+f[1]) +
             0.5 * a[3] * (f[4]+f[jmax])   
    
    for j = 3:jmax-2
       rhs[j] = 0.5 * a[1] * (f[j  ]+f[j  ]) + # Eq.(9.3)
                0.5 * a[2] * (f[j+1]+f[j-1]) + # Eq.(9.3)
                0.5 * a[3] * (f[j+2]+f[j-2]) # Eq.(9.3)
    end
   
    rhs[jmax-1] = 0.5 * a[1] * (f[jmax-1]+f[jmax-1]) + # Eq.(9.3)
                  0.5 * a[2] * (f[jmax  ]+f[jmax-2]) +
                  0.5 * a[3] * (f[1     ]+f[jmax-3])    
       
    rhs[jmax] = 0.5 * a[1] * (f[jmax]+f[jmax]) + # Eq.(9.3)
                0.5 * a[2] * (f[1]+f[jmax-1]) +
                0.5 * a[3] * (f[2]+f[jmax-2])    
                    
end

function cmpt_LHS_4th()
    M = diagm([1.0 for j =1:jmax])
    M[1,jmax] = 1/4
    M[jmax,1] = 1/4
    for j = 1:jmax-1
       M[j,j+1] = 1/4
    end
    for j = 2:jmax
       M[j,j-1] = 1/4
    end
    return M    
end

function cmpt_RHS_4th(rhs, f)
    a = 3 / 2
    
    rhs[1]    = a / (2 * dx) * (f[2] - f[jmax])
    rhs[jmax] = a / (2 * dx) * (f[1] - f[jmax-1])
    
    for j = 2: jmax-1
        rhs[j] = a / (2 * dx) * (f[j + 1] - f[j - 1])
    end
end

function cmpt_4th(x, q, c, dt, dx, nmax, enable_filter, interval)
    
    # init 
    plt = plot(push!(Array(x),x[jmax]+dx), push!(Array(q),q[1]), marker=:circle, lw=2, label="n=0") 
    
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
        if n % interval == 0
            plt = plot!(push!(Array(x),x[jmax]+dx), push!(Array(q),q[1]), marker=:circle, lw=2, label="n=$n")
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
enable_filter = false
interval = 50
cmpt_4th(x, q, c, dt, dx, nmax, enable_filter, interval) 

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
enable_filter = true
interval = 50
cmpt_4th(x, q, c, dt, dx, nmax, enable_filter, interval) 

