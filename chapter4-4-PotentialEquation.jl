nmax = 400

M = 0.1
alpha2 = 1 - M^2
Uinf = 0.1

dx = dy = 0.05

# xs, xe = -5.0, 6.0
# ys, ye = 0.0, 5.0
xs, xe = -1.0, 2.0
ys, ye = 0.0, 1.0
x_le, x_te = 0.0, 1.0

jmax = floor(Int, (xe - xs) / dx) + 1
kmax = floor(Int, (ye - ys) / dy) + 1

j_le = floor(Int, (x_le - xs) / dx)
j_te = floor(Int, (x_te - xs) / dx) + 1

x = range(xs, stop=xe, length=jmax)
y = range(ys, stop=ye, length=kmax)

phi = zeros(jmax, kmax)
u   = zeros(jmax, kmax)
v   = zeros(jmax, kmax)
dydx = [(if j_le <= j < j_te 0.4*(1.0-2.0*x[j]) else 0.0 end) for j in 1:jmax ]

# main

residual = zeros(nmax)
for n = 1:nmax
    phiold = copy(phi)
    
    # 境界条件
    phi[1, :] = 0.0
    phi[jmax, :] = 0.0
    phi[:, kmax] = 0.0
        
    for j = 1:jmax
        phi[j, 1] = phi[j, 2] - dydx[j] * dy
    end

    # Gaus Seidel法
    for k = 2:kmax-1
        for j = 2:jmax-1
            phi[j, k] = 1.0/(2.0*alpha2+2.0)*(alpha2*(phi[j-1,k]+phi[j+1,k])+phi[j,k-1]+phi[j,k+1])
        end
    end    
    residual[n]= sum(sqrt.(((phi-phiold)^2))/(jmax * kmax))
end

for j = 2:jmax-1
    u[j, :] = Uinf * (1.0 + (phi[j + 1, :] - phi[j - 1, :]) / (2 * dx))
end
u[1,:] = Uinf * (1.0 + (phi[2, :] - phi[1, :]) / dx)
u[jmax,:] = Uinf * (1.0 + (phi[jmax, :] - phi[jmax-1, :]) / dx)
    
for k = 1:kmax-1
    v[:, k] = Uinf * (phi[:, k + 1] - phi[:, k - 1]) / (2 * dy)
end
v[:,1] = Uinf * (phi[:, 2] - phi[:, 1]) / dy
v[:,kmax] = Uinf * (phi[:, kmax] - phi[:, kmax-1]) / dy
    
va = sqrt(u^2 + v^2)    

# residual plot

