using Plots, LinearAlgebra
gr()

jmax = 101

dt    = 0.002
gamma = 1.4

PI   = 1.0
RHOI = 1.0
UI   = 0.0

PE   = 0.1
RHOE = 0.1
UE   = 0.0

xmin, xmid, xmax = 0.0, 0.5, 1.0
x = range(xmin, stop=xmax, length=jmax)

dx = (xmax - xmin) / (jmax - 1)
dtdx = dt / dx

function init()
  Q = zeros(jmax,3)

  for j=1:jmax
    if x[j] <= xmid
      Q[j,1] = RHOI 
      Q[j,2] = RHOI*UI 
      Q[j,3] = PI / (gamma - 1.0) + 0.5 * RHOI * UI^2
    else
      Q[j,1] = RHOE
      Q[j,2] = RHOE*UE 
      Q[j,3] = PE / (gamma - 1.0) + 0.5 * RHOE * UE^2
    end
  end

  return Q
end

function calc_CFL(Q)
    rho  = Q[:, 1]
    rhou = Q[:, 2]
    e    = Q[:, 3]
    u = rhou ./ rho
    p = (gamma - 1.0) * (e .- 0.5.*rho.*u.^2 )
    c = sqrt.(gamma .* p ./ rho)
    sp = c .+ abs.(u)
    return maximum(sp) * dtdx   
end

function E_flux(Q, E, jmax)
    
    rho  = Q[:, 1]
    rhou = Q[:, 2]
    e    = Q[:, 3]
    
    u = rhou ./ rho
    p = (gamma - 1.0) * (e .- 0.5.*rho.*u.^2 )
    
    for j=1:jmax
      E[j,1] = rhou[j]
      E[j,2] = p[j] + rhou[j]*u[j]
      E[j,3] = (e[j]+p[j]) * u[j]
    end
end

function Roe_flux(QL, QR, E)
    
    for j = 1:jmax-1        
        dQ = zeros(3)
        
        rhoL, uL, pL = QL[  j, 1], QL[  j, 2], QL[  j, 3]
        rhoR, uR, pR = QR[j+1, 1], QR[j+1, 2], QR[j+1, 3]
        
        rhouL = rhoL * uL
        rhouR = rhoR * uR

        eL = pL / (gamma - 1.0) + 0.5 * rhoL * uL ^ 2
        eR = pR / (gamma - 1.0) + 0.5 * rhoR * uR ^ 2

        HL = (eL + pL) / rhoL
        HR = (eR + pR) / rhoR
        
        cL = sqrt((gamma - 1.0) * (HL - 0.5 * uL ^ 2))
        cR = sqrt((gamma - 1.0) * (HR - 0.5 * uR ^ 2))
                
        # Roe平均 式(6.38)        
        sqrhoL = sqrt(rhoL)
        sqrhoR = sqrt(rhoR)

        rhoAVE = sqrhoL * sqrhoR
        uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
        HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR) 
        cAVE = sqrt((gamma - 1.0)* (HAVE - 0.5 * uAVE ^ 2))
        eAVE = rhoAVE * (HAVE - cAVE ^ 2 / gamma)
        
        dQ[1] = rhoR - rhoL
        dQ[2] = rhoR * uR - rhoL * uL
        dQ[3] = eR - eL
        
        Lambda = [abs(uAVE-cAVE)  0.0        0.0
                  0.0            abs(uAVE)   0.0
                  0.0     0.0  abs(uAVE+cAVE) ]
        b1 = 0.5 * (gamma - 1.0) * uAVE ^ 2 / cAVE ^ 2
        b2 = (gamma - 1.0) / cAVE ^ 2

        R = [               1.0          1.0          1.0
                    uAVE-cAVE            uAVE    uAVE+cAVE
             HAVE - uAVE*cAVE  0.5*uAVE^2  HAVE+uAVE*cAVE]
        
        Rinv = [0.5*(b1+uAVE/cAVE)  -0.5*(b2*uAVE+cAVE)  0.5*b2
                1.0-b1                 b2*uAVE      -b2
                0.5*(b1-uAVE/cAVE)  -0.5*(b2*uAVE-cAVE)  0.5*b2]
        
      
        AQ = R * Lambda * Rinv * dQ
              
        EL = Array([rhoL * uL, pL + rhouL * uL, (eL + pL) * uL])
        ER = Array([rhoR * uR, pR + rhouR * uR, (eR + pR) * uR])
        
        E[j,:] = 0.5 .* (ER[:] .+ EL[:] .- AQ[:]) # 式(6.43)
    end
end

function minmod(x, y)
  sgn = sign(x)
  return sgn * maximum([minimum([abs(x), sgn*y]), 0.0])
end

function MUSCL(Q, order, kappa)
    # 基本変数で内挿する
    rho, rhou, e = Q[:, 1], Q[:, 2], Q[:, 3]
    
    Q[:, 2] = rhou ./ rho  # u
    Q[:, 3] = (gamma - 1.0) .* (e .- 0.5 .* rho .* Q[:, 2] .^ 2) # p
    
    if (order == 2 || order == 3)
        # 2nd / 3rd order & minmod limitter
        dQ = zeros(jmax, 3)
        for j = 1:jmax-1
          dQ[j,:] = Q[j+1,:] - Q[j,:]
        end
        b = (3.0 - kappa) / (1.0 - kappa) # 式(2.74)
        
        Dp = zeros(jmax, 3)
        Dm = zeros(jmax, 3)
        for j = 2:jmax-1
           Dp[j,:] = minmod.(dQ[j,:]  , b .* dQ[j-1,:]) # 式(2.73a)
           Dm[j,:] = minmod.(dQ[j-1,:], b .* dQ[j,:])   # 式(2.73b)
        end
        Dp[1,:] = Dp[2,:]
        Dm[1,:] = Dm[2,:]
        
        QL = copy(Q)
        QR = copy(Q)
        for j = 2:jmax-1
            QL[j,:] += 0.25 * ((1.0 - kappa) * Dp[j,:] + (1.0 + kappa) * Dm[j,:]) # 式(2.72a)
            QR[j,:] -= 0.25 * ((1.0 + kappa) * Dp[j,:] + (1.0 - kappa) * Dm[j,:]) # 式(2.72b)
        end 
    else
        # 1st order
        QL = copy(Q)
        QR = copy(Q)
    end
    return QL, QR
end

function Roe_FDS(Q, order, kappa, nmax, interval)
    E = zeros(jmax,3)

    for n = 1:nmax
        if n % interval == 0
          println("n = ",n," CFL = ", calc_CFL(Q))
        end

        Qold = copy(Q)
                
        coefs = [0.5, 1.0]
        for coef in coefs
            QL, QR = MUSCL(Qold, order, kappa)
        
            Roe_flux(QL, QR, E)
            for j = 2: jmax-1
                Qold[j,:] = Q[j,:] - coef * dtdx * (E[j,:]- E[j-1,:])
            end  
            Qold[1,:] = Q[1,:]
            Qold[jmax,:] = Q[jmax,:]
        end 
        Q[:,:] = Qold[:,:]
    end
end

nmax = 100
print_interval = 4
order = 2

# -1  = 2nd order fully upwind, 
# 0   = 2nd order upwind biased
# 1/3 = 3rd order upwind biased
kappa = 0

Q = init()
Roe_FDS(Q, order, kappa, nmax, print_interval)

## analytical solution

Pext = zeros(jmax, 3)
Qext = zeros(jmax, 3)

GUESS   = 1.0
FINC    = 0.01
itemax1 = 5000
itemax2 = 500

CI = sqrt(gamma * PI / RHOI)
CE = sqrt(gamma * PE / RHOE)
P1P5 = PI / PE

GAMI = 1.0 / gamma
GAMF = (gamma - 1.0) / (2.0 * gamma)
GAMF2 = (gamma + 1.0) / (gamma - 1.0)
GAMFI = 1.0 / GAMF

for it1 = 1:itemax1
    for it2 = 1:itemax2
        global GUESS, FINC
        SQRT1 = (gamma - 1.0) * (CE / CI) * (GUESS - 1.0)
        SQRT2 = sqrt(2.0 * gamma * (2.0 * gamma + (gamma + 1.0) * (GUESS - 1.0)))
        FUN = GUESS * (1.0 - (SQRT1 / SQRT2)) ^ (-GAMFI)
        DIF = P1P5 - FUN
        
        if abs(DIF) <= 0.000002
            break
        end
        
        if DIF >= 0.0
            GUESS += FINC
        else
            GUESS -= FINC
            FINC = 0.5 * FINC
        end
        
    end
end

P4P5 = GUESS
P4 = PE * P4P5
P3P1 = P4P5 / P1P5
P3 = P3P1 * PI

R4R5 = (1.0 + GAMF2 * P4P5) / (GAMF2 + P4P5)
RHO4 = RHOE * R4R5
U4 = CE * (P4P5 - 1.0) * sqrt(2.0 * GAMI / ((gamma + 1.0) * P4P5 + (gamma - 1.0)))
C4 = sqrt(gamma * P4 / RHO4)

R3R1 = P3P1 ^ GAMI
RHO3 = RHOI * R3R1 
U3 = 2.0 * CI / (gamma - 1.0) * (1.0 - P3P1 ^ GAMF)
C3 = sqrt(gamma * P3 / RHO3)
CS =  CE * sqrt(0.5 * ((gamma - 1.0) * GAMI + (gamma + 1.0) * GAMI * P4 / PE))

TOT = 0.0
EPST = 1.0e-14
for n = 1:nmax
    global TOT
    TOT = TOT + dt
    rad = dt / dx
    
    x1 = xmid - CI * TOT
    x2 = xmid - (CI - 0.5 * (gamma + 1.0) * U3) * TOT
    x3 = xmid + U3 * TOT
    x4 = xmid + CS * TOT
    
    for j = 1:jmax
        xx = x[j]
        if xx <= x1
            Qext[j, 1] = RHOI
            Qext[j, 2] = RHOI * UI
            Qext[j, 3] = PI / (gamma - 1.0) + 0.5 * UI * Qext[j, 1]
            Pext[j] = PI
        elseif xx <= x2
            UT = UI + (U3 - UI) / ((x2 - x1) + EPST) * ((xx - x1) + EPST)
            RTRI = (1.0 - 0.5 * (gamma - 1.0) * UT / CI) ^ (2.0 / (gamma - 1.0))
            RT = RHOI * RTRI
            PT = RTRI ^ gamma * PI
            Qext[j, 1] = RT
            Qext[j, 2] = RT * UT
            Qext[j, 3] = PT / (gamma - 1.0) + 0.5 * UT * Qext[j, 2]
            Pext[j] = PT
        elseif xx <= x3
            Qext[j, 1] = RHO3
            Qext[j, 2] = RHO3 * U3
            Qext[j, 3] = P3 / (gamma - 1.0) + 0.5 * U3 * Qext[j, 2]
            Pext[j] = P3
        elseif xx <= x4
            Qext[j, 1] = RHO4
            Qext[j, 2] = RHO4 * U4
            Qext[j, 3] = P4 / (gamma - 1.0) + 0.5 * U4 * Qext[j, 2]
            Pext[j] = P4
        else
            Qext[j, 1] = RHOE
            Qext[j, 2] = RHOE * UE
            Qext[j, 3] = PE / (gamma - 1.0) + 0.5 * UE * Qext[j, 2]
            Pext[j] = PE            
        end
   end
end

# Visualization         
plt = plot(x, Qext[:,1], linecolor=:black, linewidth=1, linestyle =:dash, label="Analytical", xlabel="x", ylabel="rho")
plt = plot!(x, Q[:,1]   , linecolor=:red, linewidth=2, label ="Numerical")
display(plt)
readline()
println("Press enter ...")

plt = plot(x, Qext[:,2]./Qext[:,1], linecolor=:black, linewidth=1, linestyle =:dash, label="Analytical", xlabel="x", ylabel="u")
plt = plot!(x, Q[:,2]./Q[:,1]   , linecolor=:red, linewidth=2, label ="Numerical")
display(plt)
readline()
println("Press enter ...")

yext = (gamma - 1.0) .* (Qext[:,3] .- 0.5 * Qext[:,2] .^2 ./ Qext[:,1])
y    = (gamma - 1.0) .* (Q[:,3]    .- 0.5 *    Q[:,2] .^2 ./ Q[:,1])
plt = plot(x, yext, linecolor=:black, linewidth=1, linestyle =:dash, label="Analytical", xlabel="x", ylabel="p")
plt = plot!(x, y   , linecolor=:red, linewidth=2, label ="Numerical")
display(plt)
readline()
println("Press enter ...")
