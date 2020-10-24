using Plots
gr()

const dt = 0.05 # time step, s
const dx = 0.1  # grid size, m

const jmax = 21 
const nmax = 6  

function init(q1, q2, dx, jmax)
  # grid generation
  x = range(0, stop=dx*(jmax-1), length=jmax)

  # initialize q using q1, q2
  q = Array{Float64}(undef, jmax)
  for j = 1:jmax
    if j < jmax/2+1
      q[j] = q1
    else
      q[j] = q2
    end
  end
  return(x, q)
end

function FTCS(q, c, dt, dx, j)  
  return 0.5 * c * (q[j + 1] + q[j])  # 式(2.29)  
end

function UPWIND1(q, c, dt, dx, j)  # allow c>0 only
  return c * q[j]   # 式(2.30)  
end

function UPWIND1_mod(q, c, dt, dx, j)  
  return 0.5 * ( c * (q[j+1] + q[j]) - abs(c) * (q[j+1] - q[j]))  # 式(2.33) 
end

function LAX(q, c, dt, dx, j)
    nu2 = 1 / (c * dt / dx)
    return 0.5 * c * ((1 - nu2) * q[j+1] + (1 + nu2) * q[j])  # 式(2.31)
end

function LAXWEN(q, c, dt, dx, j)
    nu = c * dt / dx
    return 0.5 * c * ((1 - nu) * q[j+1] +(1 + nu) * q[j])  # 式(2.32)
end

function do_computing(x ,q, c, dt, dx, nmax, ff)
  # plot of initial q
  plt = plot(x, q, marker=:circle, title="$ff", label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qold = q
    for j = 2:jmax-1
      ff1 = ff(qold, c, dt, dx, j)
      ff2 = ff(qold, c, dt, dx, j-1)
      q[j] = qold[j] - dt / dx * (ff1 - ff2)
       
    end
    # 各ステップの可視化
    if n % 2 == 0
      plt = plot!(x, q, marker=:circle, title="$ff", label="n=$n")
    end
  end 
  
  display(plt)
  readline()
end

# UPWIND1 
c = 1    
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, UPWIND1)

# modified UPWIND1
c = -1
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, UPWIND1_mod)

# FTCS
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, FTCS)

c = -1
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, FTCS)

# LAX
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, LAX)
c = -1
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, LAX)

# LAX-WENDROFF
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, LAXWEN)

c = -1
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
do_computing(x, q, c, dt, dx, nmax, LAXWEN)


