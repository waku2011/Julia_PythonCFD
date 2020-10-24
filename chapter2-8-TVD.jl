using Plots
gr()

function init(q1, q2, dx, jmax, sp::Float64=0.5)
  # grid generation
  x = range(0, stop=dx*(jmax-1), length=jmax)
  # initialize q using q1, q2
  q = Array{Float64}(undef, jmax)
  for j = 1:jmax
    if j < sp*jmax
      q[j] = q1
    else
      q[j] = q2
    end
  end
  return(x, q)
end

function UPWIND1(q, c, dt, dx, j)  
    ur = q[j+1]
    ul = q[j]
    fr = c * ur # 式(2.40a)
    fl = c * ul # 式(2.40a)
    return 0.5 * (fr + fl - abs(c) * (ur - ul)) # 式(2.39)
end

function TVD(q, c, delta, g, dt, dx, j)
    sigma = 0.5 * (abs(c) - dt / dx * c^2)
    gamma = sigma * (g[j+1] - g[j]) * delta[j] / (delta[j]^2 + 1e-16)
    phi = sigma * (g[j] + g[j + 1]) - abs(c + gamma) * delta[j]
    ur = q[j+1]
    ul = q[j]
    fr = c * ur
    fl = c * ul
    return 0.5 * (fr + fl + phi)
end

function minmod(x, y)
    sgn = sign(x)
    return sgn * max(min(abs(x), sgn * y), 0.0)
end

function do_computing(x ,q, c, dt, dx, nmax, ff, order::Int64=1, interval::Int64=2)
  # plot of initial q
  plt = plot(x, q, marker=:circle, title="$ff", label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  jmax = length(x)
  flux = zeros(Float64, jmax)
  for n = 1:nmax
    qold = q
    for j = 1:jmax-1
       flux[j] = ff(qold, c, dt, dx, j)
    end
    for j= 2:jmax-1
      q[j] = qold[j] - dt / dx * (flux[j] - flux[j-1])
    end
    # 各ステップの可視化
    if n % interval == 0
      plt = plot!(x, q, marker=:circle, title="$ff", label="n=$n")
    end
  end 
    display(plt)
  readline()
end

function do_computing2(x ,q, c, dt, dx, nmax, ff, interval::Int64=2)
  # plot of initial q
  plt = plot(x, q, marker=:circle, title="$ff", label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  jmax = length(x)
  delta = zeros(Float64, jmax)
  g     = zeros(Float64, jmax)
  flux  = zeros(Float64, jmax)
  for n = 1:nmax
    qold = q
    for j = 1:jmax-1
       delta[j] = qold[j+1]-qold[j]
    end
    for j =2:jmax-1
      g[j] = minmod(delta[j], delta[j-1])
    end
    for j = 1:jmax-1
       flux[j] = ff(qold, c,delta, g,  dt, dx, j)
    end
    for j= 2:jmax-1
      q[j] = qold[j] - dt / dx * (flux[j] - flux[j-1])
    end
    # 各ステップの可視化
    if n % interval == 0
      plt = plot!(x, q, marker=:circle, title="$ff", label="n=$n")
    end
  end 
    display(plt)
  readline()
end

# UPWIND1 
c = 1
dx = 0.1
dt = 0.05
jmax = 21
nmax = 20
q1 = 1
q2 = 0
sp = 0.05
x, q = init(q1, q2, dx, jmax, sp)
order = 1
interval = 4
do_computing(x, q, c, dt, dx, nmax, UPWIND1, order, interval)

# Yee & Hartenのnon-MUSCL型TVD
c = 1
dt = 0.05
dx = 0.1
jmax = 21
nmax = 20
q1 = 1
q2 = 0
sp = 0.05
x, q = init(q1, q2, dx, jmax, sp)
interval = 4
do_computing2(x, q, c, dt, dx, nmax, TVD, interval)

