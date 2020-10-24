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

function UPWIND1(q, c, dt, dx, j)  
    ur = q[j + 1]
    ul = q[j]
    fr = c * ur # 式(2.40a)
    fl = c * ul # 式(2.40a)
    return 0.5 * (fr + fl - abs(c) * (ur - ul)) # 式(2.39)
end

function UPWIND2(q, c, dt, dx, j)  
    ur = 1.5 * q[j+1] - 0.5 * q[j+2]
    ul = 1.5 * q[j] - 0.5 * q[j-1]
    fr = c * ur # 式(2.40b)
    fl = c * ul # 式(2.40b)
    return 0.5 * (fr + fl - abs(c) * (ur - ul)) # 式(2.39)
end

function do_computing(x ,q, c, dt, dx, nmax, ff, order, interval)
  # plot of initial q
  plt = plot(x, q, marker=:circle, title="$ff", label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qold = q
    for j = order+1:jmax-order
      ff1 = ff(qold, c, dt, dx, j)
      ff2 = ff(qold, c, dt, dx, j-1)
      q[j] = qold[j] - dt / dx * (ff1 - ff2)
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
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
order=1
interval=2
do_computing(x, q, c, dt, dx, nmax, UPWIND1, order, interval)

c = -1    
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
order=1
interval=2
do_computing(x, q, c, dt, dx, nmax, UPWIND1, order, interval)

# UPWIND2
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
order=2
interval=2
do_computing(x, q, c, dt, dx, nmax, UPWIND2, order, interval)

c = -1
q1, q2 = 0, 1
x, q = init(q1, q2, dx, jmax)
order=2
interval=2
do_computing(x, q, c, dt, dx, nmax, UPWIND2, order, interval)

