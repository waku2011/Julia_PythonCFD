using Plots
gr()

dt = 0.05 # time step, s
dx = 0.1  # grid size, m

jmax = 21 
nmax = 10 

function init(q1, q2, dx, jmax)
  # grid generation
  xs = -1.0
  x = range(xs, stop=xs+dx*(jmax-1), length=jmax)

  # initialize q using q1, q2
  q = Array{Float64}(undef, jmax)
  for j = 1:jmax
    if x[j] < 0.0
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
    fr = c * ur 
    fl = c * ul 
    return 0.5 * (fr + fl - abs(c) * (ur - ul))
end

function MC(q, c, dt, dx, j)  
    ur = q[j + 1]
    ul = q[j]
    fr = 0.5 * ur^2
    fl = 0.5 * ul^2
    c = 0.5 * (ur + ul)
    return 0.5 * (fr + fl - sign(c) * (fr - fl))
end

function GODUNOV(q, c, dt, dx, j)  
    qm = 0.5 * (q[j] + abs(q[j]))
    qp = 0.5 * (q[j+1] - abs(q[j+1]))
    return max(0.5*qm^2, 0.5*qp^2)
end

function do_computing(x ,q, c, dt, dx, nmax, ff, order::Int64=1, interval::Int64=2)
  # plot of initial q
  plt = plot(x, q, marker=:circle, title="$ff", label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qold = q
    for j = order+1:jmax-order
      ff1 = ff(qold, qold[j], dt, dx, j)
      ff2 = ff(qold, qold[j], dt, dx, j-1)
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

# MC  
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
nmax = 10
do_computing(x, q, c, dt, dx, nmax, MC)

# GODUNOV
c = 1
q1, q2 = 1, 0
x, q = init(q1, q2, dx, jmax)
nmax=10
do_computing(x, q, c, dt, dx, nmax, GODUNOV)

