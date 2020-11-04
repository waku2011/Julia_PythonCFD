using Plots; gr()

const c = 1.0   # advection speed

const dt = 0.05 # time step
const dx = 0.1  # grid size

const jmax = 21 
const nmax = 6  

function init(q1, q2, dx, jmax)
    xs = -1 # 始点
    x = range(xs, stop=xs+dx*(jmax-1), length=jmax)
    # initialize q[]
    q = Array{Float64}(undef, jmax) 
    for j = 1:jmax
      if x[j] < 0.0
        global q[j] = q1
      else
        global q[j] = q2
      end
    end
    return (q, x)
end
    
# grid generation
function do_computing_mc(x, q, c, dt, dx, nmax, interval)
  # plot of initial q
  plt = plot(x, q, marker=:circle, label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qbar = q[:]
    
    # MacCormack method
    for j = 1:jmax-1
      qbar[j] = q[j] - dt * c * (q[j+1] - q[j]) / dx
    end
    
    for j = 2:jmax-1
      q[j] = 0.5 * ((q[j] + qbar[j]) - dt * c * (qbar[j] - qbar[j-1]) / dx)
    end
    
    # 各ステップの可視化
    if n % 2 == 0
      global plt = plot!(x, q, marker=:circle, label="n=$n")
    end
  end 
end

interval = 2
q1 = 1.0
q2 = 0.0
q,x = init(q1, q2, dx, jmax)
do_computing_mc(x ,q, c, dt, dx, nmax, interval)

display(plt)

readline()

