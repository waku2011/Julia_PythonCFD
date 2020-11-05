using Plots; gr()

const a = 1.0   # diffusion coeff.

const dt = 0.01 # time step
const dx = 0.1  # grid size

const jmax = 11
const nmax = 12 

function init(jmax)
    x = range(0, stop=dx*(jmax-1), length=jmax)
    q = sin.(pi * x)
    return (x, q)
end
    
# grid generation
function do_computing(x, q, a, dt, dx, nmax, interval)
  # plot of initial q
  plt = plot(x, q, marker=:circle, label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qold = copy(q)
    
    for j = 2:jmax-1
      dq = a * dt * (qold[j+1] - 2.0 * qold[j] + qold[j-1]) / (dx^2) # 式(4.6)
      q[j] = qold[j] + dq
    end
    q[1] = 0.0
    q[jmax] = 0.0
    
    # 各ステップの可視化
    if n % 2 == 0
      global plt = plot!(x, q, marker=:circle, label="n=$n")
    end
    
  end 
end

interval = 4
x,q = init(jmax)
do_computing(x ,q, a, dt, dx, nmax, interval)

display(plt)

readline()

