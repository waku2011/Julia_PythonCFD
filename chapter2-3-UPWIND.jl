using Plots
gr()

const dt = 0.05 # time step, s
const dx = 0.1  # grid size, m

const jmax = 21 
const nmax = 6  

# grid generation
x = range(0, stop=dx*(jmax-1), length=jmax)

function do_computing(x ,q, c, dt, nmax, dx)
  # plot of initial q
  plt = plot(x, q, marker=:circle, label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  for n = 1:nmax
    qold = q
    for j = 2:jmax-1
      c2 = (c + abs(c)) / 2
      c3 = (c - abs(c)) / 2
      global q[j] = qold[j] - dt * ( c2 * (qold[j] - qold[j-1]) / dx 
                                   + c3 * (qold[j+1] - qold[j]) / dx) #式(2.23)
    end
    # 各ステップの可視化
    if n % 2 == 0
      global plt = plot!(x, q, marker=:circle, label="n=$n")
    end
  end 
end

c = 1    # positive advection velocity, m/s

# initialize q[]
q = Array{Float64}(undef, jmax)
for j = 1:jmax
  if j < jmax/2
    global q[j] = 1.0
  else
    global q[j] = 0.0
  end
end

do_computing(x ,q, c, dt, nmax, dx)

display(plt)

readline()

c = -1    # nagetive advection velocity, m/s

# initialize q[]
q = Array{Float64}(undef, jmax)
for j = 1:jmax
  if j < jmax/2+1
    global q[j] = 0.0
  else
    global q[j] = 1.0
  end
end

do_computing(x ,q, c, dt, nmax, dx)

display(plt)

readline()
