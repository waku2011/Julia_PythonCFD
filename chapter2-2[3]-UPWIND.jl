using Plots
gr()

const c  = 1    # advection velocity, m/s
const dt = 0.05 # time step, s
const dx = 0.1  # grid size, m

const jmax = 21 
const nmax = 6  

# grid generation
x = range(0, stop=dx*(jmax-1), length=jmax)

# initialize q[]
q = Array{Float64}(undef, jmax)
for j = 1:jmax
  if j < jmax/2+1
    global q[j] = 1.0
  else
    global q[j] = 0.0
  end
end
plt = plot(x, q, marker=:circle, label="n=0",xlabel="x", ylabel="q")

# time stepping
for n = 1:nmax
  qold = copy(q)
  for j = 2:jmax  # Not 2:jmax-1
     global q[j] = qold[j] - dt * c *(qold[j] - qold[j-1])/dx  # 式(2.9)
  end
  # 各ステップの可視化
  if n % 2 == 0
    global plt = plot!(x, q, marker=:circle, label="n=$n")
  end
end

display(plt)

readline()
