using Plots
gr()

c = 1.0
d = 1.0

dt = 0.05 
dx = 0.1  
dy = 0.1  

jmax = 41 
kmax = 41 
nmax = 30 

function init(dx, dy, jmax, kmax)
  x = range(0, dx*(jmax-1), length=jmax)
  y = range(0, dy*(kmax-1), length=kmax)
  q = init_q(x, y, jmax, kmax)
  return (x, y, q)
end

function init_q(x, y, jmax, kmax)  # simple Gaussian-hill
  q = Array{Float64}(undef, jmax, kmax)
  sigma, mu = 0.1, 1.0
  for j = 1:kmax
    for i = 1:jmax
       q[i,j] = exp( -(x[i]-mu)^2/(2.0*sigma^2) -(y[j]-mu)^2/(2.0*sigma^2) )
    end
  end
  return (q)
end

function do_computing(x, y, q, c, d, dt, dx, dy, nmax, interval::Int64=2)
  # plot of initial q
  plt = plot(x, y, q, title="2D Gauss-hill", label="n=0",xlabel="x", ylabel="y")
  
  # time stepping
  for n = 1:nmax
    qold = q
    for k = 2:kmax-1
       for j = 2:jmax-1
         q[j,k] = qold[j,k] - dt / dx *  c * (qold[j,k] - qold[j-1,k]) - dt / dy *  d * (qold[j,k] - qold[j,k-1])
       end
    end
    
    # 各ステップの可視化
    if n % interval == 0
      plt = plot!(x, y, q, title="2D Gauss-hill", label="n=$n", xlabel="x", ylabel="y")
      display(plt)
    end
    
    readline()   
  end
end

x, y, q = init(dx, dy, jmax, kmax)
interval = 2
do_computing(x, y, q, c, d, dt, dx, dy, nmax, interval)

