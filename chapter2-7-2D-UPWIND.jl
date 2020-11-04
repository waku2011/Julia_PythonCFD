using Plots
gr()

c = 1
d = 1

dt = 0.04
dx = 0.1  
dy = 0.1  

jmax = 41 
kmax = 41 
nmax = 25 

function init(dx, dy, jmax, kmax)
  x = range(0, dx*(jmax-1), length=jmax)
  y = range(0, dy*(kmax-1), length=kmax)
  q = init_q(x, y, jmax, kmax)
  return (x, y, q)
end

function init_q(x, y, jmax, kmax)  # simple Gaussian-hill
  q = Array{Float64}(undef, jmax, kmax)
  mu    = 1.0
  sigma = 0.2
  for k = 1:kmax
    for j = 1:jmax
       q[j,k] = exp( -(x[j]-mu)^2/(2.0*sigma^2)-(y[k]-mu)^2/(2.0*sigma^2) )
    end
  end
  return (q)
end

function do_computing(x, y, q, c, d, dt, dx, dy, nmax, interval::Int64=2)
  # plot of initial q
  plt = contourf(x, y, q', 
                 title="2D Gauss-hill (n=0)",
                 xlabel="x", xlims=(0,4), 
                 ylabel="y", ylims=(0,4),
                 framestyle = :box, 
                 color=:jet, alpha=1.0, linewidth=0, levels=19, clim=(0,1),
                 aspect_ratio=1, size=(800,800))
  display(plt)
  println("Press enter to continue ...")
  readline()
  
  # time stepping
  for n = 1:nmax
    qold = copy(q)
    println(maximum(q)," ", minimum(q))
    for k = 2:kmax-1
       for j = 2:jmax-1
         q[j,k] = qold[j,k]-dt/dx*c*(qold[j,k]-qold[j-1,k])-dt/dy*d*(qold[j,k]-qold[j,k-1])
       end
    end
    
    # 各ステップの可視化
    if n % interval == 0
      plt = contourf(x, y, q', title="2D Gauss-hill (n=$n)", 
                 xlabel="x", xlims=(0,4), 
                 ylabel="y", ylims=(0,4),
                 framestyle = :box, 
                 color=:jet, alpha=1.0, linewidth=0, levels=19, clim=(0,1),
                 aspect_ratio=1, size=(800,800))
      display(plt)
      println("Press enter to continue ...")
      readline()
    end
  end 
end

x, y, q = init(dx, dy, jmax, kmax)
interval = 2
do_computing(x, y, q, c, d, dt, dx, dy, nmax, interval)

