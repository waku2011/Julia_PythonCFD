using Plots; gr()

const c = 1.0   # advection speed

const dt = 0.05 # time step
const dx = 0.1  # grid size

const jmax = 70
const nmax = 50 

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
    return (x, q)
end

function UPWIND1(alf, q, c, dt, dx, jmax)
    for j = 1:jmax-1
        ur, ul = q[j+1], q[j]
        fr, fl = c * ur, c * ul
        alf[j] = 0.5 * (fr + fl - abs(c) * (ur - ul)) # 式(2.34)
    end
end

function myDiff(alf, jmax)
   tmp = zeros(jmax-1)
   for i=1:jmax-1
   	 tmp[i]=alf[i+1]-alf[i]
   end
   return(tmp)
end

# grid generation
function do_computing_LDU(x, q, c, dt, dx, nmax, ff, interval)
  # plot of initial q
  plt = plot(x, q, marker=:circle, label="n=0",xlabel="x", ylabel="q")
  
  # time stepping
  alf = zeros(jmax)
  dq = zeros(jmax)
  R = zeros(0)
    
  for n = 1:nmax
    qold = copy(q)
        
    # 近似LDU分解
    c_a = abs(c)
    c_p = 0.5 * (c + c_a)
    c_n = 0.5 * (c - c_a)
    nu_a = c_a * dt / dx
    nu_p = c_p * dt / dx
    nu_n = c_n * dt / dx
        
    ff(alf, qold, c, dt, dx, jmax)
    R = vcat(0.0, myDiff(alf,jmax)/dx)

    ## 第一スイープ
    for j = 2:jmax-1
       dq[j] = (-dt * R[j] + nu_p * dq[j - 1]) / (1 + nu_a)
    end
    
    ## 第二、第三スイープ
    for j = jmax-1:-1:1
       dq[j] = dq[j] - nu_n * dq[j + 1] / (1 + nu_a)
    end        
    
    for j = 2:jmax-1
       q[j] = qold[j] + dq[j]
    end        
    
    # 各ステップの可視化
    if n % interval == 0
      global plt = plot!(x, q, marker=:circle, label="n=$n")
    end
  end 
end

interval = 16
q1 = 1.0
q2 = 0.0
x, q = init(q1, q2, dx, jmax)
do_computing_LDU(x ,q, c, dt, dx, nmax, UPWIND1, interval)

display(plt)

readline()

