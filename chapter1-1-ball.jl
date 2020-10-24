using Plots
gr()

const g = 9.8   # 重力加速度
const v0 = 10   # 初速度
const h0 = 0    # 初期高さ
const dt = 0.1  # 時間刻み

# 解析解(Analytical Solution)
ta = range(0, stop=2*v0/g, length=100)
ha = - 0.5*g*(ta.^2)  .+ v0*ta .+ h0  # 式(1.1)
plt = plot(ta, ha, linecolor=:blue, label="Analytical",xlabel="time,s", ylabel="height,m")

# 数値解1(Numerical Solution)
t = Float64[]
h = Float64[]
t1 = 0
h1 = h0
push!(t,t1)
push!(h,h1)
while h1 >= 0
  global h1 +=  (-g * t1 + v0) * dt # 式(1.7)
  global t1 += dt 
  push!(t,t1)
  push!(h,h1)
end
plt = scatter!(t, h, label="Numarical 1")

# 数値解2(Numerical Solution)
t = Float64[]
h = Float64[]
t2 = 0
h2 = h0
push!(t,t2)
push!(h,h2)
while h2 >= 0
  global h2 +=  (-g * (t2+(t2+dt))/2 + v0) * dt # 式(1.7)
  global t2 += dt 
  push!(t,t2)
  push!(h,h2)
end
plt = scatter!(t, h, label="Numarical 2")

display(plt)

readline()
