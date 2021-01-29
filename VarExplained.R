
K = 0.01
nc = 897
nt = 1629
P = nc/(nc+nt)
r2 = 0.09

t = qnorm(1-K)
z = dnorm(t)
m = z/K

C = K^2*(1-K)^2 / (z^2*P*(1-P))
theta = m * (P-K)/(1-K) * (m*(P-K)/(1-K)-t)
r2_liab = r2*C / (1+r2*theta*C)

cat(sprintf("R^2 on the liability scale: %g\n",r2_liab))