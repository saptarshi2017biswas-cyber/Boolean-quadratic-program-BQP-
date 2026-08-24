Consider the following BQP
max
x
x
T Cx
s.t. xi ∈ {−1, +1}, i = 1, . . . , n
(3)
Here C ∈ Sn with n = 100. Generate a symmetrix matrix C by first generating a
random A ∈ R
n×n and then performing C = AAT
. Solve the above problem by casting
it as SDP, and then relaxing it. You could use either Rank-1 approximation or Gaussian
randmomization to generate boolean x solution.
