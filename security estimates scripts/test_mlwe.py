from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet

# n = 64
# l = 187
# k = 5
# eta = 2 ** 9
# q = 2 ** 32
# distr="uniform"

# n = 64
# l = 83
# k = 13
# eta = 2 ** 18
# q = 2 ** 32
# distr="uniform"

# n = 64
# l = 23
# k = 25
# eta = 2 ** 27
# q = 2 ** 32
# distr="uniform"

# n = 256
# l = 4
# k = 4
# eta = 2
# q = 8380417
# distr="uniform"

# n = 64
# l = 11
# k = 1
# eta = 1
# q = 2 ** 32
# distr="uniform"

# n = 64
# l = 11
# k = 1
# eta = 1
# q = 2 ** 32
# distr="uniform"

# n = 64
# l = 18
# k = 6
# eta = 1
# q = 2 ** 32
# distr="uniform"

n = 128
l = 12
k = 6
eta = 1
q = 2 ** 32
distr="uniform"

text_LWE = ["BKZ block-size $b$ to break LWE","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
table_LWE = [1*[0] for i in range(4)]


# ParameterSet = MSISParameterSet(n, w, h, B, q, norm)

ParameterSet = MLWEParameterSet(n, l, k, eta, q, distr)

v = MLWE_summarize_attacks(ParameterSet)

for i in range(4):
    table_LWE[i] = v[i]


for i in range(4):
    print(text_LWE[i])
    print("========================")
    print(table_LWE[i])


table_LWE = [1*[0] for i in range(4)]