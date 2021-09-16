from MSIS_security import MSIS_summarize_attacks, MSISParameterSet

# n = 128
# w = 2 ** 11
# h = 2
# B = 16
# q = 2 ** 13
# norm="linf"

# n = 128
# w = 36864
# h = 1
# B = 1
# q = 2 ** 6
# norm="linf"

# n = 128
# w = 1536   
# h = 1
# B = 1
# q = 2 ** 6
# norm="linf"

# n = 128
# w = 12
# h = 1
# B = 1
# q = 2**6
# norm="linf"

# n = 128
# w = 192
# h = 3
# B = 1032
# q = 2**32
# norm="linf"

# n = 128
# w = 96
# h = 8
# B = 2**20
# q = 2**32
# norm="linf"

# n = 128
# w = 48
# h = 15
# B = 2**30
# q = 2**32
# norm="linf"

#PARAMS II
# n = 256
# w = 9
# h = 5
# B = 322560
# q = 3870721
# norm="linf"

#PARAMS I
# n = 256
# w = 7   
# h = 4
# B = 168448
# q = 2021377
# norm="linf"


#PARAMS III
# n = 256
# w = 11
# h = 6
# B = 322560
# q = 3870721
# norm="linf"

# n = 64
# w = 12
# h = 2
# B = 1
# q = 2**6
# norm="linf"

# n = 64
# w = 192
# h = 5
# B = 2**9
# q = 2**32
# norm="linf"

# n = 64
# w = 96
# h = 13
# B = 2**18
# q = 2**32
# norm="linf"

# n = 64
# w = 48
# h = 25
# B = 2**27
# q = 2**32
# norm="linf"

# n = 64
# w = 128
# h = 4
# B = 625602
# q = 2**32
# norm="linf"

# n = 64
# w = 24
# h = 2
# B = 1
# q = 2**6
# norm="linf"


# n = 64
# w = 15
# h = 9
# B = 8475971
# q = 2**32
# norm="l2"

# n = 64
# w = 12
# h = 1
# B = 19954835
# q = 2**32
# norm="l2"

# n = 64
# w = 160
# h = 10
# B = 57465
# q = 2**32
# norm="linf"

# n = 64
# w = 24
# h = 2
# B = 1
# q = 2**6
# norm="linf"

# n = 64
# w = 384
# h = 30
# B = 2**30
# q = 2**32
# norm="linf"

# n = 64
# w = 192
# h = 33
# B = 2**36
# q = 2**40
# norm="linf"

# n = 64
# w = 96
# h = 36
# B = 2**42
# q = 2**48
# norm="linf"

# n = 64
# w = 192
# h = 31
# B = 2**36
# q = 2**42
# norm="linf"

# n = 64
# w = 192
# h = 27
# B = 2**36
# q = 2**48
# norm="linf"

# n = 64
# w = 96
# h = 35
# B = 2**42
# q = 2**50
# norm="linf"


# n = 64
# w = 16
# h = 1
# B = 1
# q = 2**8
# norm="linf"

# n = 64
# w = 20
# h = 1
# B = 1
# q = 2**10
# norm="linf"

# n = 64
# w = 384
# h = 29
# B = 2**37
# q = 2**48
# norm="linf"

# n = 64
# w = 192
# h = 53
# B = 2**59
# q = 2**60
# norm="linf"

# n = 64
# w = 96
# h = 100
# B = 2**75
# q = 2**80
# norm="linf"

# n = 256
# w = 6
# h = 1
# B = 1
# q = 2**3
# norm="linf"

# n = 64
# w = 384
# h = 20
# B = 2**25
# q = 2**32
# norm="linf"

# n = 64
# w = 384
# h = 39
# B = 2**44
# q = 2**48
# norm="linf"

# n = 64
# w = 192
# h = 43
# B = 2**50
# q = 2**55
# norm="linf"

# n = 64
# w = 96
# h = 50
# B = 2**56
# q = 2**58
# norm="linf"

# n = 64
# w = 96
# h = 35
# B = 2**43
# q = 2**50
# norm="linf"

n = 128
w = 12
h = 1
B = 1
q = 2**6
norm="linf"

text_SIS = ["BKZ block-size $b$ to break SIS","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
table_SIS = [1*[0] for i in range(4)]


ParameterSet = MSISParameterSet(n, w, h, B, q, norm)

v = MSIS_summarize_attacks(ParameterSet)

for i in range(4):
    table_SIS[i] = v[i]


for i in range(4):
    print(text_SIS[i])
    print("========================")
    print(table_SIS[i])


table_SIS = [1*[0] for i in range(4)]