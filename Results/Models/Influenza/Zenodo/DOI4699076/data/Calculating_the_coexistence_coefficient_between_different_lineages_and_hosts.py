"""
Calculating the coexistence coefficient between different lineages and hosts
"""

import pandas as pd

# IO
h1_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\H1_lineage.csv",header=0)
h2_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\H2_lineage.csv",header=0)
s1_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\S1_lineage.csv",header=0)
s2_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\S2_lineage.csv",header=0)
s3_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\S3_lineage.csv",header=0)
a1_f = pd.read_csv(r"G:\BACK_UP_RT\Infleunza2019_host\data\H1N1\clusters\data\A1_lineage.csv",header=0)

#
h1_region = h1_f["region"]
h1_time = h1_f["time"]
h2_region = h2_f["region"]
h2_time = h2_f["time"]
s1_region = s1_f["region"]
s1_time = s1_f["time"]
s2_region = s2_f["region"]
s2_time = s2_f["time"]
s3_region = s3_f["region"]
s3_time = s3_f["time"]
a_region = a1_f["region"]
a_time = a1_f["time"]


# gen total site
def combine(sp1,sp2):
    cb = []
    for i in range(len(sp1)):
        cb = cb + ["%s_%s"%(round(sp1[i]),sp2[i])]
    return cb


h1 = combine(h1_time, h1_region)
h2 = combine(h2_time, h2_region)
s1 = combine(s1_time, s1_region)
s2 = combine(s2_time, s2_region)
s3 = combine(s3_time, s3_region)
a = combine(a_time, a_region)
h = h1+h2
s = s1 + s2 + s3

h1_h2 = list(set(h1+h2))
s1_s2 = list(set(s1+s2))
s1_s3 = list(set(s1+s3))
s2_s3 = list(set(s2+s3))

s_a = list(set(s+a))
h_s = list(set(s+h))
h_a = list(set(h+a))

# CI
h1_h2_ci = 0
s1_s2_ci = 0
s1_s3_ci = 0
s2_s3_ci = 0
s_a_ci = 0
h_s_ci = 0
h_a_ci = 0

# by lineages
for site in h1_h2:
    h1_h2_ci = h1_h2_ci + abs(h1.count(site) / len(h1) - h2.count(site)/ len(h2))

for site in s1_s2:
    s1_s2_ci = s1_s2_ci + abs(s1.count(site) / len(s1) - s2.count(site)/ len(s2))

for site in s1_s3:
    s1_s3_ci = s1_s3_ci + abs(s1.count(site) / len(s1) - s3.count(site)/ len(s3))

for site in s2_s3:
    s2_s3_ci = s2_s3_ci + abs(s2.count(site) / len(s2) - s3.count(site)/ len(s3))

# by hosts
for site in h_s:
    h_s_ci = h_s_ci + abs(h.count(site) / len(h) - s.count(site)/ len(s))

for site in h_a:
    h_a_ci = h_a_ci + abs(h.count(site) / len(h) - a.count(site)/ len(a))

for site in s_a:
    s_a_ci = s_a_ci + abs(s.count(site) / len(s) - a.count(site)/ len(a))

h1_h2_ci = 1 - 0.5*h1_h2_ci  # H1 V.S. H2
s1_s2_ci = 1 - 0.5*s1_s2_ci  # S1 V.S. S2
s1_s3_ci = 1 - 0.5*s1_s3_ci  # S1 V.S. S3
s2_s3_ci = 1 - 0.5*s2_s3_ci  # S2 V.S. S3
h_s_ci = 1 - 0.5*h_s_ci  # Human V.S. Swine
h_a_ci = 1 - 0.5*h_a_ci  # Human V.S. Avian
s_a_ci = 1 - 0.5*s_a_ci  # Swine V.S. Avian

print(h1_h2_ci,s1_s2_ci,s1_s3_ci,s2_s3_ci,h_s_ci,h_a_ci,s_a_ci)