from sys import stdin



prev_s, prev_e = 0,0
prev_entropy, prev_overhang = 0, 0
for l in stdin:
    contig, s, e, strand, entropy, min_overhang = l.rstrip().split()
    entropy = float(entropy)
    min_overhang = int(min_overhang)
    s,e  = int(s), int(e)
    if min_overhang >= 6 and entropy >=2:
        #print contig, s, e, strand, entropy, min_overhang
        if (prev_s == s and strand == "+") or (prev_e == e and strand=="-"):
            if strand == "+":
                d = abs(prev_e-e)
                print contig, prev_s, prev_e, s, e, strand, d, entropy, prev_entropy
            else:
                d = abs(prev_s-s)
                print contig, s, e, prev_s, prev_e, strand, d, prev_entropy, entropy

        prev_s, prev_e = s,e 
        prev_entropy = entropy
        prev_overhang = min_overhang
