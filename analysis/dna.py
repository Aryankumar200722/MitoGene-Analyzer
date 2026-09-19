import pandas as pd

def nucleotide_stats(seq):
    a,t,g,c = seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")
    n = len(seq)
    return {
        "A": a, "T": t, "G": g, "C": c,
        "GC%": ((g+c)/n)*100 if n else 0,
        "AT%": ((a+t)/n)*100 if n else 0,
        "GC skew": ((g-c)/(g+c)) if (g+c) else 0,
        "AT skew": ((a-t)/(a+t)) if (a+t) else 0,
    }

def sliding_gc(seq, window=100, step=25):
    rows = []
    for start in range(0, max(1, len(seq)-window+1), step):
        chunk = seq[start:start+window]
        if len(chunk) < window:
            break
        gc = (chunk.count("G")+chunk.count("C"))/len(chunk)*100
        rows.append({"Position": start+1, "GC%": gc})
    return pd.DataFrame(rows)
