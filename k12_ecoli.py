
archivo = "ecoli-k12-ref.fna"
kmer = 12
 
# Tabla de complemento para A,C,G,T
comp_tab = str.maketrans("ACGT", "TGCA")

def only(s):
    for ch in s:
        if ch not in "ACGT":
            return False
    return True

def procesar(seq, k, counts):
    n = len(seq)
    if n < k:
        return
    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if not only(kmer):
            continue
        rev = kmer.translate(comp_tab)[::-1]
        canon = kmer if kmer < rev else rev
        counts[canon] = counts.get(canon, 0) + 1

def count_from_fasta_gz(path, k):
    counts = {}
    seq_buf = []
    with open(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if seq_buf:
                    seq = "".join(seq_buf).upper().replace("\r", "")
                    procesar(seq, k, counts)
                    seq_buf = []
            else:
                seq_buf.append(line.strip())

        if seq_buf:
            seq = "".join(seq_buf).upper().replace("\r", "")
            procesar(seq, k, counts)
    return counts

def build_hist(counts, cap):
    hist = {}
    for cnt in counts.values():
        b = cnt if cnt <= cap else cap
        hist[b] = hist.get(b, 0) + 1
    return hist

counts = count_from_fasta_gz(archivo, kmer)
hist = build_hist(counts, 100)

items = sorted(hist.items(), key=lambda kv: (-kv[1], kv[0]))

for a, b in items[:10]:
    print(f"{a}\t{b}")
