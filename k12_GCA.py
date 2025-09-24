
archivo = "GCA_000001405.15_GRCh38_full_analysis_set.fna"
kmer = 12
 
# Tabla de complemento para A C T G
complemento = str.maketrans("ACGT", "TGCA")

# Ver que solo haya A C T G,
def solo_bases(string):
    for caracter in string:
        if caracter not in "ACGT":
            return False
    return True

def procesar(sequencia, k, conteo):
    n = len(sequencia)
    if n < k:
        return
    for i in range(n - k + 1):
        # Para poder sacar las líneas de tamaño K
        kmer = sequencia[i:i+k]
        # Solo consideramos las bases A C T G, si no, no >:(
        if not solo_bases(kmer):
            continue
        # Revertirlo por el tema de los complementos
        reverso = kmer.translate(complemento)[::-1]
        canon = min(kmer, reverso)
        conteo[canon] = conteo.get(canon, 0) + 1

# Ya, entonces ahora se hace el conteo y junte.
def contar_kmer(link, k):
    dicc = {}
    seq_linea = []
    with open(link, "rt") as archivo:
        for linea in archivo:
            if linea and linea[0] == ">":
                if seq_linea:
                    seq = ""
                    for l in seq_linea:
                        seq += l
                    seq = seq.upper().replace("\r", "")
                    procesar(seq, k, dicc)
                    seq_linea = []
            else:
                seq_linea.append(linea.strip())

        if seq_linea:
            seq = ""
            for l in seq_linea:
                seq += l
            seq = seq.upper().replace("\r", "")
            procesar(seq, k, dicc)
    return dicc

# Para hacer el histograma
def histograma(dicc, tope = 100):
    hist = {}
    for valor in dicc.values():
        canon = min(valor, tope)
        hist[canon] = hist.get(canon, 0) + 1
    return hist

conteo = contar_kmer(archivo, kmer)
hist = histograma(conteo, 100)

items = sorted(hist.items(), key=lambda kv: (-kv[1], kv[0]))

# Printear
for i, j in items:
    print(f"{i} {j}")
items = sorted(hist.items(), key=lambda kv: (-kv[1], kv[0]))
