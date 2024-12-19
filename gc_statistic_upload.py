import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

# Función para calcular el contenido GC
def gc_content(sequence):
    # Limpiar la secuencia y pasarla a mayúsculas
    sequence = sequence.rstrip("\n").upper()
    # Contar G y C, calcular el porcentaje
    gc_count = sequence.count("G") + sequence.count("C")
    return round((gc_count / len(sequence)) * 100, 2)

# Leer secuencia del archivo FASTA
def read_fasta(filename, max_lines=100):
    sequence = ""
    with open(filename) as file:
        next(file)  # Saltar la primera línea (cabecera del FASTA)
        for i, line in enumerate(file):
            if i >= max_lines:
                break
            sequence += line.strip()
    return sequence

# Configuración inicial
filename = "B_anthracis.fasta"
sequence = read_fasta(filename)

# Calcular el contenido GC de la secuencia original
original_gc = gc_content(sequence)
print(f"GC Content (original): {original_gc}%")

# Generar secuencias aleatorias y calcular su contenido GC
def generate_random_gc(sequence, replicates=100, subsequence_length=10):
    random_sequence = list(sequence)
    random_gc_list = []

    for _ in range(replicates):
        random.shuffle(random_sequence)
        subsequence = random_sequence[:subsequence_length]
        subsequence = "".join(subsequence)
        random_gc_list.append(gc_content(subsequence))

    return pd.Series(random_gc_list)

# Parámetros ajustables
replicates = 1000
subsequence_length = 10
random_gc_list = generate_random_gc(sequence, replicates, subsequence_length)

# Calcular p-value
statistic = 50  # Valor estadístico para comparar
p_value = sum(random_gc_list >= statistic) / replicates
print(f"P-value: {p_value}")

# Gráficos
plt.figure(figsize=(10, 6))
sns.histplot(random_gc_list, kde=True, color="skyblue", bins=20, label="GC Content (Random)")
plt.axvline(statistic, color="red", linestyle="--", label=f"Statistic ({statistic}%)")
plt.axvline(original_gc, color="green", linestyle="--", label=f"Original GC Content ({original_gc}%)")
plt.title("Distribution of GC Content in Randomized Subsequences")
plt.xlabel("GC Content (%)")
plt.ylabel("Frequency")
plt.legend()
plt.tight_layout()

# Mostrar el gráfico
plt.show()