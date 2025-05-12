import random
import textwrap
import sys
from pathlib import Path

DNA_ALPHABET = "ACGT"
SIGNATURE = "Max"

def generate_random_dna(length: int) -> str:

    return "".join(random.choices(DNA_ALPHABET, k=length))

def insert_signature(sequence: str, signature: str) -> str:
    idx = random.randint(0, len(sequence))
    return sequence[:idx] + signature + sequence[idx:]

def calc_statistics(sequence: str) -> dict:
    length = len(sequence)
    counts = {nuc: sequence.count(nuc) for nuc in DNA_ALPHABET}
    percentages = {nuc: round(counts[nuc] / length * 100, 2) for nuc in DNA_ALPHABET}
    at_sum = counts["A"] + counts["T"]
    cg_sum = counts["C"] + counts["G"]
    cg_at_ratio = round(cg_sum / at_sum, 3) if at_sum else None
    return {
        "length": length,
        "counts": counts,
        "percentages": percentages,
        "cg_at_ratio": cg_at_ratio,
    }

def save_fasta(filename: Path, header: str, sequence: str, width: int = 80) -> None:
    with filename.open("w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for line in textwrap.wrap(sequence, width):
            f.write(f"{line}\n")

def main() -> None:
    try:
        length = int(input("Podaj długość sekwencji (liczba dodatnia): "))
        if length <= 0:
            raise ValueError("Długość musi być dodatnia!")
    except ValueError as err:
        print("Błąd:", err)
        sys.exit(1)

    seq_id = input("Podaj ID sekwencji: ").strip()
    description = input("Podaj opis sekwencji: ").strip()

    dna_seq = generate_random_dna(length)
    stats = calc_statistics(dna_seq)

    final_seq = insert_signature(dna_seq, SIGNATURE)

    fasta_file = Path(f"{seq_id}.fasta")
    save_fasta(fasta_file, f"{seq_id} {description}", final_seq)

    print("\nPlik FASTA zapisano jako:", fasta_file)
    print("Statystyki sekwencji (bez podpisu):")
    for nuc in DNA_ALPHABET:
        print(f"  {nuc}: {stats['percentages'][nuc]}% ({stats['counts'][nuc]} nt)")
    if stats['cg_at_ratio'] is not None:
        print(f"Stosunek (C+G)/(A+T): {stats['cg_at_ratio']}")
    else:
        print("Nie można obliczyć stosunku – brak A i T w sekwencji.")

if __name__ == "__main__":
    main()
