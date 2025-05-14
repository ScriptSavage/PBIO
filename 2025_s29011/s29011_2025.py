"""
program tworzy losowe sekwencje DNA, zapisuje je w formacie FASTA oraz
wylicza podstawowe statystyki nukleotydowe, które zapisuje do pliku CSV oraz FASTA

Narzedzie moze sluzyc do :
1) prostych analiz w bioinformatyce
2) szybkiego generowania przykladowych danych
"""

import random                       # standardowy PRNG Pythona
import textwrap                     # zawijanie długich linii sekwencji w FASTA
import sys                          # wyjście awaryjne z kodem błędu
from pathlib import Path            # wygodna obsługa ścieżek plików
from collections import Counter     # szybkie zliczanie elementów w sekwencji
from dataclasses import dataclass    # deklaratywne klasy danych (Stats)
import csv                          # zapis wyników w formacie CSV
from typing import Optional         # aliasy typów dla parametrów opcjonalnych
import re                           # walidacja identyfikatora sekwencji


DNA_ALPHABET = "ACGT"              # alfabet nukleotydów dla DNA
MyName = "Maksymilian"
# ORIGINAL:
# def generate_random_dna(length: int) -> str:
#     return "".join(random.choices(DNA_ALPHABET, k=length))
# MODIFIED (dodano *rng* – łatwiejsza testowalność kodu dla różnych przypadków)
#Generuje losową sekwencję DNA o zadanej długości
def generate_random_dna(length: int, *, rng: Optional[random.Random] = None) -> str:
    r = rng or random
    return "".join(r.choices(DNA_ALPHABET, k=length))


def insert_signature(sequence: str, signature: str, *, rng: Optional[random.Random] = None) -> str:
    """Zwraca sekwencję z wstawionym *signature* w losowym miejscu."""
    r = rng or random
    idx = r.randint(0, len(sequence))
    return sequence[:idx] + signature + sequence[idx:]

# Dataclass przechowywanie statystyk sekwencji
@dataclass
class Stats:
    length: int                     # całkowita długość sekwencji (nt)
    counts: dict[str, int]          # ile razy występuje każdy nukleotyd
    percentages: dict[str, float]   # odsetek każdego nt w %
    cg_at_ratio: Optional[float]    # (C+G)/(A+T) lub None jeśli A+T==0


# ORIGINAL:
# def calc_statistics(sequence: str) -> dict:
#     length = len(sequence)
#     counts = {nuc: sequence.count(nuc) for nuc in DNA_ALPHABET}
#     percentages = {nuc: round(counts[nuc] / length * 100, 2) for nuc in DNA_ALPHABET}
#     at_sum = counts["A"] + counts["T"]
#     cg_sum = counts["C"] + counts["G"]
#     cg_at_ratio = round(cg_sum / at_sum, 3) if at_sum else None
#     return {"length": length, "counts": counts, "percentages": percentages, "cg_at_ratio": cg_at_ratio}
# MODIFIED (Counter + Stats: jedno przejście po sekwencji, szybszy i czytelniejszy wynik):
def calc_statistics(sequence: str) -> Stats:
    """Zwraca statystyki bazowe dla podanej sekwencji DNA."""
    length = len(sequence)                              # długość sekwencji
    counts = Counter(sequence)                          # liczenie wszystkich znaków
    # tworzymy słownik z pelnym alfabetem – Counter zwraca 0 gdy klucza brak
    counts_full = {n: counts[n] for n in DNA_ALPHABET}
    percentages = {n: round(counts_full[n] / length * 100, 2) for n in DNA_ALPHABET}
    at_sum = counts_full['A'] + counts_full['T']        # suma A+T
    cg_sum = counts_full['C'] + counts_full['G']        # suma C+G
    ratio = round(cg_sum / at_sum, 3) if at_sum else None
    return Stats(length, counts_full, percentages, ratio)



# === ZAPIS FASTA ===
# ORIGINAL:
# def save_fasta(filename: Path, header: str, sequence: str, width: int = 80) -> None:
#     with filename.open("w", encoding="utf-8") as f:
#         f.write(f">{header}\n")
#         for line in textwrap.wrap(sequence, width):
#             f.write(f"{line}\n")
# MODIFIED (ladniejszy docstring):
#zapis do pliku fasta
def save_fasta(filename: Path, header: str, sequence: str, width: int = 80) -> None:
    seq_lines = textwrap.wrap(sequence, width) if width else [sequence]
    with filename.open("w", encoding="utf-8") as f:
        f.write(f">{header}\n")                    # nagłówek FASTA
        for line in seq_lines:
            f.write(f"{line}\n")                  # kolejne fragmenty sekwencji

# ORIGINAL ( Brak zapisu do CSV )
# MODIFIED ( dodanie zapisu do CSV )
def save_stats_csv(stats: Stats, csv_path: Path, seq_id: str) -> None:
    file_exists = csv_path.exists()
    with csv_path.open("a", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:                            # nagłówek tylko raz
            writer.writerow([
                "ID", "Length", "A", "C", "G", "T",
                "%A", "%C", "%G", "%T", "(C+G)/(A+T)"
            ])
        # wiersz danych dla jednej sekwencji
        writer.writerow([
            seq_id,
            stats.length,
            stats.counts['A'], stats.counts['C'], stats.counts['G'], stats.counts['T'],
            stats.percentages['A'], stats.percentages['C'], stats.percentages['G'], stats.percentages['T'],
            stats.cg_at_ratio if stats.cg_at_ratio is not None else "NA",
        ])

def main() -> None:
    """Pobiera dane od użytkownika, generuje sekwencję, zapisuje FASTA i CSV."""
    try:
        length = int(input("Podaj długość sekwencji (liczba dodatnia): "))
        if length <= 0:
            raise ValueError
    except ValueError:
        sys.exit("BŁĄD: Długość musi być dodatnia!")

    # — ID sekwencji —
    # ORIGINAL:
    # seq_id = input("Podaj ID sekwencji: ").strip()
    # MODIFIED (walidacja ID):
    seq_id = input("Podaj ID sekwencji: ").strip()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", seq_id):
        sys.exit("ERROR: ID może zawierać tylko litery, cyfry, podkreślenia, kropki i myślniki.")

    description = input("Podaj opis sekwencji: ").strip()

    # -- generowanie i statystyki --------------------------------------------
    dna_seq = generate_random_dna(length)              # losowy łańcuch nukleotydów
    stats = calc_statistics(dna_seq)                   # obliczenia Counter‐based

    # wstawiamy podpis do sekwencji zapisywanej w FASTA
    final_seq = insert_signature(dna_seq, MyName)

    # -- zapis plików ---------------------------------------------------------
    fasta_path = Path(f"{seq_id}.fasta")               # nazwa pliku FASTA
    save_fasta(fasta_path, f"{seq_id} {description}", final_seq)

    csv_path = Path(f"{seq_id}_stats.csv")            # nazwa pliku CSV
    save_stats_csv(stats, csv_path, seq_id)            # dopisanie wiersza

    # -- finalny rezultat
    print("\n⚙  Plik FASTA zapisano jako:", fasta_path)
    print("Plik CSV  zapisano jako:", csv_path)
    print("\nStatystyki sekwencji:")
    for nuc in DNA_ALPHABET:
        print(f"  {nuc}: {stats.percentages[nuc]}% ({stats.counts[nuc]} nt)")
    if stats.cg_at_ratio is not None:
        print(f"  Stosunek (C+G)/(A+T): {stats.cg_at_ratio}")
    else:
        print("  Nie można obliczyć stosunku – brak A i T w sekwencji.")

if __name__ == "__main__":
    main()
