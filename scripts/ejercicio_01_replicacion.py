"""
Ejercicio 1: Replicación del ADN y Extensión con Biopython.

Este script define la secuencia de ADN dada, calcula su cadena complementaria
y su reversa complementaria empleando Biopython, y verifica los resultados
frente a las resoluciones teóricas manuales.
"""

from Bio.Seq import Seq


def main() -> None:
    # 1. Secuencia molde / parental dada (5' -> 3')
    seq_5_to_3_str = "ATGCCGTTAGCT"
    dna_seq = Seq(seq_5_to_3_str)

    print("=" * 60)
    print("EJERCICIO 1: REPLICACIÓN DEL ADN - EXTENSIÓN CON BIOPYTHON")
    print("=" * 60)
    print(f"Hebra parental original (5' -> 3'): {dna_seq}")

    # 2. Generación de la hebra complementaria
    # En Biopython, Seq.complement() sustituye A<->T y C<->G posición a posición,
    # por lo que el resultado mantiene la orientación relativa 3' -> 5'.
    dna_complement = dna_seq.complement()
    print(f"Hebra complementaria (3' -> 5'):    {dna_complement}")

    # En biología molecular estandarizada (dirección 5' -> 3'):
    dna_rev_complement = dna_seq.reverse_complement()
    print(f"Hebra complementaria (5' -> 3'):    {dna_rev_complement}")
    print("-" * 60)

    # 3. Comprobación contra los resultados obtenidos manualmente
    manual_complement_3_5 = "TACGGCAATCGA"
    manual_rev_complement_5_3 = "AGCTAACGGCAT"

    assert str(dna_complement) == manual_complement_3_5, (
        f"Discrepancia en complementaria (3'->5'): esperado {manual_complement_3_5}, obtenido {dna_complement}"
    )
    assert str(dna_rev_complement) == manual_rev_complement_5_3, (
        f"Discrepancia en reversa complementaria (5'->3'): esperado {manual_rev_complement_5_3}, obtenido {dna_rev_complement}"
    )

    print(">> VERIFICACIÓN MANUAL vs BIOPYTHON: EXITOSA (COINCIDENCIA EXACTA)")
    print("-" * 60)

    # 4. Esquematización de la replicación semiconservativa (1 ronda)
    print("\nEsquema de Replicación Semiconservativa:")
    print("Dúplex parental:")
    print(f"  5'- {dna_seq} -3'")
    print(f"  3'- {dna_complement} -5'")
    print("\nTras una ronda de replicación se originan 2 moléculas dúplex hijas:")
    print("  Molécula Hija 1:")
    print(f"    [Parental]    5'- {dna_seq} -3'")
    print(f"    [Nueva hebra] 3'- {dna_complement} -5'")
    print("  Molécula Hija 2:")
    print(f"    [Nueva hebra] 5'- {dna_seq} -3'")
    print(f"    [Parental]    3'- {dna_complement} -5'")
    print("=" * 60)


if __name__ == "__main__":
    main()
