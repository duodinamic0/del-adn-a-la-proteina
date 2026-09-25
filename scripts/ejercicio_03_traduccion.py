"""
Ejercicio 3: Traducción del ARNm a proteína y Extensión con Biopython.

Este script traduce un transcrito de ARNm dado a su secuencia aminoacídica
mediante el módulo Bio.Seq de Biopython, compara el resultado con la resolución
manual teórica y simula computacionalmente el impacto de mutaciones puntuales
en el codón de inicio (AUG -> GUG) y mutaciones de pérdida de parada (stop-loss).
"""

from Bio.Seq import Seq


def main() -> None:
    # 1. Transcrito de ARNm dado (5' -> 3')
    mrna_str = "AUGUAUGCUUAA"
    mrna_seq = Seq(mrna_str)

    print("=" * 70)
    print("EJERCICIO 3: TRADUCCIÓN DEL ARNm A PROTEÍNA - BIOPYTHON")
    print("=" * 70)
    print(f"Transcrito de ARNm (5' -> 3'): {mrna_seq}")
    print(f"Longitud: {len(mrna_seq)} ribonucleótidos (4 codones)")

    codones = [mrna_str[i:i+3] for i in range(0, len(mrna_str), 3)]
    print(f"Codones identificados: {codones}")
    print(f"  - Codón de inicio: {codones[0]} (Met / M)")
    print(f"  - Codones internos: {codones[1:-1]}")
    print(f"  - Codón de paro:   {codones[-1]} (Stop / *)")
    print("-" * 70)

    # 2. Traducción con Biopython
    # Traducción estándar (el codón de paro se representa convencionalmente como '*')
    protein_with_stop = mrna_seq.translate()
    # Traducción deteniéndose estrictamente antes del codón de parada
    protein_to_stop = mrna_seq.translate(to_stop=True)

    print("TRADUCCIÓN AUTOMATIZADA CON BIOPYTHON:")
    print(f"  Secuencia traducida (con símbolo Stop '*'): {protein_with_stop}")
    print(f"  Péptido maduro funcional (to_stop=True):    {protein_to_stop}")
    print("-" * 70)

    # 3. Comprobación contra la resolución teórica manual
    manual_protein_with_stop = "MYA*"
    manual_peptide = "MYA"

    assert str(protein_with_stop) == manual_protein_with_stop, (
        f"Discrepancia con stop: esperado {manual_protein_with_stop}, obtenido {protein_with_stop}"
    )
    assert str(protein_to_stop) == manual_peptide, (
        f"Discrepancia sin stop: esperado {manual_peptide}, obtenido {protein_to_stop}"
    )
    print(">> VERIFICACIÓN MANUAL vs BIOPYTHON: EXITOSA (COINCIDENCIA EXACTA)")
    print("   Metionina (M) - Tirosina (Y) - Alanina (A)")
    print("-" * 70)

    # 4. Simulación computacional de mutaciones
    print("EXPERIMENTO: SIMULACIÓN DE MUTACIONES EN EL MARCO DE TRADUCCIÓN:")

    # 4.1. Mutación en el codón de inicio: AUG -> GUG
    mut_inicio_str = "GUGUAUGCUUAA"
    mut_inicio_seq = Seq(mut_inicio_str)
    # Traducción canónica estándar
    prot_mut_inicio = mut_inicio_seq.translate(to_stop=True)
    print("\n   a) Mutación en el codón de inicio (AUG -> GUG):")
    print(f"      ARNm mutado:     5'- {mut_inicio_seq} -3'")
    print(f"      Péptido (código estándar): {prot_mut_inicio} (Val - Tyr - Ala)")
    print("      Nota biológica: En eucariotas anula la iniciación (requiere AUG).")
    print("      En bacterias puede actuar como inicio alternativo incorporando fMet,")
    print("      pero con una eficiencia sustancialmente reducida (10-30%).")

    # 4.2. Mutación de pérdida del codón de paro (Stop-loss): UAA -> CAA (Glutamina / Q)
    # Se añade una secuencia downstream hipotética correspondiente a la región 3'-UTR
    mut_stoploss_str = "AUGUAUGCUCAAGCCUGA"  # UAA mutó a CAA, y se continúa hasta el siguiente stop UGA
    mut_stoploss_seq = Seq(mut_stoploss_str)
    prot_mut_stoploss = mut_stoploss_seq.translate(to_stop=True)
    print("\n   b) Mutación de pérdida de parada / Stop-loss (UAA -> CAA):")
    print(f"      ARNm mutado (con 3'-UTR): 5'- {mut_stoploss_seq} -3'")
    print(f"      Péptido extendido anómalo:  {prot_mut_stoploss} (Met-Tyr-Ala-Gln-Ala)")
    print("      Nota biológica: El ribosoma continúa traduciendo hacia la región 3'-UTR,")
    print("      generando extensiones C-terminales anómalas, agregación proteica")
    print("      o activación de la vía de vigilancia 'Non-stop Decay' (NSD).")
    print("=" * 70)


if __name__ == "__main__":
    main()
