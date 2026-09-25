"""
Ejercicio 2: Transcripción del ADN a ARN y Extensión con Biopython.

Este script lee una secuencia de ADN desde un archivo en formato FASTA,
obtiene su correspondiente transcrito de ARNm mediante Biopython,
y realiza experimentos sistemáticos alterando la orientación de la hebra
para evaluar el impacto biológico y bioinformático en la secuencia transcrita.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main() -> None:
    # Localizar archivo FASTA de entrada
    base_dir = Path(__file__).resolve().parent.parent
    fasta_path = base_dir / "data" / "ejercicio_02_adn.fasta"

    if not fasta_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo FASTA en: {fasta_path}")

    print("=" * 70)
    print("EJERCICIO 2: TRANSCRIPCIÓN DEL ADN A ARN - EXTENSIÓN CON BIOPYTHON")
    print("=" * 70)

    # 1. Lectura del archivo FASTA con Bio.SeqIO
    record = SeqIO.read(fasta_path, "fasta")
    print(f"Archivo FASTA cargado:  {fasta_path.name}")
    print(f"Identificador (ID):     {record.id}")
    print(f"Descripción:            {record.description}")
    print(f"Longitud de secuencia:  {len(record.seq)} nucleótidos")
    print(f"Secuencia ADN (5'->3'): {record.seq}")
    print("-" * 70)

    # 2. Transcripción directa (convención estándar de Biopython)
    # En Biopython, se asume que la secuencia provista en el FASTA corresponde a la
    # hebra codificante (sentido, 5' -> 3'). Por ende, .transcribe() sustituye T por U.
    mrna_direct = record.seq.transcribe()
    print("1. TRANSCRIPCIÓN DIRECTA (Hebra codificante 5' -> 3'):")
    print(f"   ADN codificante: 5'- {record.seq} -3'")
    print(f"   ARNm transcrito: 5'- {mrna_direct} -3'")
    print(f"   Codones:         {list(map(''.join, zip(*[iter(str(mrna_direct))]*3)))}")
    print("-" * 70)

    # 3. Experimento: Cambio de orientación de la hebra
    print("2. EXPERIMENTO: ALTERACIÓN DE LA ORIENTACIÓN DE LA HEBRA:")
    
    # 3.1. Caso Reversa Complementaria (Hebra antisentido 5' -> 3')
    # Representa la hebra opuesta leída en su propio sentido 5' -> 3'.
    rev_comp_dna = record.seq.reverse_complement()
    mrna_rev_comp = rev_comp_dna.transcribe()
    print("   a) Transcripción desde la hebra inversa-complementaria (5' -> 3'):")
    print(f"      ADN (rev-comp):  5'- {rev_comp_dna} -3'")
    print(f"      ARNm resultante: 5'- {mrna_rev_comp} -3'")
    print(f"      Codones:         {list(map(''.join, zip(*[iter(str(mrna_rev_comp))]*3)))}")

    # 3.2. Caso Inversión pura (Secuencia leída de 3' a 5' sin complementar)
    reversed_dna = record.seq[::-1]
    mrna_reversed = reversed_dna.transcribe()
    print("\n   b) Transcripción tras invertir el orden físico (3' -> 5' tratada como 5' -> 3'):")
    print(f"      ADN invertido:   5'- {reversed_dna} -3'")
    print(f"      ARNm resultante: 5'- {mrna_reversed} -3'")
    print(f"      Codones:         {list(map(''.join, zip(*[iter(str(mrna_reversed))]*3)))}")
    print("-" * 70)

    # 4. Conclusiones biológicas del experimento
    print("3. ANÁLISIS DEL EXPERIMENTO:")
    print("   - La ARN polimerasa requiere una polaridad estricta (5' -> 3').")
    print("   - Invertir la orientación o seleccionar la hebra contraria cambia radicalmente")
    print("     la pauta de lectura (marco de lectura), la composición de codones y el")
    print("     péptido codificado, eliminando los codones de inicio o generando codones aberrantes.")
    print("-" * 70)

    # 5. Exportar el ARNm primario a un nuevo archivo FASTA
    output_fasta = base_dir / "data" / "ejercicio_02_arnm.fasta"
    mrna_record = SeqRecord(
        mrna_direct,
        id=f"{record.id}_mRNA",
        description="ARNm transcrito primario (5'->3') obtenido con Biopython",
    )
    SeqIO.write(mrna_record, output_fasta, "fasta")
    print(f">> Archivo de ARNm generado exitosamente en: {output_fasta.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()
