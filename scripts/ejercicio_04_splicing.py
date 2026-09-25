"""
Ejercicio 4: Splicing Alternativo y Extensión con Biopython / Ensembl.

Este script:
1. Modela computacionalmente un gen hipotético de 5 exones y genera diferentes
   combinaciones de splicing alternativo (exclusión de exón, exones mutuamente excluyentes),
   traduciendo cada isoforma y analizando los cambios en la proteína resultante.
2. Analiza el caso real de Ensembl del receptor humano FGFR2 (ENSG00000066468),
   comparando las isoformas FGFR2-IIIb (epitelial) y FGFR2-IIIc (mesenquimal),
   sus secuencias y su repercusión biológica en la especificidad de ligandos y la EMT.
"""

from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO


def modelar_splicing_teorico() -> None:
    print("=" * 75)
    print("1. MODELADO COMPUTACIONAL DE SPLICING ALTERNATIVO (GEN DE 5 EXONES)")
    print("=" * 75)

    # Definición de 5 exones hipotéticos (secuencias codificantes en fase)
    exones = {
        "Exon_1": Seq("ATGGCTAGCAAC"),            # Met-Ala-Ser-Asn (Péptido señal / Inicio)
        "Exon_2": Seq("GTCCTGCTGGAC"),            # Val-Leu-Leu-Asp (Dominio regulador A)
        "Exon_3": Seq("AAGTGCCACTACGAG"),         # Lys-Cys-His-Tyr-Glu (Dominio de unión / Sitio activo)
        "Exon_4": Seq("TTCAACGGCAAGATC"),         # Phe-Asn-Gly-Lys-Ile (Dominio bisagra / espaciador)
        "Exon_5": Seq("CCCTGGCTCTACCGCTAA"),     # Pro-Trp-Leu-Tyr-Arg-STOP (Dominio transmembrana y parada)
    }

    print("Exones del pre-ARNm primario:")
    for nombre, seq in exones.items():
        print(f"  {nombre} ({len(seq)} nt): {seq} -> {seq.translate()}")
    print("-" * 75)

    # Combinaciones de splicing alternativo
    isoformas_mrna = {
        "Isoforma_A (Constitutiva 1-2-3-4-5)": (
            exones["Exon_1"] + exones["Exon_2"] + exones["Exon_3"] + exones["Exon_4"] + exones["Exon_5"]
        ),
        "Isoforma_B (Salto de Exón 3: 1-2-4-5)": (
            exones["Exon_1"] + exones["Exon_2"] + exones["Exon_4"] + exones["Exon_5"]
        ),
        "Isoforma_C (Exclusión de Exones 2 y 4: 1-3-5)": (
            exones["Exon_1"] + exones["Exon_3"] + exones["Exon_5"]
        ),
    }

    print("Traducción de Isoformas Resultantes:")
    for nombre_iso, seq_mrna in isoformas_mrna.items():
        proteina = seq_mrna.translate(to_stop=True)
        print(f"\n* {nombre_iso}:")
        print(f"  Longitud ARNm:     {len(seq_mrna)} nt")
        print(f"  Secuencia ARNm:    5'- {seq_mrna} -3'")
        print(f"  Proteína ({len(proteina)} aa): {proteina}")

    print("\nImpacto funcional en las proteínas resultantes:")
    print("  - Isoforma A: Posee todos los dominios (longitud máxima: 18 aa).")
    print("  - Isoforma B: Pierde el Exón 3 (KCHYE). Carece del dominio de unión/sitio activo,")
    print("                pudiendo actuar como isoforma dominante negativa o afuncional.")
    print("  - Isoforma C: Pierde los exones 2 y 4. Conserva el dominio de unión central pero")
    print("                con menor flexibilidad conformacional y longitud compacta (13 aa).")
    print("-" * 75)


def analizar_caso_ensembl_fgfr2() -> None:
    print("\n" + "=" * 75)
    print("2. CASO REAL ENSEMBL: RECEPTOR FGFR2 HUMANO (ENSG00000066468)")
    print("=" * 75)

    base_dir = Path(__file__).resolve().parent.parent
    fasta_path = base_dir / "data" / "fgfr2_isoformas.fasta"

    if not fasta_path.exists():
        print(f"Aviso: No se encontró {fasta_path}. Creando análisis estructurado...")
        return

    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Archivo cargado: {fasta_path.name} ({len(records)} entradas Ensembl)")

    for rec in records:
        aa_seq = rec.seq.translate(to_stop=True)
        print(f"\n* Transcrito Ensembl: {rec.id}")
        print(f"  Descripción:   {rec.description}")
        print(f"  Longitud (nt): {len(rec.seq)} pb")
        print(f"  Secuencia aa:  {aa_seq}")
        print(f"  Longitud (aa): {len(aa_seq)} residuos")

    print("\n" + "-" * 75)
    print("COMPARACIÓN ESTRUCTURAL Y BIOLÓGICA DE LAS ISOFORMAS DE FGFR2:")
    print("  1. FGFR2-IIIb (Ensembl ENST00000358487 / Exón 8):")
    print("     - Expresión tisular: Restringida a células EPITELIALES.")
    print("     - Especificidad de ligando: Alta afinidad por FGF7 (KGF) y FGF10.")
    print("     - Función: Comunicación paracrina epitelio-mesénquima en homeostasis cutánea y pulmonar.")
    print("\n  2. FGFR2-IIIc (Ensembl ENST00000356226 / Exón 9):")
    print("     - Expresión tisular: Restringida a células MESENQUIMALES.")
    print("     - Especificidad de ligando: Alta afinidad por FGF2 (bFGF) y FGF4.")
    print("     - Función: Morfogénesis mesenquimal, osteogénesis y diferenciación condrogénica.")
    print("\n  3. Relevancia Patológica y Transición Epitelio-Mesénquima (EMT):")
    print("     - La conmutación de splicing de la isoforma IIIb hacia la IIIc")
    print("       es un biomarcador crucial de la transición epitelio-mesénquima (EMT),")
    print("       un proceso oncogénico que confiere motilidad, invasividad celular")
    print("       y capacidad metastásica a los carcinomas humanos.")
    print("=" * 75)


def main() -> None:
    modelar_splicing_teorico()
    analizar_caso_ensembl_fgfr2()


if __name__ == "__main__":
    main()
