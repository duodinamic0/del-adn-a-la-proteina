"""
Ejercicio 6: Actividad integradora: del ADN a la proteína.
Pipeline bioinformático interactivo del Dogma Central de la Biología Molecular.

Este script ejecuta y monitoriza paso a paso los tres procesos fundamentales:
1. Replicación Semiconservativa del ADN (dúplex parental -> 2 dúplex hijas).
2. Transcripción de ADN a ARNm (hebra molde 3'->5' -> ARNm 5'->3').
3. Traducción de ARNm a Proteína (marco abierto de lectura -> polipéptido maduro).

Utiliza la secuencia oficial de la región codificante (CDS) del gen de la Insulina
humana (NCBI RefSeq NM_000207.3 / NP_000198.1), informando detalladamente
en cada etapa de los eventos mecanísticos y moleculares subyacentes.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class DogmaCentralPipeline:
    """Orquestador bioinformático del Dogma Central con registro detallado de eventos."""

    def __init__(self, fasta_input: Path) -> None:
        self.fasta_input = fasta_input
        self.record_dna = None
        self.dna_seq_5_3 = None
        self.dna_template_3_5 = None
        self.mrna_seq_5_3 = None
        self.protein_seq = None

    def log_step(self, step_number: int, title: str) -> None:
        print("\n" + "=" * 80)
        print(f"PASO {step_number}: {title.upper()}")
        print("=" * 80)

    def cargar_secuencia_dna(self) -> None:
        self.log_step(0, "Carga e Inspección de la Secuencia Biológica en FASTA")
        print(f"[INFO] Leyendo archivo de entrada: {self.fasta_input.name}")
        self.record_dna = SeqIO.read(self.fasta_input, "fasta")
        self.dna_seq_5_3 = self.record_dna.seq

        print(f"[METADATOS] ID de secuencia:     {self.record_dna.id}")
        print(f"[METADATOS] Descripción:        {self.record_dna.description}")
        print(f"[MOLÉCULA]  Longitud:           {len(self.dna_seq_5_3)} pares de bases (pb)")
        print(f"[MOLÉCULA]  Extremo 5' inicial: {self.dna_seq_5_3[:30]}...")
        print(f"[MOLÉCULA]  Extremo 3' final:   ...{self.dna_seq_5_3[-30:]}")

        # Comprobar marco abierto de lectura (ORF)
        assert str(self.dna_seq_5_3).startswith("ATG"), "Aviso: la secuencia no comienza con ATG"
        print("[VERIFICACIÓN] Codón de inicio ATG detectado en la posición 1-3.")
        print("[ESTADO] Cadena codificante parental cargada satisfactoriamente en memoria.")

    def replicar_adn(self) -> None:
        self.log_step(1, "Replicación Semiconservativa del ADN")
        print("[MECANISMO MOLECULAR]")
        print("  1. La ADN Helicasa desenrolla el dúplex parental rompiendo los puentes de H.")
        print("  2. Las proteínas SSB estabilizan las hebras sencillas separadas.")
        print("  3. La Primasa sintetiza cebadores de ARN iniciadores.")
        print("  4. La ADN Polimerasa sintetiza las nuevas hebras complementarias en sentido 5' -> 3'.")
        print("-" * 80)

        # Generación de la hebra complementaria (que actúa como molde natural orientada 3' -> 5')
        self.dna_template_3_5 = self.dna_seq_5_3.complement()
        # Hebra complementaria escrita en sentido estándar 5' -> 3'
        dna_rev_comp_5_3 = self.dna_seq_5_3.reverse_complement()

        print("[PROCESAMIENTO BIOINFORMÁTICO]")
        print(f"  > Hebra parental (codificante 5'->3'):   {self.dna_seq_5_3[:40]}...-3'")
        print(f"  > Hebra complementaria (molde 3'->5'):  {self.dna_template_3_5[:40]}...-5'")
        print(f"  > Hebra complementaria (reversa 5'->3'): {dna_rev_comp_5_3[:40]}...-3'")
        print("\n[RESULTADO DE LA REPLICACIÓN (1 Ronda Semiconservativa)]")
        print("  Se originan dos moléculas dúplex hijas idénticas al dúplex parental:")
        print("  • Molécula Hija 1:")
        print(f"      [Parental]    5'- {self.dna_seq_5_3[:30]}... -3'")
        print(f"      [Neosíntesis] 3'- {self.dna_template_3_5[:30]}... -5'")
        print("  • Molécula Hija 2:")
        print(f"      [Neosíntesis] 5'- {self.dna_seq_5_3[:30]}... -3'")
        print(f"      [Parental]    3'- {self.dna_template_3_5[:30]}... -5'")
        print("[ESTADO] Replicación completada con 100% de fidelidad de bases.")

    def transcribir_a_arnm(self, output_fasta: Path) -> None:
        self.log_step(2, "Transcripción de ADN a ARN Mensajero (ARNm)")
        print("[MECANISMO MOLECULAR]")
        print("  1. La ARN Polimerasa II dependiente de ADN se posiciona sobre el promotor.")
        print("  2. La enzima reconoce la hebra molde (3' -> 5') y avanza sintetizando en 5' -> 3'.")
        print("  3. Cada adenina del molde incorpora un uracilo (U), y cada timina una adenina (A).")
        print("-" * 80)

        # En Biopython, transcribe() toma la hebra codificante (5'->3') y sustituye T por U
        self.mrna_seq_5_3 = self.dna_seq_5_3.transcribe()

        print("[PROCESAMIENTO BIOINFORMÁTICO]")
        print(f"  > Cadena molde leída (3'->5'):  {self.dna_template_3_5[:40]}...")
        print(f"  > Transcrito primario (5'->3'): {self.mrna_seq_5_3[:40]}...")
        print(f"  > Longitud del ARNm:            {len(self.mrna_seq_5_3)} ribonucleótidos")

        # Comprobación de sustitución T -> U
        assert "T" not in str(self.mrna_seq_5_3), "Error: Se detectaron residuos de Timina en el ARN"
        assert str(self.mrna_seq_5_3).startswith("AUG"), "Error: El transcrito no inicia en AUG"

        # Guardar en archivo FASTA de salida
        mrna_record = SeqRecord(
            self.mrna_seq_5_3,
            id=f"{self.record_dna.id}_mRNA",
            description="ARNm de Insulina humana generado mediante pipeline Biopython",
        )
        SeqIO.write(mrna_record, output_fasta, "fasta")
        print(f"[SALIDA] Archivo de ARNm guardado en: {output_fasta.name}")
        print("[ESTADO] Transcripción finalizada exitosamente.")

    def traducir_a_proteina(self, output_fasta: Path) -> None:
        self.log_step(3, "Traducción de ARNm a Proteína Funcional")
        print("[MECANISMO MOLECULAR]")
        print("  1. La subunidad pequeña 40S con Met-tRNA escanea el 5'-UTR hasta hallar el codón AUG.")
        print("  2. Se acopla la subunidad 60S formando el ribosoma 80S competente.")
        print("  3. En cada ciclo de elongación (sitios A -> P -> E) se añade un aminoácido.")
        print("  4. Al encontrar un codón de parada (UAG/UAA/UGA), factores de liberación liberan el péptido.")
        print("-" * 80)

        # Traducción con Biopython
        self.protein_seq = self.mrna_seq_5_3.translate(to_stop=True)

        print("[PROCESAMIENTO BIOINFORMÁTICO]")
        print(f"  > Péptido traducido: {self.protein_seq}")
        print(f"  > Longitud:          {len(self.protein_seq)} aminoácidos")
        print(f"  > Extremo N-term:    {self.protein_seq[0]} (Metionina iniciadora)")
        print(f"  > Extremo C-term:    {self.protein_seq[-1]} (Asparagina terminal de la cadena A)")

        # Comprobación con la proteína canónica humana NP_000198.1 (110 aa)
        secuencia_esperada_insulina = (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTR"
            "REAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        )
        assert str(self.protein_seq) == secuencia_esperada_insulina, (
            "Discrepancia en la traducción con la preproinsulina humana canónica NP_000198.1"
        )
        print("[VALIDACIÓN] ¡Concordancia exacta al 100% con la Preproinsulina humana (NP_000198.1)!")

        # Estructura de dominios biológicos de la preproinsulina
        print("\n[DOMINIOS BIOLÓGICOS DE LA PREPROINSULINA]")
        print(f"  • Péptido Señal (aa 1-24):   {self.protein_seq[0:24]} (Direccionamiento al RE)")
        print(f"  • Cadena B (aa 25-54):       {self.protein_seq[24:54]} (Cadena B madura)")
        print(f"  • Péptido C conector (55-89): {self.protein_seq[54:89]} (Escindido por prohormona convertasas)")
        print(f"  • Cadena A (aa 90-110):      {self.protein_seq[89:110]} (Cadena A madura unida por disulfuros)")

        # Guardar en archivo FASTA de salida
        prot_record = SeqRecord(
            self.protein_seq,
            id="NP_000198.1_Insulina_Humana",
            description="Preproinsulina humana madura traducida mediante pipeline Biopython",
        )
        SeqIO.write(prot_record, output_fasta, "fasta")
        print(f"[SALIDA] Archivo de Proteína guardado en: {output_fasta.name}")
        print("[ESTADO] Traducción completada con éxito.")

    def ejecutar_pipeline_completo(self, output_mrna: Path, output_prot: Path) -> None:
        print("#" * 80)
        print("INICIO DEL PIPELINE BIOPYTHON: DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR")
        print("#" * 80)
        self.cargar_secuencia_dna()
        self.replicar_adn()
        self.transcribir_a_arnm(output_mrna)
        self.traducir_a_proteina(output_prot)
        print("\n" + "#" * 80)
        print("PIPELINE EJECUTADO CON ÉXITO: ADN -> ARN -> PROTEÍNA")
        print("#" * 80)


def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    fasta_input = base_dir / "data" / "insulina_humana_cds.fasta"
    output_mrna = base_dir / "data" / "ejercicio_06_insulina_arnm.fasta"
    output_prot = base_dir / "data" / "ejercicio_06_insulina_proteina.fasta"

    pipeline = DogmaCentralPipeline(fasta_input)
    pipeline.ejecutar_pipeline_completo(output_mrna, output_prot)


if __name__ == "__main__":
    main()
