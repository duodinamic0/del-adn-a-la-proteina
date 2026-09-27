"""
Ejercicio 5: Introducción a las proteínas - Secuencia, Estructura y Función.

Este script:
1. Analiza el péptido problema (Met-Ile-Ser-Gly-Val-Lys-His), identificando sus
   extremos N y C terminales, sus propiedades fisicoquímicas, masa molecular y pI.
2. Analiza la estructura tridimensional de una proteína modelo del Protein Data Bank
   (Lisozima de Gallus gallus, PDB ID: 1AKI) mediante Bio.PDB de Biopython, cuantificando
   la distribución de hélices alfa, láminas beta y bucles.
3. Evalúa la distribución radial de residuos hidrofóbicos vs. hidrofílicos respecto al
   centro de masa para ilustrar el efecto hidrofóbico en el empaquetamiento del núcleo.
4. Modela el impacto biofísico de mutaciones puntuales desestabilizadoras (hidrofóbico -> hidrofílico).
"""

from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def analizar_peptido_problema() -> None:
    print("=" * 75)
    print("1. ANÁLISIS DEL PÉPTIDO PROBLEMA (Met - Ile - Ser - Gly - Val - Lys - His)")
    print("=" * 75)

    seq_1letter = "MISGVKH"
    nombres_3letter = ["Met", "Ile", "Ser", "Gly", "Val", "Lys", "His"]

    print(f"Secuencia: {' - '.join(nombres_3letter)} (Código 1 letra: {seq_1letter})")
    print(f"Extremo N-terminal (Amino-terminal libre, H3N+-): {nombres_3letter[0]} (Metionina)")
    print(f"Extremo C-terminal (Carboxilo-terminal libre, -COO-): {nombres_3letter[-1]} (Histidina)")
    print(f"Estructura química: H3N+-Met-Ile-Ser-Gly-Val-Lys-His-COO-")
    print("-" * 75)

    # Clasificación fisicoquímica de los aminoácidos
    propiedades = [
        ("Met", "Apolar alifático (azufrado)", "Iniciador, apolar hidrofóbico"),
        ("Ile", "Apolar alifático ramificado", "Altamente hidrofóbico, empaquetamiento"),
        ("Ser", "Polar sin carga (hidroxilado)", "Formador de puentes de H, fosforilable"),
        ("Gly", "Apolar especial (sin C-beta)", "Aquiral, máxima flexibilidad conformacional"),
        ("Val", "Apolar alifático ramificado", "Hidrofóbico, favorece láminas beta"),
        ("Lys", "Básico con carga positiva (+1)", "Hidrofílico superficial, puentes salinos"),
        ("His", "Básico aromático (anillo imidazol)", "pKa ~ 6.0, sensor de pH y catálisis"),
    ]

    print("Propiedades fisicoquímicas de los residuos del heptapéptido:")
    for aa, tipo, rol in propiedades:
        print(f"  - {aa}: {tipo:<35} | {rol}")

    # Análisis con Bio.SeqUtils.ProtParam
    analisis = ProteinAnalysis(seq_1letter)
    print("\nParámetros calculados con Biopython (ProtParam):")
    print(f"  - Masa molecular estimada: {analisis.molecular_weight():.2f} Da")
    print(f"  - Punto isoeléctrico teórico (pI): {analisis.isoelectric_point():.2f}")
    print(f"  - Índice de hidropatía promedio (GRAVY): {analisis.gravy():.2f}")
    print("-" * 75)


def analizar_estructura_pdb() -> None:
    print("\n" + "=" * 75)
    print("2. ANÁLISIS ESTRUCTURAL DE LISOZIMA (PDB ID: 1AKI) CON Bio.PDB")
    print("=" * 75)

    base_dir = Path(__file__).resolve().parent.parent
    pdb_path = base_dir / "data" / "1aki.pdb"

    if not pdb_path.exists():
        print(f"Error: no se localizó el archivo PDB en {pdb_path}")
        return

    # 1. Parsear estructura cristalográfica con PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1AKI", str(pdb_path))
    model = structure[0]
    chain_a = model["A"]

    residues = [r for r in chain_a.get_residues() if r.id[0] == " "]
    total_res = len(residues)
    print(f"Proteína analizada: Lisozima de gallina (Gallus gallus)")
    print(f"Cadena: A | Longitud total: {total_res} aminoácidos")

    # 2. Extracción de elementos secundarios de los registros HELIX y SHEET del archivo PDB
    helices_res = set()
    sheets_res = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HELIX"):
                start_res = int(line[21:25].strip())
                end_res = int(line[33:37].strip())
                for r in range(start_res, end_res + 1):
                    helices_res.add(r)
            elif line.startswith("SHEET"):
                start_res = int(line[22:26].strip())
                end_res = int(line[33:37].strip())
                for r in range(start_res, end_res + 1):
                    sheets_res.add(r)

    n_helix = len(helices_res)
    n_sheet = len(sheets_res)
    n_loop = total_res - (n_helix + n_sheet)

    print("\nComposición de Estructura Secundaria:")
    print(f"  - Hélices alfa y 3_10: {n_helix} residuos ({n_helix / total_res * 100:.1f}%)")
    print(f"  - Láminas beta:        {n_sheet} residuos ({n_sheet / total_res * 100:.1f}%)")
    print(f"  - Bucles y giros:      {n_loop} residuos ({n_loop / total_res * 100:.1f}%)")
    print("-" * 75)

    # 3. Análisis del Centro de Masa y Distribución Radial (Efecto Hidrofóbico)
    ca_coords = []
    hydrophobic_dists = []
    hydrophilic_dists = []

    res_hidrofobicos = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"}
    res_hidrofilicos_cargados = {"ARG", "LYS", "ASP", "GLU"}

    for r in residues:
        if "CA" in r:
            ca_coords.append(r["CA"].coord)

    center_of_mass = np.mean(ca_coords, axis=0)

    for r in residues:
        if "CA" in r:
            dist = np.linalg.norm(r["CA"].coord - center_of_mass)
            res_name = r.get_resname()
            if res_name in res_hidrofobicos:
                hydrophobic_dists.append(dist)
            elif res_name in res_hidrofilicos_cargados:
                hydrophilic_dists.append(dist)

    mean_dist_hydrophobic = np.mean(hydrophobic_dists)
    mean_dist_hydrophilic = np.mean(hydrophilic_dists)

    print("Distribución Espacial de Residuos (Efecto Hidrofóbico):")
    print(f"  - Centro de masa geométrico C-alfa: {np.round(center_of_mass, 2)}")
    print(f"  - Distancia promedio de residuos hidrofóbicos (Leu, Ile, Val, Phe, Met, Trp): {mean_dist_hydrophobic:.2f} Å")
    print(f"  - Distancia promedio de residuos cargados hidrofílicos (Arg, Lys, Asp, Glu):   {mean_dist_hydrophilic:.2f} Å")
    print(f"  >> Diferencia radial: {mean_dist_hydrophilic - mean_dist_hydrophobic:.2f} Å hacia la superficie")
    print("  >> Conclusión: Los residuos hidrofóbicos se concentran preferentemente en el")
    print("     núcleo interno anhidro, mientras los polares/cargados se proyectan hacia el solvente.")
    print("-" * 75)

    # 4. Simulación de mutación puntual interna: Ile98 -> Glu98
    print("3. SIMULACIÓN BIOFÍSICA DE MUTACIÓN PUNTUAL EN EL NÚCLEO:")
    print("   Caso: Isoleucina 98 (residuo hidrofóbico enterrado en la hélice alpha 5 de 1AKI)")
    r98 = chain_a[98]
    dist_r98 = np.linalg.norm(r98["CA"].coord - center_of_mass)
    print(f"   - Posición 98: {r98.get_resname()}98 a {dist_r98:.2f} Å del centro de masa (núcleo interno)")
    print("   - Mutación hipotética: Ile98Glu (Apolar hidrofóbico -> Polar ácido cargado negativamente)")
    print("   - Consecuencias moleculares esperadas:")
    print("     1. Penalización termodinámica de desolvatación (~15 a 20 kcal/mol) al enterrar un carboxilo (-COO-).")
    print("     2. Ruptura del empaquetamiento estérico de van der Waals de las cadenas laterales circundantes.")
    print("     3. Desestabilización drástica de la energía libre de plegamiento (Delta-Delta G < 0).")
    print("     4. Colapso del plegamiento nativo, desnaturalización parcial y agregación citotóxica.")
    print("=" * 75)


def main() -> None:
    analizar_peptido_problema()
    analizar_estructura_pdb()


if __name__ == "__main__":
    main()
