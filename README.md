# Del ADN a la Proteína: Fundamentos de Biología Molecular y Biopython

Este repositorio compila una serie de ejercicios teórico-prácticos centrados en los procesos fundamentales del dogma central de la biología molecular: **replicación del ADN**, **transcripción** y **traducción**, complementados con análisis bioinformáticos automatizados mediante la librería [Biopython](https://biopython.org/).

---

## 🛠️ Configuración del Entorno y Reproducibilidad

El proyecto gestiona sus dependencias y entornos virtuales mediante [`uv`](https://docs.astral.sh/uv/), un gestor ultrarrápido de paquetes para Python.

### Requisitos previos
- Tener instalado `uv` en el sistema.

### Instalación y ejecución
Para sincronizar las dependencias e inicializar el entorno:
```bash
uv sync
```

Para ejecutar cualquiera de los scripts implementados:
```bash
# Ejercicio 1: Replicación
uv run scripts/ejercicio_01_replicacion.py

# Ejercicio 2: Transcripción
uv run scripts/ejercicio_02_transcripcion.py

# Ejercicio 3: Traducción
uv run scripts/ejercicio_03_traduccion.py

# Ejercicio 4: Splicing alternativo
uv run scripts/ejercicio_04_splicing.py
```

---

## Índice de Ejercicios
1. [Ejercicio 1: Replicación del ADN](#ejercicio-1-replicación-del-adn)
2. [Ejercicio 2: Transcripción del ADN a ARN](#ejercicio-2-transcripción-del-adn-a-arn)
3. [Ejercicio 3: Traducción del ARNm a proteína](#ejercicio-3-traducción-del-arnm-a-proteína)
4. [Ejercicio 4: Splicing alternativo](#ejercicio-4-splicing-alternativo)

---

## Ejercicio 1: Replicación del ADN

### 1. Objetivo
Comprender en profundidad el mecanismo **semiconservativo** de la replicación del ADN, la polaridad y complementariedad de las hebras, la función de las enzimas centrales del replisoma y el impacto biológico de las mutaciones puntuales derivadas de fallos de fidelidad en la síntesis.

---

### 2. Planteamiento y Análisis de la Secuencia

Se suministra el siguiente dúplex de ADN bicatenario:

$$\begin{aligned}
\text{Hebra superior (sentido / directa):} \quad & 5'\text{ – ATG CCG TTA GCT – }3' \\
\text{Hebra inferior (molde / complementaria):} \quad & 3'\text{ – TAC GGC AAT CGA – }5'
\end{aligned}$$

Ambas cadenas se mantienen unidas mediante puentes de hidrógeno específicos entre pares de bases nitrogenadas según las reglas de complementariedad de Watson y Crick (adenina con timina mediante 2 puentes de H; citosina con guanina mediante 3 puentes de H) y presentan orientación **antiparalela**.

---

### 3. Ronda de Replicación: Mecanismo Semiconservativo

De acuerdo con el modelo demostrado experimentalmente por **Matthew Meselson y Franklin Stahl (1958)**, la replicación del ADN es **semiconservativa**: durante la duplicación celular, las dos cadenas del dúplex parental se desenrollan y se separan, sirviendo cada una como molde directriz para la polimerización de una nueva hebra complementaria.

Al finalizar una ronda completa de replicación, se obtienen **dos moléculas hijas bicatenarias idénticas**, cada una compuesta por una hebra parental original y una hebra de nueva síntesis (neosintetizada).

#### Identificación de las nuevas hebras y dúplex resultantes

1. **Molécula Hija 1 (originada a partir de la hebra parental $5' \to 3'$):**
   - **Hebra parental (molde):** $5'\text{ – ATG CCG TTA GCT – }3'$
   - **Nueva hebra neosintetizada:** $3'\text{ – TAC GGC AAT CGA – }5'$
   *(Nótese que, escrita bajo la convención estándar internacional $5' \to 3'$, la hebra neosintetizada corresponde a $5'\text{ – AGC TAA CGG CAT – }3'$).*

2. **Molécula Hija 2 (originada a partir de la hebra parental $3' \to 5'$):**
   - **Hebra parental (molde):** $3'\text{ – TAC GGC AAT CGA – }5'$
   - **Nueva hebra neosintetizada:** $5'\text{ – ATG CCG TTA GCT – }3'$

#### Esquema del resultado tras una ronda replicativa:

```text
Dúplex Parental Original:
  5'- ATG CCG TTA GCT -3'
  3'- TAC GGC AAT CGA -5'
         │
         ▼  (Desenrollamiento y síntesis por ADN polimerasa)
  ┌─────────────────────────────────┬─────────────────────────────────┐
  │         Molécula Hija 1         │         Molécula Hija 2         │
  ├─────────────────────────────────┼─────────────────────────────────┤
  │ [Parental] 5'-ATG CCG TTA GCT-3'│ [Neosíntesis] 5'-ATG CCG TTA GCT-3'│
  │ [Neosínt.] 3'-TAC GGC AAT CGA-5'│ [Parental]    3'-TAC GGC AAT CGA-5'│
  └─────────────────────────────────┴─────────────────────────────────┘
```

Ambos dúplex conservan con total fidelidad la información génica original.

---

### 4. Función de las Enzimas Implicadas en el Replisoma

La replicación es un proceso finamente orquestado por un complejo macromolecular denominado **replisoma**. A continuación se detalla el papel fisiológico de las enzimas clave:

| Enzima | Clasificación funcional | Función molecular en el proceso |
| :--- | :--- | :--- |
| **ADN Helicasa** | ATPasa dependiente de ADN (*e.g.*, DnaB en bacterias, complejo MCM2-7 en eucariotas) | Se une al origen de replicación y se desplaza rompiendo los puentes de hidrógeno entre las bases nitrogenadas. Su acción desestabiliza y desnaturaliza localmente la doble hélice, generando y haciendo avanzar la **horquilla de replicación**. |
| **ADN Primasa** | ARN polimerasa dependiente de ADN (*e.g.*, DnaG en bacterias, complejo Pol $\alpha$/Primasa en eucariotas) | Las ADN polimerasas son incapaces de iniciar la síntesis *de novo* porque requieren obligatoriamente un extremo $3'\text{-OH}$ libre preexistente. La primasa sintetiza un oligonucleótido corto de ARN (**cebador** o *primer*, de $\sim 10\text{-}12$ nucleótidos) complementario a la hebra molde, facilitando dicho extremo nucleofílico iniciador. |
| **ADN Polimerasa** | Nucleotidiltransferasa (*e.g.*, Pol III / Pol I en procariotas; Pol $\delta$ / $\epsilon$ en eucariotas) | Cataliza la adición secuencial de desoxirribonucleótidos trifosfato (dNTPs) complementarios a la hebra molde en sentido estrictamente **$5' \to 3'$**, liberando pirofosfato inorgánico ($\text{PP}_i$). Además, posee una actividad catalítica intrínseca de **exonucleasa $3' \to 5'$ de corrección de pruebas** (*proofreading*), que le permite retroceder, escindir nucleótidos mal apareados y reanudar la síntesis correcta. |
| **ADN Ligasa** | Ligasa dependiente de ATP o $\text{NAD}^+$ | Sella las discontinuidades covalentes monofilamentarias (*nicks*) en el esqueleto fosfodiéster. Concretamente, cataliza la formación del enlace fosfodiéster entre el grupo $3'\text{-OH}$ libre de un fragmento contiguo y el grupo $5'\text{-fosfato}$ del siguiente, siendo indispensable para unir los **fragmentos de Okazaki** en la hebra retardada y sellar el ADN tras la sustitución de los cebadores de ARN por ADN. |

---

### 5. Reflexión Biológica: Consecuencias de un Error No Corregido por la ADN Polimerasa

#### Mecanismos de salvaguarda de la fidelidad
En condiciones fisiológicas normales, la ADN polimerasa exhibe una tasa intrínseca de incorporación errónea de aproximadamente $10^{-4}$ a $10^{-5}$ por nucleótido incorporado. Sin embargo, gracias a su actividad **exonucleasa $3' \to 5'$ de autocorrección (*proofreading*)**, esta tasa desciende a $\sim 10^{-7}$. Posteriormente, el sistema enzimático posreplicativo de reparación de desapareamientos (**MMR**, *Mismatch Repair*) detecta y subsana la mayoría de las distorsiones residuales en la hebra hija recién sintetizada, alcanzando una tasa de error global final de aproximadamente $10^{-9}$ a $10^{-10}$ mutaciones por par de bases y por división celular.

#### ¿Qué sucede si la polimerasa comete un error y no se corrige?

Si tanto la actividad exonucleasa correctora como el sistema MMR fallan en corregir una base mal apareada (por ejemplo, incorporando una Guanina frente a una Timina):

1. **Aparición de un desapareamiento transitorio (*mismatch*):**
   Durante la ronda actual, existirá una distorsión estructural en el dúplex de ADN ($\text{T-G}$ en lugar de $\text{T-A}$).

2. **Fijación permanente de la mutación tras la siguiente ronda de replicación:**
   Al dividirse nuevamente la célula, las hebras del dúplex con el error se separarán como moldes independientes:
   - La hebra que portaba la base original ($\text{T}$) servirá de molde para sintetizar una molécula hija normal con el par $\text{T-A}$.
   - La hebra que incorporó la base errónea ($\text{G}$) servirá de molde para incorporar una Citosina ($\text{C}$), dando lugar a una molécula hija portadora de un par permanente $\text{G-C}$. En este punto, **el error se ha fijado covalentemente como una mutación puntual** y ya no es reconocible como un daño por los sistemas de reparación.

3. **Repercusiones biológicas y celulares:**
   - **En regiones no codificantes / intergénicas:** Puede carecer de consecuencias fenotípicas (mutación neutra), o bien alterar motivos de unión de factores de transcripción o regiones promotoras afectando la expresión génica.
   - **En regiones codificantes (exones):**
     - **Mutación silenciosa (sinónima):** Debido a la degeneración del código genético, el nuevo codón puede codificar para el mismo aminoácido, no alterando la proteína.
     - **Mutación de sentido erróneo (*missense*):** Sustituye un aminoácido por otro. Dependiendo de la naturaleza fisicoquímica del residuo y de si afecta al sitio activo o al plegamiento terciario, la proteína resultante puede conservar actividad parcial, perder completamente su función o adquirir propiedades deletéreas.
     - **Mutación sin sentido (*nonsense*):** La alteración puede transformar un codón con sentido en un codón de parada prematuro ($\text{UAA}$, $\text{UAG}$ o $\text{UGA}$ en el ARNm), truncando la traducción y generando proteínas afuncionales que suelen ser degradadas.
   - **Patología y evolución:** En células somáticas de organismos pluricelulares, la acumulación de mutaciones en protooncogenes (como *RAS*) o genes supresores de tumores (como *TP53*) conduce a la transformación oncogénica y desarrollo tumoral (como ocurre en el síndrome de Lynch, causado precisamente por defectos en los genes del sistema MMR: *MSH2*, *MLH1*). Por el contrario, en células de la línea germinal, una mutación puede transmitirse a la descendencia constituyendo el sustrato primario de la variabilidad genética y la evolución biológica.

---

### 6. Extensión Práctica con Biopython (Ejercicio 1)

Para complementar el análisis teórico, se implementó el script [`scripts/ejercicio_01_replicacion.py`](scripts/ejercicio_01_replicacion.py) utilizando el módulo `Bio.Seq` de Biopython.

#### Distinción bioinformática: `.complement()` vs. `.reverse_complement()`
- El método `dna_seq.complement()` reemplaza cada nucleótido por su homólogo de Watson-Crick ($\text{A} \leftrightarrow \text{T}$, $\text{C} \leftrightarrow \text{G}$), conservando los índices en la misma posición relativa; por ende, representa la hebra orientada en sentido **$3' \to 5'$**:
  $$\text{Original } (5'\to 3'): \text{ATG CCG TTA GCT} \quad \implies \quad \text{Complementaria } (3'\to 5'): \text{TAC GGC AAT CGA}$$
- El método `dna_seq.reverse_complement()` realiza tanto la complementariedad de bases como la inversión de la cadena, presentando la secuencia en el sentido biológico estándar **$5' \to 3'$**:
  $$\text{Reverse Complement } (5'\to 3'): \text{AGC TAA CGG CAT}$$

#### Verificación algorítmica
El script ejecuta aserciones directas (`assert`) que comparan los resultados teóricos obtenidos manualmente con las salidas computacionales de Biopython:
```python
from Bio.Seq import Seq

dna_seq = Seq("ATGCCGTTAGCT")
manual_complement_3_5 = "TACGGCAATCGA"
manual_rev_complement_5_3 = "AGCTAACGGCAT"

assert str(dna_seq.complement()) == manual_complement_3_5
assert str(dna_seq.reverse_complement()) == manual_rev_complement_5_3
```
La ejecución arroja coincidencia exacta al $100\%$, confirmando la validez del modelado teórico.

---

## Ejercicio 2: Transcripción del ADN a ARN

### 1. Objetivo
Comprender los principios mecanísticos de la transcripción génica catalizada por la ARN polimerasa dependiente de ADN, discriminar inequívocamente entre la hebra codificante y la hebra molde en función de la direccionalidad enzimática y el marco de lectura, delimitar estructuralmente las regiones promotoras y codificantes, y evaluar bioinformáticamente el impacto de la polaridad de las cadenas en la síntesis del ARN mensajero (ARNm).

---

### 2. Planteamiento de la Secuencia y Discriminación de Cadenas

Se parte del siguiente fragmento bicatenario de ADN:

$$\begin{aligned}
\text{Hebra Superior:} \quad & 5'\text{ – ATG CCT GAA TGC – }3' \\
\text{Hebra Inferior:} \quad & 3'\text{ – TAC GGA CTT ACG – }5'
\end{aligned}$$

#### Identificación de la cadena molde y la cadena codificante

1. **Restricción catalítica de la ARN Polimerasa:**
   La ARN polimerasa sintetiza cadenas de polirribonucleótidos exclusivamente en sentido **$5' \to 3'$**, incorporando ribonucleótidos trifosfato ($\text{NTPs}$: $\text{ATP}$, $\text{UTP}$, $\text{CTP}$, $\text{GTP}$) mediante ataque nucleofílico del grupo $3'\text{-OH}$ del nucleótido previo sobre el fosfato $\alpha$ del nucleótido entrante. En consecuencia, la enzima debe deslizarse sobre una hebra molde que discurra en dirección **antiparalela**, esto es, estrictamente **$3' \to 5'$**.

2. **Correspondencia con el codón canónico de inicio:**
   En el código genético universal, la traducción se inicia a partir del codón $\text{AUG}$ ($5'\text{-AUG-}3'$), el cual codifica para N-formilmetionina (procariotas) o metionina (eucariotas). Para que un transcrito primario comience en $5'\text{-AUG-}3'$, la cadena de ADN molde sobre la que polimeriza la enzima debe contener necesariamente el triplete complementario orientado en sentido $3' \to 5'$:
   $$3'\text{-TAC-}5' \quad \xrightarrow{\text{Transcripción}} \quad 5'\text{-AUG-}3'$$

3. **Determinación:**
   - **Cadena Molde (*Template strand* / Hebra no codificante / Antisentido / Hebra $-$):**
     Es la **hebra inferior** ($3'\text{ – TAC GGA CTT ACG – }5'$). Sirve de sustrato físico directo para el apareamiento transitorio de bases con los ribonucleótidos entrantes.
   - **Cadena Codificante (*Coding strand* / Hebra sentido / Hebra $+$):**
     Es la **hebra superior** ($5'\text{ – ATG CCT GAA TGC – }3'$). Su secuencia nucleotídica es cualitativamente idéntica a la del ARNm maduro generado, con la salvedad de que contiene Desoxitimidina ($\text{dT}$) en lugar de Uridina ($\text{U}$).

---

### 3. Obtención del Transcrito de ARN Mensajero (ARNm)

Aplicando las reglas de complementariedad canónicas del ARN frente a la hebra molde de ADN:
- $\text{Adenina de ADN (A)} \implies \text{Uracilo de ARN (U)}$
- $\text{Timina de ADN (T)} \implies \text{Adenina de ARN (A)}$
- $\text{Citosina de ADN (C)} \implies \text{Guanina de ARN (G)}$
- $\text{Guanina de ADN (G)} \implies \text{Citosina de ARN (C)}$

La ARN polimerasa realiza la lectura en dirección $3' \to 5'$ y sintetiza en dirección $5' \to 3'$:

$$\begin{array}{rcccccccccccc}
\text{Molde ADN } (3'\to 5'): & \text{T} & \text{A} & \text{C} & \text{G} & \text{G} & \text{A} & \text{C} & \text{T} & \text{T} & \text{A} & \text{C} & \text{G} \\
& \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow & \updownarrow \\
\text{ARNm } (5'\to 3'): & \mathbf{A} & \mathbf{U} & \mathbf{G} & \mathbf{C} & \mathbf{C} & \mathbf{U} & \mathbf{G} & \mathbf{A} & \mathbf{A} & \mathbf{U} & \mathbf{G} & \mathbf{C}
\end{array}$$

Por tanto, el transcrito primario de ARNm resultante es:

$$\mathbf{5'\text{ – AUG CCU GAA UGC – }3'}$$

Agrupado en sus marcos de lectura de tripletes (codones):
- **Codón 1:** $5'\text{-AUG-}3' \implies \text{Metionina (Met / M)}$ — *Codón de Inicio*
- **Codón 2:** $5'\text{-CCU-}3' \implies \text{Prolina (Pro / P)}$
- **Codón 3:** $5'\text{-GAA-}3' \implies \text{Ácido Glutámico (Glu / E)}$
- **Codón 4:** $5'\text{-UGC-}3' \implies \text{Cisteína (Cys / C)}$

---

### 4. Arquitectura Génica: Región Promotora vs. Región Codificante

Un gen funcional no consiste únicamente en la secuencia que codifica un polipéptido; posee una organización modular altamente estructurada:

```text
       [ Río arriba / Upstream ]                      [ Río abajo / Downstream ]
  5' ───[ REGIÓN PROMOTORA ]───┬──────────────[ REGIÓN CODIFICANTE (CDS) ]──────────────┬─── 3'
                               │                                                         │
                        Sitio de Inicio                                            Sitio de Fin
                       Transcripción (+1)                                        Transcripción
```

#### A. Región Promotora
- **Ubicación topológica:** Se encuentra situada *río arriba* (*upstream*, hacia el extremo $5'$ de la hebra codificante y el extremo $3'$ de la hebra molde), adyacente pero previa al sitio de inicio de la transcripción (designado formalmente como la posición nucleotídica $+1$).
- **Características estructurales y motivos conservados:**
  - **En organismos procariotas:** Presenta dos secuencias consenso canónicas: la caja $-35$ ($5'\text{-TTGACA-}3'$) y la caja $-10$ o caja de Pribnow ($5'\text{-TATAAT-}3'$), separadas por una distancia óptima de $16\text{-}18$ pares de bases.
  - **En organismos eucariotas:** El promotor mínimo o basal suele albergar elementos como la caja TATA (alrededor de la posición $-25$ a $-30$, consenso $5'\text{-TATAAA-}3'$), el elemento de reconocimiento TFIIB (BRE) y el iniciador (Inr).
- **Función fisiológica:**
  1. **Reclutamiento de la holoenzima:** Es reconocida específicamente por el factor sigma ($\sigma$) en bacterias o por los factores generales de transcripción (complejo TFIID/TBP) en eucariotas.
  2. **Orientación y polaridad:** Al ser una secuencia asimétrica, orienta espacialmente a la ARN polimerasa, determinando sin ambigüedad cuál de las dos cadenas de la doble hélice funcionará como molde y la dirección unidireccional de avance de la transcripción.
  3. **Apertura de la burbuja transcripcional:** Facilita la transición del "complejo cerrado" al "complejo abierto" mediante la fusión desnaturalizante de pares de bases $\text{A-T}$, permitiendo el acceso de la polimerasa al molde monocatenario sin necesidad de cebador.

> [!NOTE]
> La región promotora **no forma parte de la secuencia codificante** y **no es traducida a proteína**. Cumple una función estrictamente reguladora del flujo de expresión génica.

#### B. Región Codificante (*Coding Sequence* / CDS / Marco Abierto de Lectura - ORF)
- **Ubicación topológica:** Se localiza *río abajo* (*downstream*) del promotor y de la región líder $5'$ no traducida ($5'\text{-UTR}$). Comienza puntualmente en el codón de inicio de traducción ($\text{AUG}$) y se extiende hasta alcanzar un codón de terminación ($\text{UAA}$, $\text{UAG}$ o $\text{UGA}$).
- **Características estructurales:** Está constituida por una serie ininterrumpida y colineal de tripletes de nucleótidos (codones). El código se lee de forma no solapante y sin signos de puntuación ("sin comas"). En eucariotas superiores, la región codificante primaria suele estar fragmentada en módulos de **exones** separados por secuencias no codificantes intercaladas denominadas **intrones**, los cuales son eliminados postranscripcionalmente mediante el proceso de corte y empalme (*splicing*).
- **Función fisiológica:** Dicta directamente la secuencia primaria de aminoácidos del polipéptido resultante a través del apareamiento codón-anticodón en los ribosomas. En el fragmento bajo estudio, la región codificante analizada corresponde a la secuencia $5'\text{– ATG CCT GAA TGC – }3'$.

---

### 5. Extensión Práctica con Biopython (Ejercicio 2)

Para la automatización de la transcripción y el análisis del impacto de la polaridad de las hebras, se diseñó el script modular [`scripts/ejercicio_02_transcripcion.py`](scripts/ejercicio_02_transcripcion.py), el cual procesa archivos en formato estándar **FASTA** almacenados en el directorio `data/`.

#### Lectura FASTA y modelo bioinformático de Biopython
En los repositorios biológicos internacionales (NCBI, Ensembl, DDBJ), las secuencias de ADN se depositan universalmente representadas como la **hebra codificante en orientación $5' \to 3'$**. Bajo este estándar:
- El método `Seq.transcribe()` de Biopython asume que la entrada es la cadena codificante y efectúa el reemplazo estequiométrico de $\text{T} \to \text{U}$:
  ```python
  from Bio import SeqIO

  record = SeqIO.read("data/ejercicio_02_adn.fasta", "fasta")
  mrna_seq = record.seq.transcribe()
  # Resultado: Seq('AUGCCUGAAUGC')
  ```

#### Experimento computacional: Alteración de la orientación de la hebra
En el script se ejecutó un experimento comparativo para verificar qué ocurre si se modifica la orientación o la naturaleza de la hebra suministrada:

1. **Hebra Codificante Canónica ($5' \to 3'$):**
   - Entrada: `5'- ATGCCTGAATGC -3'`
   - ARNm: `5'- AUGCCUGAAUGC -3'`
   - Codones: `['AUG', 'CCU', 'GAA', 'UGC']` $\implies$ Péptido: `Met - Pro - Glu - Cys`

2. **Hebra Reversa-Complementaria / Hebra Antisentido ($5' \to 3'$ de la otra cadena):**
   - Generada mediante `record.seq.reverse_complement()`: `5'- GCATTCAGGCAT -3'`
   - Transcrito: `5'- GCAUUCAGGCAU -3'`
   - Codones: `['GCA', 'UUC', 'AGG', 'CAU']` $\implies$ Péptido: `Ala - Phe - Arg - His`

3. **Inversión Lineal Directa ($3' \to 5'$ leída en sentido contrario):**
   - Generada mediante inversión de índices `record.seq[::-1]`: `5'- CGTAAGTCCGTA -3'`
   - Transcrito: `5'- CGUAAGUCCGUA -3'`
   - Codones: `['CGU', 'AAG', 'UCC', 'GUA']` $\implies$ Péptido: `Arg - Lys - Ser - Val`

#### Conclusión del experimento bioinformático:
La direccionalidad química de los ácidos nucleicos es biológicamente unívoca. Cambiar la orientación o seleccionar la hebra no canónica destruye por completo el marco abierto de lectura (ORF), alterando la secuencia polipeptídica codificada, suprimiendo los sitios de iniciación canónicos o introduciendo codones de parada prematuros.

El script genera de manera autónoma el archivo de salida [`data/ejercicio_02_arnm.fasta`](data/ejercicio_02_arnm.fasta) con el transcrito primario estandarizado.

---

## Ejercicio 3: Traducción del ARNm a proteína

### 1. Objetivo
Aplicar los principios del código genético universal para traducir con precisión una pauta de lectura de ARN mensajero a una secuencia polipeptídica funcional, identificar los tripletes reguladores de inicio y terminación, y evaluar mecanísticamente las consecuencias moleculares, celulares y evolutivas de mutaciones puntuales en el codón de inicio y de mutaciones de pérdida de parada (*stop-loss*).

---

### 2. Planteamiento y Análisis de la Pauta de Lectura

Se proporciona el siguiente transcrito maduro de ARN mensajero:

$$\mathbf{5'\text{ – AUG UAU GCU UAA – }3'}$$

El ARN mensajero es una molécula monocatenaria que se lee de manera polarizada y unidireccional en sentido **$5' \to 3'$** por el complejo ribosómico. La secuencia consta de 12 ribonucleótidos organizados en cuatro tripletes continuos o codones no solapantes:

$$\begin{array}{rcccc}
\text{Posición del Codón:} & \text{Codón 1} & \text{Codón 2} & \text{Codón 3} & \text{Codón 4} \\
\text{Secuencia de Ribonucleótidos } (5'\to 3'): & \mathbf{AUG} & \mathbf{UAU} & \mathbf{GCU} & \mathbf{UAA}
\end{array}$$

#### Identificación del codón de inicio y del codón de paro:

1. **Codón de Inicio:**
   - **Triplete:** $\mathbf{5'\text{-AUG-}3'}$ (Codón 1).
   - **Función molecular:** Establece la fase o marco abierto de lectura (*Open Reading Frame*, ORF). Es reconocido en el sitio P (*peptidil*) del ribosoma por el complejo de preiniciación y el ARNt iniciador especializado:
     - En eucariotas: $\text{Met-tRNA}_i^{\text{Met}}$, incorporando **Metionina (Met / M)**.
     - En bacterias: $\text{fMet-tRNA}_i^{\text{fMet}}$, incorporando **N-formilmetionina (fMet)**.
2. **Codón de Paro / Terminación:**
   - **Triplete:** $\mathbf{5'\text{-UAA-}3'}$ (Codón 4, históricamente denominado *Ochre*).
   - **Función molecular:** Señaliza el cese absoluto de la elongación polipeptídica. A diferencia de los codones con sentido, los codones de parada ($\text{UAA}$, $\text{UAG}$ y $\text{UGA}$) **no corresponden a ningún aminoacil-ARNt**. En su lugar, cuando este triplete se ubica en el sitio A (*aminoacil*) del ribosoma, es reconocido por factores proteicos de liberación (*Release Factors*: RF1/RF2 en bacterias; eRF1 en eucariotas). Estos factores inducen una actividad peptidil-transferasa hidrolítica que escinde el enlace éster entre la cadena polipeptídica y el último ARNt, liberando la proteína naciente y promoviendo el reciclaje y disociación de las subunidades ribosómicas.

---

### 3. Traducción a Secuencia Polipeptídica

Consultando la tabla del código genético estándar universal:

| Codón ($5'\to 3'$) | Aminoácido codificado | Código de 3 letras | Código de 1 letra | Propiedades fisicoquímicas del residuo |
| :--- | :--- | :--- | :--- | :--- |
| **AUG** | Metionina | $\text{Met}$ | $\text{M}$ | Apolar alifático, contiene azufre en puente tioéter; aminoácido iniciador canónico. |
| **UAU** | Tirosina | $\text{Tyr}$ | $\text{Y}$ | Aromático, polar no cargado; posee un grupo hidroxilo fenólico ($-\text{OH}$) reactivo. |
| **GCU** | Alanina | $\text{Ala}$ | $\text{A}$ | Apolar alifático, cadena lateral metilo ($-\text{CH}_3$); conformacionalmente versátil. |
| **UAA** | *Codón de Parada (Stop)* | $\text{Stop}$ | $*$ | Señal de disociación y terminación de la traducción; no aporta aminoácido al péptido. |

Por consiguiente, la cadena polipeptídica madura sintetizada a partir de este transcrito es un oligopéptido de 3 aminoácidos (tripéptido):

$$\mathbf{\text{H}_2\text{N – Met – Tyr – Ala – COOH}} \quad \text{o bien en notación estándar: } \mathbf{\text{Met-Tyr-Ala}} \quad (\mathbf{\text{MYA}})$$

---

### 4. Reflexión Biológica: Impacto Molecular de Mutaciones Críticas

#### A. ¿Qué ocurriría si el codón de inicio mutara de $\text{AUG}$ a $\text{GUG}$?

El codón canónico de inicio es $\text{AUG}$. Si se produce una transición $\text{A} \to \text{G}$ en la primera posición ($5'\text{-GUG-}3'$):

1. **En fase de elongación:**
   En el código genético universal estándar, $\text{GUG}$ codifica canónicamente para **Valina (Val / V)** cuando se encuentra en el cuerpo del marco de lectura.

2. **Impacto en la fase de iniciación según el dominio biológico:**
   - **En células eucariotas:**
     El complejo de preiniciación $43\text{S}$ realiza el escaneo unidireccional de la región $5'\text{-UTR}$ (*ribosome scanning*) reconociendo el triplete $\text{AUG}$ mediante el apareamiento exacto de bases de Watson y Crick con el bucle anticodón ($3'\text{-UAC-}5'$) del $\text{Met-tRNA}_i^{\text{Met}}$, estabilizado por la secuencia consenso de Kozak ($(\text{gcc})\text{gccRccAUGG}$). La sustitución por $\text{GUG}$ produce un desapareamiento estérico y termodinámico en la primera base que **reduce drásticamente la eficiencia de iniciación en más de un $95\%$**. Las consecuencias fisiológicas son:
     - **Ausencia de síntesis proteica:** Si no existen otros codones de inicio, el ARNm no se traduce eficazmente (alelo nulo o pérdida total de función proteica). En humanos, mutaciones en el $\text{AUG}$ de inicio son causantes de enfermedades genéticas severas, como ciertas formas graves de talasemia $\beta$ (mutaciones en el gen *HBB*).
     - **Iniciación aberrante alternativa *downstream*:** El ribosoma puede continuar el escaneo hasta encontrar un codón $\text{AUG}$ críptico o secundario más adelante. Si dicho $\text{AUG}$ está fuera del marco de lectura original (*out-of-frame*), se sintetizará un polipéptido aberrante, afuncional o prematuramente truncado.
   - **En bacterias y arqueas:**
     El inicio no depende de un escaneo desde el extremo $5'$, sino de la interacción directa entre la secuencia de **Shine-Dalgarno** del ARNm y el ARN ribosómico $16\text{S}$. En procariotas, $\text{GUG}$ actúa como un **codón de inicio alternativo funcional** (presente de forma natural en $\sim 8\text{-}14\%$ de los genes bacterianos, como en el operón *lac* de *E. coli*). Cuando el complejo se ensambla con la ayuda del factor IF-2 y el $\text{fMet-tRNA}_i^{\text{fMet}}$, se incorpora igualmente **N-formilmetionina (fMet)** como primer residuo, pero con una afinidad y tasa de iniciación significativamente menores ($\sim 10\text{-}30\%$ respecto a $\text{AUG}$).

---

#### B. ¿Qué ocurriría si el codón de paro desapareciera por mutación?

La conversión de un codón de terminación en un codón con sentido debido a una mutación puntual (por ejemplo, $\text{UAA} \to \text{CAA}$ [Glutamina], $\text{UAA} \to \text{AAA}$ [Lisina] o $\text{UAA} \to \text{UAC}$ [Tirosina]) se clasifica como una **mutación de pérdida de parada (*Stop-loss* o *Non-stop mutation*)**.

Las consecuencias moleculares y celulares son de gran relevancia fisiopatológica:

```text
Transición normal:
  ARNm: ... [GCU] [UAA] ─── (Fin de síntesis: liberación normal del péptido Met-Tyr-Ala)
              Ala   STOP

Mutación Stop-Loss (UAA -> CAA):
  ARNm: ... [GCU] [CAA] [GCC] [UGA] ... (Cola poli-A)
              Ala   Gln   Ala   STOP críptico
              └─────────────────────┘
              Extensión C-terminal anómala del polipéptido
```

1. **Elongación aberrante y lectura continua en la región $3'\text{-UTR}$ (*Translational Readthrough*):**
   Al carecer de señal de terminación en el sitio A, el ribosoma transloca y continúa incorporando aminoácidos leyendo a través de la región no traducida $3'\text{-UTR}$ (*untranslated region*), la cual carece de selección evolutiva para codificar polipéptidos funcionales.
2. **Síntesis de colas de polilisina y estancamiento ribosómico (*Ribosome Stalling*):**
   Si no se encuentra un codón de parada críptico en fase dentro de la región $3'\text{-UTR}$, el ribosoma alcanza la cola de poliadenilación ($\text{poli(A)}$). Al traducir los codones $\text{AAA}$, incorpora tramos continuos de **polilisina básica** ($[\text{Lys}]_n$). Las cargas electrostáticas fuertemente positivas de la polilisina interaccionan con las paredes electrostáticamente negativas del túnel de salida ribosómico (conformadas por ARNr), provocando el **estancamiento físico irreversible de la maquinaria ribosómica**.
3. **Disparo de los sistemas de vigilancia celular y control de calidad:**
   - **Vía *Non-stop Decay* (NSD) en eucariotas:** Las subunidades ribosómicas estancadas son detectadas por el complejo proteico de rescate **Pelota/Hbs1**. El ARNm es canalizado hacia el complejo ribonucleolítico del **exosoma** para su degradación acelerada en dirección $3' \to 5'$. A su vez, el polipéptido extendido anómalo es marcado covalentemente con ubiquitina por la ligasa **Ltn1** (complejo RQC) y degradado en el proteasoma $26\text{S}$.
   - **Sistema del ARNtm (*tmRNA* / SsrA) en procariotas:** En bacterias, el ARN de transferencia-mensajero rescata al ribosoma detenido en el extremo $3'$ de un ARNm truncado o sin parada, añadiendo una etiqueta peptídica terminal ($\text{AANDENYALAA}$) que actúa como diana específica para las proteasas citosólicas ClpXP y Lon.
4. **Impacto biológico y patológico:**
   Cualquier proteína que logre eludir la degradación presentará una extensión C-terminal aberrante que alterará su superficie hidrofóbica, desestabilizará su plegamiento terciario o cuaternario, causará agregación tóxica intracelular y conllevará la pérdida total de actividad biológica.

---

### 5. Extensión Práctica con Biopython (Ejercicio 3)

Se implementó el script ejecutable [`scripts/ejercicio_03_traduccion.py`](scripts/ejercicio_03_traduccion.py) para automatizar la traducción del ARNm utilizando el módulo `Bio.Seq.Seq`.

#### Métodos de traducción computacional: `translate()` vs. `to_stop=True`
- **Traducción integral (`seq.translate()`):** Biopython procesa la secuencia en tripletes y traduce cada codón canónico según la tabla estándar del código genético (`NCBI Table 1`). Los codones de parada se representan universalmente mediante el carácter de asterisco (`*`):
  ```python
  from Bio.Seq import Seq

  mrna = Seq("AUGUAUGCUUAA")
  prot_full = mrna.translate()
  # Resultado: Seq('MYA*')
  ```
- **Traducción del péptido maduro (`seq.translate(to_stop=True)`):** Interrumpe la elongación inmediatamente antes de traducir el primer codón de parada, simulando la acción biológica de los factores de liberación ribosómicos:
  ```python
  peptide = mrna.translate(to_stop=True)
  # Resultado: Seq('MYA')
  ```

#### Simulación computacional de mutaciones
El script simula formalmente ambas condiciones mutagénicas analizadas teóricamente:
1. **Mutación en el inicio ($\text{AUG} \to \text{GUG}$):**
   ```python
   mut_inicio = Seq("GUGUAUGCUUAA").translate(to_stop=True)
   # Resultado: Seq('VYA')
   ```
   Demuestra cómo el código genético estándar decodifica $\text{GUG}$ como Valina ($\text{V}$) en lugar de Metionina ($\text{M}$).
2. **Mutación de pérdida de parada (*Stop-loss* $\text{UAA} \to \text{CAA}$):**
   ```python
   # Simulación con lectura hacia la región 3'-UTR hipotética
   mut_stoploss = Seq("AUGUAUGCUCAAGCCUGA").translate(to_stop=True)
   # Resultado: Seq('MYAQA')
   ```
   Evidencia la elongación aberrante con incorporación de Glutamina ($\text{CAA} \to \text{Q}$) y Alanina ($\text{GCC} \to \text{A}$) hasta alcanzar un codón de parada secundario ($\text{UGA}$).

Las aserciones algorítmicas (`assert`) implementadas en el script confirman una concordancia del $100\%$ entre la predicción manual y la ejecución bioinformática.

---

## Ejercicio 4: Splicing alternativo

### 1. Objetivo
Comprender los mecanismos moleculares del procesamiento postranscripcional de maduración del pre-ARNm en eucariotas, analizar la acción catalítica del espliceosoma y la arquitectura modular de los exones, diseñar combinaciones viables de corte y empalme alternativo (*alternative splicing*) y evaluar cómo este mecanismo expande exponencialmente la diversidad del proteoma sin requerir un incremento proporcional en el número de genes codificantes, ilustrándolo con el caso paradigmático de Ensembl del receptor humano **FGFR2** (`ENSG00000066468`).

---

### 2. Planteamiento y Diseño de Combinaciones de Splicing (Gen de 5 Exones)

Se considera un gen eucariota hipotético estructurado en **5 exones** intercalados por intrones no codificantes en su transcrito primario (pre-ARNm):

```text
Transcrito primario (pre-ARNm):
  5'-[CAP]───[Exón 1]───(Intrón 1)───[Exón 2]───(Intrón 2)───[Exón 3]───(Intrón 3)───[Exón 4]───(Intrón 4)───[Exón 5]───[poli-A]-3'
```

El **espliceosoma** (complejo ribonucleoproteico nuclear dinámico constituido por las partículas snRNP $\text{U1}$, $\text{U2}$, $\text{U4/U6}$ y $\text{U5}$, junto con más de 100 factores proteicos accesorios) reconoce de manera altamente regulada las secuencias consenso dadoras de splicing en $5'$ ($\text{GU}$), aceptoras en $3'$ ($\text{AG}$) y el punto de ramificación adenina ($\text{A}$).

Mediante la modulación del reconocimiento de estos sitios por factores reguladores (*e.g.*, proteínas promotoras SR y represoras hnRNP), se pueden ensamblar diferentes combinaciones de exones en el ARNm maduro. Se diseñan las siguientes isoformas:

```text
1. Isoforma Canónica / Constitutiva (1-2-3-4-5):
   5'-[Exón 1]───[Exón 2]───[Exón 3]───[Exón 4]───[Exón 5]-3'

2. Isoforma Alternativa A (Salto de Exón / Exon Skipping: 1-2-4-5):
   5'-[Exón 1]───[Exón 2]──────────────[Exón 4]───[Exón 5]-3'   (Exclusión del Exón 3)

3. Isoforma Alternativa B (Exclusión Doble / Modular: 1-3-5):
   5'-[Exón 1]──────────────[Exón 3]──────────────[Exón 5]-3'   (Exclusión de Exones 2 y 4)
```

---

### 3. Diferencias Estructurales y Funcionales en las Proteínas Resultantes

En las proteínas eucariotas, los exones suelen codificar **dominios estructurales o módulos funcionales discretos** (modelo de *exon shuffling*). Por tanto, la inclusión o exclusión de exones específicos genera diferencias cuantitativas y cualitativas cruciales:

1. **Alteración en la longitud y peso molecular:**
   - La **Isoforma Constitutiva (1-2-3-4-5)** genera la proteína de máxima longitud ($23\text{ aa}$ en nuestro modelo), conservando la integridad de todos los dominios codificados.
   - La **Isoforma 1-2-4-5** produce un polipéptido acortado ($18\text{ aa}$) por la eliminación selectiva de los aminoácidos del Exón 3.
   - La **Isoforma 1-3-5** produce una variante compacta ($14\text{ aa}$) con un tamaño molecular sustancialmente inferior.

2. **Ganancia, pérdida o modulación de dominios funcionales:**
   - **Exón 1 (Péptido señal e iniciación):** Conservado en todas las isoformas, garantiza la translocación hacia el retículo endoplásmico o el inicio correcto de la síntesis citosólica.
   - **Exón 3 (Dominio catalítico / Sitio activo / Unión a ligando):**
     - En la **Isoforma 1-2-3-4-5**, el dominio catalítico o de interacción está presente y enmarcado en su contexto terciario natural.
     - En la **Isoforma 1-2-4-5**, la pérdida del Exón 3 genera una proteína que carece del sitio de unión o actividad enzimática. Con frecuencia, estas variantes actúan biológicamente como **isoformas dominantes-negativas**: conservan la capacidad de dimerizar o anclarse a la membrana pero son catalíticamente inertes, compitiendo con la isoforma completa y regulando negativamente la vía de señalización.
   - **Exones 2 y 4 (Módulos espaciadores o de regulación alostérica):**
     - En la **Isoforma 1-3-5**, la eliminación simultánea de los módulos conectores 2 y 4 genera un receptor o enzima rígida, donde el dominio activo se asocia directamente a la región C-terminal o transmembrana (Exón 5), pudiendo provocar activación constitutiva independiente de efectores o cambios drásticos en la estabilidad proteica.

3. **Mantenimiento del marco abierto de lectura (ORF) vs. Salto de fase:**
   - Si los exones excluidos poseen una longitud nucleotídica múltiplo exacto de 3 (como en nuestro diseño de 12, 15 y 18 nt), la lectura downstream en los exones 4 y 5 permanece en fase (*in-frame*), dando lugar a isoformas proteicas viables.
   - Si un exón alternativo posee un número de nucleótidos no divisible por 3, su exclusión o inclusión altera el marco de lectura (*frameshift*), introduciendo un codón de parada prematuro (**PTC**, *Premature Termination Codon*) en los exones siguientes. Dicho transcrito es detectado y degradado por la vía de **degradación mediada por mutaciones sin sentido (*Nonsense-Mediated mRNA Decay*, NMD)**, funcionando el splicing alternativo en este caso como un interruptor molecular que apaga la expresión génica a nivel postranscripcional (mecanismo *AS-NMD*).

---

### 4. Reflexión Biológica: Diversidad del Proteoma sin Aumento de Genes

#### La "Paradoja del Número de Genes" y la Complejidad Biológica
Durante el Proyecto Genoma Humano, se estimaba inicialmente que la complejidad fisiológica, cognitiva e inmunológica humana requeriría entre $100.000$ y $150.000$ genes codificantes. El resultado real reveló que el genoma humano posee únicamente alrededor de **$\sim 20.000$ genes codificantes de proteínas**, una cifra prácticamente equivalente a la del nematodo milimétrico *Caenorhabditis elegans* ($\sim 19.000\text{ genes}$) y notablemente inferior a la de plantas como el arroz (*Oryza sativa*, $>30.000\text{ genes}$).

¿Cómo se explica esta discrepancia biológica? La respuesta radica fundamentalmente en el **splicing alternativo**:

```text
             Genoma Humano (~20.000 genes)
                          │
                          ▼ (Transcripción masiva)
             Pre-ARNm no procesados
                          │
                          ▼ (Splicing alternativo en >95% de genes multiexónicos)
             Transcriptoma (>100.000 - 150.000 transcritos de ARNm maduros)
                          │
                          ▼ (Traducción + Modificaciones Postraduccionales)
             Proteoma (>500.000 - 1.000.000 de proteoformas funcionales)
```

#### Ventajas evolutivas y bioenergéticas del mecanismo:
1. **Economía de espacio genómico y bajo coste energético:**
   Mantener, replicar y reparar covalentemente cromosomas con cientos de miles de genes individuales exigiría un consumo bioenergético metabólico descomunal, retrasaría la cinética de división celular y multiplicaría el blanco mutacional genómico. El splicing alternativo comprime la información: un único locus genómico actúa como una matriz multifuncional.
2. **Modularidad combinatoria (*Exon Shuffling*):**
   Los exones codifican motivos estructurales discretos (dominios SH2, dedos de zinc, dominios quinasa, etc.). Al combinar estos módulos como bloques de construcción, un solo gen puede codificar variantes solubles, unidas a membrana, de alta afinidad o de baja afinidad.
3. **Regulación espacio-temporal ultraprecisa:**
   Permite expresar isoformas especializadas según el tipo de tejido (isoformas cardíacas vs. neuronales vs. hepáticas) o según la etapa del desarrollo embrionario, orquestadas por concentraciones relativas de factores de corte y empalme tisulares sin necesidad de activar promotores independientes.

---

### 5. Extensión con Biopython y Bases de Datos (Ensembl: Caso *FGFR2*)

Para aterrizar el modelado teórico en la genómica funcional real, se analizó el locus del gen humano **FGFR2** (*Fibroblast Growth Factor Receptor 2*, Ensembl ID: **`ENSG00000066468`**, localizado en el cromosoma 10: 121,478,332-121,598,458) mediante el script [`scripts/ejercicio_04_splicing.py`](scripts/ejercicio_04_splicing.py) y las secuencias extraídas en [`data/fgfr2_isoformas.fasta`](data/fgfr2_isoformas.fasta).

#### Caso paradigmático: Exones mutuamente excluyentes IIIb y IIIc
El pre-ARNm de *FGFR2* sufre un evento clásico de **splicing alternativo mutuamente excluyente** en el dominio extracelular:

```text
Locus FGFR2 pre-ARNm:
  ───[Exón 7 (IgIII N-term)]───┬───[Exón 8 / IIIb]───┬───[Exón 9 / IIIc]───┬───[Exón 10 (Transmembrana)]───
                               │                     │                     │
      Splicing en Epitelio     └─────────────────────┴──────(Excluido)─────┘ ───> Isoforma FGFR2-IIIb (ENST00000358487)
                               │                     │                     │
      Splicing en Mesénquima   └─────(Excluido)──────┴─────────────────────┘ ───> Isoforma FGFR2-IIIc (ENST00000356226)
```

Ambos exones codifican la mitad C-terminal del tercer bucle similar a inmunoglobulina (dominio $\text{IgIII}$), el cual delimita físicamente la cavidad de reconocimiento y unión al ligando FGF.

#### Comparación de transcritos e isoformas de Ensembl:

| Propiedad | Isoforma Epitelial: **FGFR2-IIIb** | Isoforma Mesenquimal: **FGFR2-IIIc** |
| :--- | :--- | :--- |
| **Transcrito Ensembl** | `ENST00000358487` (FGFR2-201) | `ENST00000356226` (FGFR2-202) |
| **Exón variable incorporado** | Exón 8 (IIIb) | Exón 9 (IIIc) |
| **Expresión tisular** | Estrictamente células **epiteliales** | Estrictamente células **mesenquimales** |
| **Péptido del dominio variable** | `IESSNKINLNVSFNNVTWLEDAGNYTCLAGNSIGISFHSAWLTVL` | `FKCPSSGTPNPTLRWLKNGKEFKPDHRIGGYKVRYATWSIIMDSV` |
| **Afinidad de Ligando** | Alta afinidad por **FGF7 (KGF)** y **FGF10** | Alta afinidad por **FGF2 (bFGF)**, **FGF4**, **FGF6** y **FGF9** |
| **Función fisiológica** | Comunicación paracrina unidireccional epitelio-mesénquima en la piel y el epitelio pulmonar. | Morfogénesis mesenquimal, desarrollo esquelético y osteogénesis. |

#### Relevancia médica y Transición Epitelio-Mesénquima (EMT)
La conmutación o intercambio de splicing alternativo de la isoforma **FGFR2-IIIb** a la **FGFR2-IIIc** es un marcador fisiopatológico capital de la **Transición Epitelio-Mesénquima (EMT)**. 
- Durante la progresión del cáncer epitelial (carcinomas), factores de transcripción como SNAIL, SLUG y TWIST silencian las proteínas de unión al ARN reguladoras epiteliales (como ESRP1 y ESRP2), provocando el salto a la isoforma mesenquimal IIIc.
- La expresión anómala de FGFR2-IIIc en células tumorales permite que respondan autocrinamente al ligando FGF2, desencadenando la pérdida de polaridad apical-basal, la disolución de uniones intercelulares de cadherina y una agresiva capacidad de invasión celular y metástasis a distancia.

El script `scripts/ejercicio_04_splicing.py` reproduce computacionalmente tanto el modelado teórico de las combinaciones de exones como el análisis secuencial comparativo de los dos transcritos de Ensembl.


