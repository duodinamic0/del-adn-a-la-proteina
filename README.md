# Del ADN a la Proteína: Dogma Central y Análisis Bioinformático

---

## Ejercicio 1. Replicación del ADN

**1. Dúplex parental y mecanismo semiconservativo:** 
A partir del dúplex parental `5'-ATG CCG TTA GCT-3'` / `3'-TAC GGC AAT CGA-5'`, la replicación semiconservativa (Meselson y Stahl) separa ambas hebras para sintetizar cadenas complementarias antiparalelas, rindiendo dos moléculas hijas bicatenarias idénticas:
- *Hija 1:* Molde `5'-ATG CCG TTA GCT-3'` + hebra neosintetizada `3'-TAC GGC AAT CGA-5'` (orientación estándar: `5'-AGC TAA CGG CAT-3'`).
- *Hija 2:* Molde `3'-TAC GGC AAT CGA-5'` + hebra neosintetizada `5'-ATG CCG TTA GCT-3'`.

**2. Maquinaria enzimática:** 
La *helicasa* desenrolla la doble hélice rompiendo puentes de hidrógeno con gasto de ATP; la *primasa* sintetiza cebadores de ARN aportando el extremo $3'\text{-OH}$ libre; la *ADN polimerasa* elonga en sentido $5' \to 3'$ y corrige errores mediante exonucleasa $3' \to 5'$ (*proofreading*); y la *ligasa* sella los enlaces fosfodiéster en las muescas y fragmentos de Okazaki.

**3. Reflexión (Error no corregido):** 
Si un desapareamiento (*mismatch*) elude el *proofreading* y el sistema MMR, en la siguiente ronda de replicación la base anómala actúa como molde legítimo fijando covalentemente una mutación permanente. En regiones codificantes puede ser silenciosa, de cambio de sentido (*missense*) o sin sentido (*nonsense*, codón de paro prematuro), con riesgo oncogénico en células somáticas o anomalías hereditarias en línea germinal.

**4. Extensión con Biopython:** 
En el Jupyter Notebook `notebooks/ejercicio_01_replicacion.ipynb` se valida la secuencia con `Bio.Seq`: `dna.complement()` rinde la orientación `3'->5'` (`TACGGCAATCGA`) y `dna.reverse_complement()` la estándar `5'->3'` (`AGCTAACGGCAT`), coincidiendo al 100% con la deducción manual.

---

## Ejercicio 2. Transcripción del ADN a ARN

**1. Discriminación de cadenas y transcrito:** 
Dado que la ARN polimerasa sintetiza en sentido $5' \to 3'$, utiliza como molde antiparalelo la **hebra inferior** (`3'-TAC GGA CTT ACG-5'`), cuyo triplete `3'-TAC-5'` genera el codón canónico de inicio. La hebra superior (`5'-ATG CCT GAA TGC-3'`) actúa como cadena codificante. El transcrito de ARNm obtenido es:
$$\text{5'-AUG CCU GAA UGC-3' \quad (Met - Pro - Glu - Cys)}$$

**2. Arquitectura génica:** 
El *promotor* es una región reguladora no traducida situada río arriba (*upstream*, flanco 5') con motivos consenso (caja TATA o $-10$/$-35$) que recluta a la ARN polimerasa y fija el inicio de la transcripción (+1). La *región codificante (CDS/ORF)* se extiende río abajo desde el codón de inicio (AUG) hasta el de parada, dictando la secuencia colineal de aminoácidos.

**3. Extensión con Biopython:** 
En el Jupyter Notebook `notebooks/ejercicio_02_transcripcion.ipynb` se procesa el FASTA con `record.seq.transcribe()`. Al invertir la secuencia (`seq[::-1]`) o transcribir la hebra complementaria reversa, se destruye el marco de lectura dando péptidos totalmente diferentes (`Arg-Lys-Ser-Val` o `Ala-Phe-Arg-His`), demostrando que la información genética depende de forma unívoca de la polaridad química $5' \to 3'$.

---

## Ejercicio 3. Traducción del ARNm a proteína

**1. Pauta de lectura y péptido resultante:** 
En el transcrito `5'-AUG UAU GCU UAA-3'`, el codón `AUG` fija el inicio incorporando Metionina, seguido de Tirosina (`UAU`) y Alanina (`GCU`). El codón de terminación `UAA` no codifica aminoácidos sino que recluta factores proteicos de liberación (RF/eRF1) para hidrolizar el polipéptido maduro: $\mathbf{\text{Met-Tyr-Ala}}$ (**MYA**).

**2. Reflexión (Impacto de mutaciones críticas):** 
- *Mutación de inicio (**AUG** => **GUG**):* En eucariotas impide el reconocimiento por el complejo de preiniciación 43S, suprimiendo la traducción ($>90\%$) o forzando inicios aberrantes secundarios fuera de marco. En bacterias actúa como inicio alternativo funcional guiado por la secuencia Shine-Dalgarno con menor afinidad.
- *Pérdida de parada (*Stop-loss*, ej. $\text{UAA} \to \text{CAA}$):* El ribosoma invade la región $3'\text{-UTR}$ (*readthrough*) y al llegar a la cola poli-A traduce polilisinas básicas que colapsan el túnel ribosómico (*ribosome stalling*), activando la degradación del ARNm por Non-stop Decay (NSD) y proteólisis de la proteína en el proteasoma 26S (RQC/Ltn1).

**3. Extensión con Biopython:** 
En el Jupyter Notebook `notebooks/ejercicio_03_traduccion.ipynb`, `mrna.translate(to_stop=True)` valida la síntesis de `MYA`. Se simula la traducción aberrante por mutación en el inicio (`VYA`) y la elongación anómala en *stop-loss* (`MYAQA...`), confirmando las predicciones teóricas.

---

## Ejercicio 4. Splicing alternativo

**1. Modelado en un gen de 5 exones:** 
A partir del pre-ARNm con exones 1-2-3-4-5, el espliceosoma (snRNP U1-U6) genera diversas combinaciones maduras: la *Canónica* (1-2-3-4-5, longitud completa), la *Isoforma A* (salto del Exón 3, 1-2-4-5) y la *Isoforma B* (exclusión modular de 2 y 4, 1-3-5).

**2. Consecuencias funcionales:** 
Los exones representan dominios discretos (*exon shuffling*). La pérdida del Exón 3 genera proteínas sin dominio activo que pueden actuar como **dominantes-negativas** (compitiendo por receptores sin señalizar). Si la exclusión preserva un múltiplo de 3 nt el marco permanece *in-frame*; si altera la pauta (*frameshift*) genera codones de parada prematuros (PTC), induciendo degradación por la vía **Nonsense-Mediated Decay (NMD)**.

**3. Diversidad proteica y caso FGFR2:** 
Resuelve la paradoja del genoma humano ($\sim 20.000$ genes generan $>100.000$ transcritos): optimiza el espacio cromosómico y el coste bioenergético celular coordinando perfiles tisulares específicos. En el Jupyter Notebook `notebooks/ejercicio_04_splicing.ipynb` se analiza el gen humano **FGFR2** (`ENSG00000066468`), donde la inclusión mutuamente excluyente del Exón 8 (isoforma epitelial IIIb, afín a FGF7/10) o el Exón 9 (isoforma mesenquimal IIIc, afín a FGF2) gobierna la identidad celular. El cambio patológico de IIIb a IIIc desencadena la **Transición Epitelio-Mesénquima (EMT)** e invasión en carcinomas.

---

## Ejercicio 5. Introducción a las proteínas

**1. Polaridad química:** 
En el péptido $\text{Met-Ile-Ser-Gly-Val-Lys-His}$ (`MISGVKH`), el extremo **N-terminal** ($\text{H}_3\text{N}^+-$) corresponde a la Metionina (grupo amino libre) y el **C-terminal** ($-\text{COO}^-$) a la Histidina (grupo carboxilo libre).

**2. Reflexión biofísica (Plegamiento y mutaciones):** 
Según el Principio de Anfinsen, la estructura nativa está codificada en la secuencia lineal. La alternancia de residuos apolares y polares define motivos regulares como hélices $\alpha$ (anfipáticas cada 3.6 residuos) y láminas $\beta$ (cada 2 residuos). El plegamiento está impulsado por el **efecto hidrofóbico**, empaquetando cadenas apolares en un núcleo interno anhidro. Sustituir un residuo apolar interno por uno cargado/hidrofílico (ej. $\text{Ile} \to \text{Glu}$) acarrea una penalización de desolvatación masiva ($+15 \text{ a } +20\text{ kcal/mol}$) que supera la tenue estabilidad nativa ($\Delta G \approx -5 \text{ a } -15\text{ kcal/mol}$), provocando el colapso conformacional, agregación amiloide citotóxica y degradación proteasomal.

**3. Extensión con el PDB (Lisozima `1AKI`):** 
En el Jupyter Notebook `notebooks/ejercicio_05_proteinas.ipynb` se analiza la Lisozima (129 aa, 44.2% hélices $\alpha$, 4.7% láminas $\beta$). Los residuos hidrofóbicos se concentran a una distancia media de $10.89\text{ AA}$ del centro de masa, frente a $14.83\text{ AA}$ de los hidrofílicos cargados ($\Delta r = +3.94\text{ AA}$ hacia la superficie solvente), evidenciando el núcleo apolar. La mutación simulada `Ile98Glu` (enterrada a 9.36 Å del centro) desestabilizaría letalmente la arquitectura terciaria nativa.

---

## Ejercicio 6. Actividad integradora: del ADN a la proteína

**1. Secuencia de referencia (NCBI GenBank):** 
Se procesó la CDS de la **Insulina humana (*INS*)** ([`NM_000207.3`](https://www.ncbi.nlm.nih.gov/nuccore/NM_000207.3), 333 pb / 111 codones) en `data/insulina_humana_cds.fasta`, la cual codifica la preproinsulina de 110 aa ([`NP_000198.1`](https://www.ncbi.nlm.nih.gov/protein/NP_000198.1)).

**2. Pipeline del Dogma Central:** 
- *Replicación:* La helicasa abre el dúplex parental y la polimerasa sintetiza en $5' \to 3'$, originando dos moléculas hijas bicatenarias idénticas.
- *Transcripción:* La ARN polimerasa II transcribe la hebra molde ($3' \to 5'$) en el ARNm maduro de 333 nt (`5'-AUG GCC CUG... AAC UAG-3'`), guardado en `data/ejercicio_06_insulina_arnm.fasta`.
- *Traducción:* El ribosoma genera la preproinsulina (110 aa: Péptido señal [1-24], Cadena B [25-54], Péptido C conector [55-89] escindido por convertasas PC1/3-PC2, y Cadena A [90-110] unida por 3 puentes disulfuro), guardada en `data/ejercicio_06_insulina_proteina.fasta`.

**3. Reflexión (Punto más vulnerable a errores):** 
Aunque cuantitativamente la traducción comete más errores ($\sim 10^{-3}\text{-}10^{-4}$ vs. $10^{-9}\text{-}10^{-10}$ en replicación), **la replicación del ADN es el eslabón más vulnerable del dogma central**: los errores en ARNm o proteínas son transitorios y diluidos por el recambio celular, mientras que un fallo replicativo se fija covalentemente como mutación permanente, transmitiéndose al 100% de transcritos y proteínas hijas, provocando transformación oncogénica en células somáticas o enfermedades genéticas incurables en línea germinal.

**4. Validación en Biopython:** 
El pipeline en el Jupyter Notebook `notebooks/ejercicio_06_pipeline_dogma.ipynb` reporta en tiempo real cada fase con `Bio.Seq`/`Bio.SeqIO` y valida la secuencia al 100% mediante `assert` frente a la referencia oficial del NCBI.
