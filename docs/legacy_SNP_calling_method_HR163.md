# HR163 (2502G0098 S / 2502G0099 R) - come sono stati chiamati i 54 SNP

## Tool usato: BLASTN assembly-vs-assembly. Nessun read mapping, nessun Snippy.

Catena esatta (frame 7da28d36, env `blast-dpcr`, cell 281-290):

1. Estrazione del cromosoma (contig `1`) dai due assembly ibridi:
   0098 = 5,303,302 bp; 0099 = 5,304,421 bp.
2. `makeblastdb -in chr_2502G0098.fasta -dbtype nucl -out dbchr98`
   -> il reference e' l'isolato SENSIBILE; le coordinate finali sono sul RESISTENTE (query).
3. Finestre di 5 kb sul cromosoma 0099 (1,061 finestre), ciascuna allineata al db:
   `blastn -query win99.fasta -db dbchr98 -task megablast -evalue 1e-100 -max_target_seqs 1`
   best HSP per finestra -> 1,061/1,061 finestre con hit, 5,300,821 bp allineati,
   113 mismatch, 29 gap, divergenza 0.002132%.
4. Le sole finestre con >=1 mismatch ri-allineate con output pairwise
   (`blastn -outfmt 0`) e i mismatch estratti dalle righe Query/Sbjct -> `snps_99_vs_98.tsv`
   (113 posizioni).
5. Tabella + flag array di ripetizioni in tandem (3,482,000-3,485,000):
   59 posizioni dentro l'array, **54 fuori** -> il numero riportato nella lettera.
   `HR163_chromosomal_SNPs.csv`
6. SNP folA promotore: **pos 4,489,835 sul cromosoma 0099, G(S) -> A(R)**, fuori dall'array.

Cosa NON e' stato fatto: mapping delle read, filtri su VAF/profondita', masking di profagi
(PhiSpy), masking di ricombinazione (Gubbins), core-genome alignment.
Il "54" e' quindi un conteggio di differenze fra consensus, non un set di varianti
validate sulle read.

## Precedente read-based nel progetto (coppia 137, NON HR163)
Frame 95e740da, env `bactgenomics`: MUMmer **dnadiff** sugli assembly (29 SNP) + verifica
sulle read per posizione -> `dnadiff_snp_read_reconciliation.csv`. Solo **3/29 confermati**
(VAF_R ~0.94-0.99 con VAF_S = 0); le altre 26 avevano VAF simile nei due bracci (0.08-0.73
in entrambi), cioe' rumore condiviso / multi-mapping, non differenze fissate.
Precedente utile: su una coppia clonale il passaggio dall'assembly alle read ha tagliato
il conteggio di ~10x.

## Confronto con la pipeline Posteraro 2024 (Snippy + PhiSpy + Gubbins + SNP-sites)

| | nostra | Posteraro 2024 |
|---|---|---|
| input | 2 assembly ibridi | read vs 1 reference ibrido Unicycler |
| allineamento | blastn megablast, finestre 5 kb | bwa mem (dentro Snippy) |
| chiamata | mismatch nell'HSP | freebayes con filtri Snippy |
| soglie | nessuna | default: mincov 10, minfrac 0.9, minqual 100 |
| masking | solo array tandem, manuale | profagi (PhiSpy) + ricombinazione (Gubbins), via snippy-core --mask |
| output | catalogo pairwise di differenze | core-genome SNP alignment (SNP-sites) per MST |
| scopo | quali basi differiscono in questa coppia | distanze fra molti isolati |

Punti da controllare quando arrivano gli output Snippy:

1. **Reference dichiarato.** Snippy ha un solo reference. La nostra tabella e' in coordinate
   0099 con 0098 come reference: se usi 0099 come reference i segni si invertono e le
   posizioni cambiano di ~1 kb (i due cromosomi differiscono di 1,119 bp in lunghezza).
2. **`minfrac 0.9`.** In un ceppo eteroresistente una variante sottoclonale viene scartata
   di default. Se il conteggio crolla, riesegui con `--minfrac 0.5` e guarda le VAF in
   `snps.vcf` (campo AO/RO) prima di concludere che una posizione non esiste.
3. **La SNP folA e' intergenica.** In `snps.tab` comparira' come `intergenic_region`, non
   come `missense`. Non filtrare per effetto sul CDS o la perdi.
4. **Verifica che 4,489,835 (coord. 0099) non cada in una regione mascherata** da PhiSpy o
   Gubbins: se il masking la rimuove, snippy-core non la riporta pur essendo reale.
   Con due soli genomi Gubbins e' comunque poco informativo (serve piu' di una coppia per
   stimare la ricombinazione) - per la coppia userei Snippy senza Gubbins e riserverei
   PhiSpy+Gubbins al contesto multi-isolato.
5. **Attesa sul numero.** Le 59 posizioni nell'array tandem quasi certamente scompaiono
   (coverage anomalo / multi-mapping). Delle 54 restanti, aspettati un sottoinsieme:
   il precedente della coppia 137 dice che le differenze fra consensus non regolate su
   read si riducono molto. Il numero da riportare nella lettera va aggiornato a quello
   Snippy, dichiarando il tool e le soglie.
6. **Indel.** La nostra conta esclude i 29 gap. Snippy li riporta (`snps.tab` tipo `del`/`ins`):
   vanno tenuti separati dal conteggio di sostituzioni.
