# Annotation of Eukaryotic Genomes: Arabidposis thaliana accession
This course directly follows the course on Genome and Transcriptome assembly. Here, we are going to annotate the HiFiasm assembly produced during the previous class. We are also going to use the Trinity transcriptome assembly. 

## Transposable Elements (TE) Annotation 
TEs are moving elements of great importance in genome evolution, and overly represented in plant genomes. Inside one genome, one can derive orders, superfamilies, families and clades of TEs forming their own phylogeny through structure similarity. 

They are broadly categorized into Class I DNA retrotransposons and Class II DNA transposons due to their difference in transposition mechanisms:
- Class I TEs have a copy-paste mecahnism, using an RNA intermediate to increase their copy number in the genome. 
- Class II TEs have a cut-adn-past mechanism, transposing directly via DNA. 

In this part, we are going to use Extensive de novc TE Annotator EDTA to annotate TEs in our Arabidopsis accession and classify them using TEsorter. Since eukaryotic genomes are repeat-rich, it is necessary to mask them up front such that downstream tools are not misinterpreting homologies and falsy alignemnts that can arise from them. This allows for better homology and orthology annotation later, but the biology of TEs themselves can be very informative, notably through their effect on genome size and the recentness of their activity in genome evolution. 

### Annotating TEs with EDTA
To use EDTA for TE annotation:

Run script 01a_run_EDTA.sh 

Here, we are using the options: 
--genome: the HiFiasm genome that was previously assembled
--species others: a generic model 
--step all: to run the full EDTA pipeline
--sensitive 1: to use RepeatModeler which allows for higher sensitivity in TE discovery
--cds: the coding sequences of Arabidopsis thaliana reference to avoid labeling true genes as TEs.

Running the EDTA pipeline produces several files, of which 5 are going to be used in downstream analysis: 
- A .mod.EDTA.TElib.fa file: FASTA file with the TE classified up to the superfamily, where each sequence represents a TE family.
- A .mod.EDTA.TEanno.gff3 file: GFF file with both intact and fragmented TE annotations across the genome. 
- A .mod.EDTA.intact.gff3 file: Similarly, a GFF3 only containing intact TE annotations, which is helpful to identify recently active or well-conserved TEs. 
- A .mod.EDTA.TEanno.sum file: Summary file with the annotated TEs, their counts and number of bp in the genome per superfamily and family
- A .genome.mod.out inside the .mod.EDTA.anno subfolder: Output from RepeatMasker with the details information on TE copies, and namely their percentage of diversity compared to the reference sequences which is useful for dating TE insertions. 

#### Annotation of LTR-RT candidates
Once EDTA has finised the identification of raw LTR candidates, it is possible to directly classify these into known clades and superfamilies using TEsorter/REXdb. 
An LTR-RT inserts with two identical LTRs on 5' and 3'. Over time, these can accumulate substitutions, such that two LTRs diverge. 
The EDTA pipeline contains discovery tools, such as LTRharvest and LTR_FINDER to find candidate LTR-RTs based on the similar direct repeats flanking an internal region. The, it goes through a curation step with LTR_retriever to validate candidates and provides an identity measure by aligning the two LTR sequences for each element, computing the identity as a measure of the number of matchin aligned bases divided by the alignment length. Then, EDTA writes it as ltr_identity in the GFF, scoring from 0 to 1, where 1 is a very recent insertion, and below 0.95 is considered an older insertion as more substitutions have accumulated. 
To study LTR-RTs more in detail, we can also assign superfamily and clades using TEsorter: 

Run 01b_run_TEsorter.sh

TEsorter takes as input the internal region of the LTR, and searched conserved protein domains such as GAG, proteases, RT, RH and such usinf HMMs profiles from REXdb. The best-scoring domain set and domain order give the superfamily assigment (Copia vs Gypsy), while more finer HMMs and motifs allow for clade classification. If the domains are too partial or degraded, they are assigned as unknown or NA. 

The output files are canonical classification table .cls.tsv, and per-superfamily, per-clade counts that can be used for plotting bar charts by running 01c_full_length_LTRs_identity.R script provided with the course material. The script was slightly modified to match our working environment. To launch it, use the helper script:

Run 01c_run_clade_plots.sh

It produces a series of histograms of LTR pair identity on the x-axis, for both Copia and Gypsy clades. y-axis shows counts per bin, capped at 8 (there might by more copies in total). An identity of 1.00 signifies very recent insertions, whereas LTR-RTs with identities below 0.95 signify older insertions. It is to be noted that quick cout reports showed 60 Copia and 85 Gypsy LTRs, but also 436 NA, which explains the empty space as many LTRs lacked a clade label. 
We can appreciate that, for our TE set, there are clear differences with Gypsy clades more present than Copia. Nonetheless, both show many clades (Ale, Ivana, SIRE and Tork for Copia, Athila, Reina, Retand and Tekay for Gypsy) with very recent insertions (identity at max 1.00) and less clades contained older insertions (notably Bianca for Copia and Retand, Reina and Athila for Gypsy, with identities between 0.85-0.90). This suggests that the accession's genome is not fixed, but rather ungoing TE-host dynamics, with lineage heterogeneity as shown by the difference in age and abundance profiles between Copia and Gypsy adn well as clades within them. 

#### Visualizing TE superfamilies across contigs
Once EDTA has ran completely, the .mod.EDTA.TEanno.sum is outputed and gives indications as per the total TE content of our genome's accession. we can see that there is a total of 22528520 bp that are interspersed TEs which represents 15.14% of the genome (148809274 bp). The number of hits was 32572. These are all the TEs EDTA has masked for downstream analysis. 
The most abundant superfamily is LTR unknown category, which represent 6.29% of the total genome. Helitrons (3.21%) and Gypsy (2.05%) are next, followed by Copia (0.64%). This is consistent with previous LTR analysis showing a higher abundance of Gypsy clades compared to Copia. 

To visualize the distirbution of TEs across our contigs, we use the circlize R package. The following wrapper scripts creates the FASTA indexation of our assembly such that it can be plotted through the R script 02b_circlize_te_density.R

Run 02a_run_circlize.sh

The output plots show slices representing the 10 longest scaffolds clockwise, where the ticks mark the Mb. Stacked tracks or coloured rings represent a different TE superfamily. The height of a ring at a given position is its density in fraction of basis in that 100 kbo window annotated a superfamily, the taller the denser. Density is relative per window and different rings are not normalized to the same maximum, so only shapes and peaks should be compared at given positions rather than absolute heights across tracks. 
We can clearly appreciate scaffold sectors where serveral rings flare together, some overlaps and differences between the different TE families distributions across the contigs. 

#### Clade-level classifiction based on protein domains 
Using TEsorter, we can refine our analysis on LTR-RTs Copia and Gypsy using the .mod.EDTA.TElib.fa consensus TE library file outputed by EDTA, which contains all the family representative with mxied classes of TE. 

The first step is to subset Copy and Gypsy squences to speed the TEsorter as well as avoid domains hits from non-LTR elements, and make clade calls interpretable. 

Run 03a_extract_copia_gypsy.sh
This extract consensus sequences FASTA files for Copia and Gypsy superfamilies.

Run 03b_tesorter_classification_count.sh
Taking the Copia and Gypsy consensus FASTAS previously extracted from the EDTA's TE library, we run TEsorter Rexdb-plant to assign clade-level labels o those sequences. Then, we normalize heterogeneous TEsorter TSV into a TE/SUPERFAMILY/CLADE hierarchy table and map thes clade labels bock onto our genome-annotation GFF to count how many annotated elements fall into each clade.
The main run produces .rexdb-plant-dom-faa that contains annotated protein sequences that can be used for TE phylogenetic analysis, as well as rexdb-plant.cls.tsv that contains the actual classification of TEs into their classes, orders and families, among other files. Moreover, we parse these files to produce other useful TSVs like tesorter_all_class.tsv serving as master classification table, and count files: of conesnesus sequences per clade, per superfamily and clade pair and of annotated fenome elements per clade. 

This time, only 35 Gypsy clades and 26 Copia elements were found. Here, using the complete TEanno library gives a measure of abundance, whereas previous LTR library can tell about recency.  

#### Dating and estimatio of insertion age of TE
The age of a TE, which refers to the time of its insertion, can be estimated by measuring the divergence of its sequence from a consensus one, under the assumtion that copies of the same TE family originated from a single active element and that at the time of insertion, all copies are identical but over time, they accumulate mutations. The more divergent the copy becomes, the older it can be. 
For this, we are going to exploit the RepeatMasker's table of TE that were found in our genome by EDTA. 

Run script 04a_parseRM.sh

This script is a wrapper that allows to launch the 05-parseRM.pl file given in the course material. It allows binning 0-50% divergence in 1% steps. 
Perl reads each TE copy from the .out file and normalizes divergene from the consensus using CpG hypermutability model. Each copy is bucketed into divergence bins from 0 to 50%, and are aggregated by superfamily. This allows for the TE landscape table creation that can the be plotted. 

Run 04b_plot_TE_landscape.sh

This script is a wrapper that launches the R script that was provided through course material to plot the TE divergence binning produced by the previous step. It has been slighlty modified to accomodat our specific working environement as well as to include dating estimates by using the estimated substitution rate for Brassicaceae as r=8.22*10^-9 substitutions per synonymous site per year and the insertion time as T = K/2r, K=sequence divergence. 
Recent peaks can my appreciated, with Copia spiking at 0.9 Myr. Gypsy appears to dominate between 2-10 Myr, with broad peaking and lots of mass. Overall, Gypsy dynamics are more prominent than Copia, with Copia haveing a very recent burst and a substantial ancient load while Gypsy being more broad, mid-age. Recent Copia activity is thus detectable, but small, most Copia copies are old and ypsy is the main actor, while also accumulating older copies. 

## Annotation of genes with the MAKER pipeline
MAKER integrates multiple source of evidence for eukaryotic genome annotation. Namely, it contains ab initio prediction models using tools like Augustus and GeneMark to prdict gene structures based on sequence patterns, but also accepts RNA-seq data which provides evidence of gene expression as well as protein homology-based annotation which uses known protein sequences from related species to predict gene structures with, for example, BLAST to identify conserved genes through alignement between transcripts and known proteins to genome. 
To begin, we need to generate MAKER control files: 

Run 05a_maker_setup_ctl.sh

This scripts creates an annotation working directory and runs maker -CTL to generate the control files. Then, it edits maker_opts.ctl to point to our genome and Trinity transcripts, as well as protien evidence from TAIR and UniProt databases. It disables DFam masking and uses the EDTA TElib for RepeatMasker, and enraples RepearRunner TE proteins. It also uses the pre-trained arabidopsis model from the Augustus database for de novo protein prediction. Since we are going to run MAKER with MPI, the CPUs are kept to 1. 
Once the control files are generated, we can run MAKER:

Run 05b_run_MAKER.sh

The script runs MAKER under MPI, a way to run many separate processes by sending messages; here, we are launching 50 MAKER processes passing their results via files.
Then, MAKER splits the assembly into manageable chunks per contig and built a datastore. Those are hashed buckets, each containing a run directory for the contigs.
For each contig, it:
- Masked repeats using the EDTA TE library
- Searched RNA-seq transcripts and protein homology (BLAST/Exonerate) to anchor evidence
- Ran ab initio preictors using Augustus with arabidopsis species model 
- Combined evidence and predictions into gene models, even alternatice splicing when supported. 
- Wrote per-contig outputs (GFFs and per-contig FASTAS) inside datastore. 
However, MAKER does not concatenate per-contig results: we need to merge them at the end using a helper script as well as the master index that records, for every contig chunk, the path to its run directory as well as its status. 

Run 05c_maker_merge.sh

The merge script will read the master index as well as 
- gff3_merge: stitches all per-contig GFF2 into a single, consistent GFF3 with options: 
    - -s: include embedded FASTA at the end of the GFF which is handy but can produce big files
    - -n: noseq, GFF without that embedded FASTA which is sometime preferable. 
- fasta_merge: grabs predicted transcripts and proteins from all contigs into two FASTAs.

There are 31767 gene loci predicted by maker (features of type gene) and 37761 mRNA models (highly suppoorted: 36322 AED ≤ 0.5, 96.2%) which suggests that there are about 1.19 isoforms/gene. (Cross-checked with FASTA files: 37761 proteines and 37761 mRNAs)
TAIR10 protein coding gene are around 27k; we are about 18% above the reference. Inflations drivers can by allele duplication, fragmented model split by repeats, and cross-species evidence from UniProt Virididplantae.
Next, we are going to refine the gene models using InterProScan to annotate domains from Pfam, superfamilies and mapping these IDs to GO terms. 
With this, we can flag TE-like proteins that contain RT, integrases and transposes to purge them further, give support to weak AED models with conserved models and get function labels. 
First, we need to create an ID map from the NOSEQ GFF:

Run 05d_maker_rename_ids.sh

Running this sscript allows to create connsistent IDs looking like IST0000001 across GFF and FASTAs that allows InterProScan to annotate the protein FASTA and later inject those results back into the GFF, which requires matching IDs!

Run 05e_run_interproscan.sh 

After gene domain annotation,  29975 / 37761 (79.4%) proteins got at least one hit and are annotated. mRNAs with any InterPro/Pfam tag in the GFF summed to 29975, and 18869 had mRNAs with GO terms added. 

Now to assess the quality of our predictions and annotations, we can measure the Anootation Edit Distance AED, that ranges from 0 to 1 where 0 is perfect evidnece support and 1 no support. This measure is computed per mRNA by comparing the predicted features to aligned evidence from RNA-seq assemblies and proteins. It is an overlap of sensitivity (SN), the fraction of evidence covered by the model and precision (PPV), the fraction of the model covered by evidence such taht AED = 1- (SN + PPV)/2.

Run 05f_run_AED.sh

This script uses AED_cdf_generator.pl on our functional, renamed GFF to produce a CDF file. Then, it extract raw values from mRNA feature attributes and prints highly supported AEDs (≤ 0.5).96.2% (36322/37761) of proteins were well supported. To visualize the AED distirbution:

Run 05f_plot_aed.R

This produces an histogram showing the distribution from 0 to 1 of AEDs, which is pleasingly left-skewed. The ECDF shows a steep increase, supporting our models. 
Now that we have AED values, we can filter our GFF file for quality by leaving out poorly supported gene models: 

Run 05g_filter_gff.sh

The script retains transcripts with AED < 1 and/or Pfam domain through quality_filter.pl -s to get a clean annotation. 

Finally, we need to keep FASTA entries for mRNAs that passed the quality filter we applied above, such that we have clean transcript/protein sets for downstream analysis. 

Run 05h_filter_FASTA.sh

Now the FASTA file matches the filtered GFF content.

To further refine and validate MAKER genes, we can use AED to filter out genes that are not well-supported ≥ 0.5, remove TE-like models through InterPro/ Pfam annotation and keep models with functional domains and evidence. Further, cross-checking with BUSCO cna be informatice of genome completeness and duplication/ fragmentation.
AED allows to quantify the agreement between the predicted model and the evidence. The lower the AED score, the better, and generally ≤ 0.5 is accepted as good support.
InterPro/ GO give domain and GO hits that provide biological plausaibility of predicted proteins, helping to flag TE-like seuqences for removal and giving functional summaries such as pathways and families that can also be used for comparative analyses. 

## Quality Assessment of Gene Annotations

### BUSCO: Quality Asssessment of Gene Annotations
To assess the quality og our gene annotation so far, we are going to use BUSCO (Benchmarking Universal Single-Copy Orthologs) to infer the completness of our annotation theough the presence/ absence of highly conserved orthologs across a range of taxa. 

Since BUSCO expects one ortholog per gene, keeping multiple isoforms per gene would inflate the Duplicated category falsely. Thus, we are going to first extract the longest, assumed most complete coding isoform from our filtered FASTA files. 

Run 06a_extract_longest_isoforms.sh

Once the longest isoforms for proteins and transcripts have been extracted, 

Run 06b_run_BUSCO.sh

This allows to quantify completeness as BUSCO check our longest isoforms against a curated set of single-copy orthologs for Brassicales. We run both protein and transcriptome mode; protein mode amtches HMM profiles and domains and is insensitive to UTRs, giving an idea on the annotation quality of the coding sequences. The transcriptome mode works on nucleotide transcripts and tuhs toeralets missing ORFs, but is sensitive to UTRs. 
To get the associated bar plots, we adapted the R script from the BUSCO documentation:

Run 06b_busco_plot.R

Comparisons can be made before and after annotation of the assemblies (both hifiasm and trinity). 
BUSCO metrics are C (Complete), split into S (single-copy) and D (Duplicated), F (Fragmented), M (Missing) and n (total BUSCOs in lineage dataset). Generally, we want high C and low F/M, and high D shows redundancy due to polyploidy or mRNA expression in the Trinity assembly. 

After this, useful summary statistics regarging the quality of our annotation can be retrieved with AGAT, Another GFF/GTF Analysis Toolkit, which summarizes annotation features like genes, mRNA, exons, CDS, UTRs, introns lengths, counts and distirbution. It is useful for QC the structure of the MAKER gene set, quantify its composition and also compare versions pre/ post-filtering, for instance. 

Run 06c_run_AGAT.sh

Now, it is possible to add a track for the genes that we have annotated on top of our circos plots: 

Run 06d_circlize_QC_annotation.R

To grasp what we have been doing so far, the assembly and its annotations can be visualized through a genome browser like Geneious. Download the FASTA file of the HiFiasm assembly (istisu1.bp.p_ctg.fa) as well as the annotated GFF files (filtered.genes.renamed.gff3 and istisu1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3). Load the fasta file into Geneious, select all the contigs to group them as list and drop the GFF files on top. Additionally, the evidence used for gene prediciton can also be tracked (istisu1.bp.p_ctg.all.noseq.gff.renamed.functional_iprscan_quality_filtered.gff). It is also possible to visualize the transcript data: 

Run 06e_align_transcript.sh

This produces a BAM file that can be added to Geneious track.









## Functional Annotation
From the predicted proteins from MAKER and the filtered GFF describing gene models, which are structurally correct and well-supported provided their AEDs, the aim is to attach biological meaning by comparing predicted proteins to curated references. 
We are using two complementary references: 
- UniProt (Viridiplantae, reviewed): contains high-quality, manually curated proteins with names, functions and GO terms. 
- TAIR10 representative models: gold-standard gene set, great for closest Arabidopses otholog



First, we are using BLASTP to compare our protein sequences against an UniProt curated database. BLASTP will find sequence similarity: hits with strong scores (which have low E-value and high bitscore) indicate homology to known proteins. If they are similar, we can transfer putative function to our model. This only tells us about the likely function of an anonymous ORF. 
From the hits returned by BLAST, we are sorting by query ID and by bitscore in descending order. The first per query is kept, since it is the most confident homolog. We want to keep on name per protein for downstream labelling for the standard summary. 
To inject the function back into our files, we are using two MAKER helpers: maker_functional_fasta which takes UniProt FASTA, the besthits and our protein FASTA and outputs a new protein FASTA where headers include UniProt names and accession to see the putative functions in the header, and maker_functional_gff which takes UniProt FASTA, the full BLAST table and our GFF and outputs a GFF with extra attributes, reflecting UniProt info so genome browsers display functions alongside gene models. 

