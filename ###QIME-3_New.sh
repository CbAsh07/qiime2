###QIME-3_New
#! /bin/bash

## Check inputs
if [[ $1 == "" ]]; then
	echo -e "Usage: eo_qiime2_UAB [manifest_file]"
	exit 1
else
	manifest_file=$1
fi

## Load Conda virtual env
echo -e "$(date -R): Activating conda env \"qiime2-amplicon-2024.2\""
eval "$(conda shell.bash hook)"
source activate /home/ethan/miniconda3/envs/qiime2-amplicon-2024.2

## Set global params
# cpus=$(( "$(nproc --all)" / 2))
# mem=$(awk '/MemFree/ { printf "%.0f \n", $2/1024/1024 }' /proc/meminfo)
# dir=$(pwd -P)

## Set script specific params

## Function definitions
# Steps here are based on "Moving Pictures" tutorial with some adjustment to 
# account for data specificity

qiime_import () {
# Import sequences to qiime. Manifest file should have: sample-id, 
# forward-absolute-filepath, reverse-absolute-filepath all tab separated

	echo -e "$(date -R): Running qiime import"

	qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path "${manifest_file}" \
	--output-path paired-end-demux.qza \
	--input-format PairedEndFastqManifestPhred33V2

	echo -e "$(date -R): Running qiime import visualization"

	qiime demux summarize \
	--i-data paired-end-demux.qza \
	--o-visualization demux.qzv

	echo -e "$(date -R): Done"

}

qiime_import_visualize () {
# Visualize quality scores, needed for dada2 quality control

	echo -e "$(date -R): Running qiime import visualization"

	qiime demux summarize \
	--i-data paired-end-demux.qza \
	--o-visualization demux.qzv

	echo -e "$(date -R): Done"

}

qiime_dada2 () {
# Sequence quality control, set --p-trunc-len-f/r (truncate position for forward
# and reverse reads) based on qiime2view of demux.qzv above

	echo -e "$(date -R): Review \"demux.qzv\" and select forward and reverse \
	truncation positons..."
	read -r -p "$(date -R): Enter forward read truncation position: " trunc_f
	read -r -p "$(date -R): Enter reverse read truncation position: " trunc_r
	
	# Default truncation positions based on manual analysis
	trunc_f=${trunc_f:-250}
	trunc_r=${trunc_r:-161}

	echo -e "$(date -R): Running qiime dada2 with forward truncation position: \
	${trunc_f} and reverse truncation postion ${trunc_r}"

	qiime dada2 denoise-paired \
	--i-demultiplexed-seqs paired-end-demux.qza \
	--p-trunc-len-f "${trunc_f}" \
	--p-trunc-len-r "${trunc_r}" \
	--o-table table-dada2.qza \
	--o-representative-sequences rep-seqs-dada2.qza \
	--o-denoising-stats stats-dada2.qza

	echo -e "$(date -R): Running qiime dada2 visualization"
	
	qiime metadata tabulate \
	--m-input-file stats-dada2.qza \
	--o-visualization stats-dada2.qzv

	echo -e "$(date -R): Done"
	
	echo -e "$(date -R): Constucting feature table with dada2-processed data"

	cp table-dada2.qza table.qza
	cp rep-seqs-dada2.qza rep-seqs.qza

	qiime feature-table summarize \
	--i-table table.qza \
	--o-visualization table.qzv \
	--m-sample-metadata-file sample-metadata.tsv
	qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv

	echo -e "$(date -R): Done"

}

qiime_dada2_visualize () {
# Visualized dada2 output

	echo -e "$(date -R): Running qiime dada2 visualization"
	
	qiime metadata tabulate \
	--m-input-file stats-dada2.qza \
	--o-visualization stats-dada2.qzv

	echo -e "$(date -R): Done"

}

qiime_dada2_feature_table () {
# Construct a feature table from denoised data to explore data

	cp table-dada2.qza table.qza
	cp rep-seqs-dada2.qza rep-seqs.qza

	qiime feature-table summarize \
	--i-table table.qza \
	--o-visualization table.qzv \
	--m-sample-metadata-file sample-metadata.tsv
	qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv

	echo -e "$(date -R): Done"

}

#####################################################################################

qiime_phylo_tree () {
# Align sequences to calculate a phylogenetic tree used for alpha/beta diversity

	echo -e "$(date -R): Running qiime phylo tree"

	qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences rep-seqs.qza \
	--o-alignment aligned-rep-seqs.qza \
	--o-masked-alignment masked-aligned-rep-seqs.qza \
	--o-tree unrooted-tree.qza \
	--o-rooted-tree rooted-tree.qza

	echo -e "$(date -R): Done"
}

qiime_table_filter () {
# Filter feature table to remove Odoribacter samples

	echo -e "$(date -R): Filtering qiime feature table"

	qiime feature-table filter-samples \
	--i-table table.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-where 'NOT [group]="Odoribacter"' \
	--o-filtered-table filtered-table.qza

	echo -e "$(date -R): Done"
}

qiime_core_metrics () {
# Calculate alpha and beta diversity. Sampling depth should be chosen based on 
# table.qzv above: want as high as possible while excluding as few samples as possible

	echo -e "$(date -R): Running qiime core metrics"

	qiime diversity core-metrics-phylogenetic \
	--i-phylogeny rooted-tree.qza \
	--i-table filtered-table.qza \
	--p-sampling-depth 50000 \
	--m-metadata-file sample-metadata.tsv \
	--output-dir core-metrics-results

	echo -e "$(date -R): Done"
}

qiime_alpha_diversity () {
# Calculate alpha diversity metrics for sample group comparisons

	echo -e "$(date -R): Running qiime alpha diversity"
	
	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/shannon_vector.qza \
	--m-metadata-file sample-metadata.tsv \
	--o-visualization core-metrics-results/shannon-group-significance.qzv

	# qiime diversity alpha-group-significance \
	# --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	# --m-metadata-file sample-metadata.tsv \
	# --o-visualization core-metrics-results/faith-pd-group-significance.qzv

	# qiime diversity alpha-group-significance \
	# --i-alpha-diversity core-metrics-results/evenness_vector.qza \
	# --m-metadata-file sample-metadata.tsv \
	# --o-visualization core-metrics-results/evenness-group-significance.qzv

	echo -e "$(date -R): Done"
}

qiime_beta_diversity () {	
# Calculate beta diversity metrics for sample group comparisons

	echo -e "$(date -R): Running qiime beta diversity on jaccard distance matrix"
	
	qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
	--m-metadata-file sample-metadata.tsv \
	--m-metadata-column sample-group \
	--o-visualization core-metrics-results/jaccard-sample-group-significance.qzv \
	--p-pairwise
	# qiime diversity beta-group-significance \
	# --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
	# --m-metadata-file sample-metadata.tsv \
	# --m-metadata-column cohort \
	# --o-visualization core-metrics-results/jaccard-cohort-significance.qzv \
	# --p-pairwise
	# qiime diversity beta-group-significance \
	# --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
	# --m-metadata-file sample-metadata.tsv \
	# --m-metadata-column alt-sample-group \
	# --o-visualization core-metrics-results/jaccard-alt-sample-group-significance.qzv \
	# --p-pairwise

	echo -e "$(date -R): Done"
}

qiime_alpha_rarefaction () {
# Check sampling depth for alpha diversity calculations

	echo -e "$(date -R): Running qiime alpha rarefaction"

	qiime diversity alpha-rarefaction \
	--i-table table.qza \
	--i-phylogeny rooted-tree.qza \
	--p-max-depth 10000 \
	--m-metadata-file sample-metadata.tsv \
	--o-visualization alpha-rarefaction.qzv

	echo -e "$(date -R): Done"

}

qiime_classifier () {
# Classify sequences using Greengenes 2022 10 database

	classifier="gg_2022_10_backbone_full_length.nb.qza"

	echo -e "$(date -R): Running qiime classifier with ${classifier}"

	qiime feature-classifier classify-sklearn \
	--i-classifier ${classifier} \
	--i-reads rep-seqs.qza \
	--o-classification taxonomy.qza

	qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv

	qiime taxa barplot \
	--i-table table.qza \
	--i-taxonomy taxonomy.qza \
	--m-metadata-file sample-metadata.tsv \
	--o-visualization taxa-bar-plots.qzv

	echo -e "$(date -R): Done"

}

qiime_ANCOM_BC_genus () {
# Detects differential abundance using ANCOM-BC method

	echo -e "$(date -R): Running qiime differential abundance with ANCOM-BC at Genus level"
	echo -e "$(date -R): Comparing sample-group"

	qiime taxa collapse \
	--i-table table.qza \
	--i-taxonomy taxonomy.qza \
	--p-level 6 \
	--o-collapsed-table table-genus.qza

	qiime composition ancombc \
	--i-table table-genus.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula 'group' \
	--o-differentials ancombc-sample-group.qza \
	qiime composition da-barplot \
	 --i-data ancombc-sample-group.qza \
	 --p-significance-threshold 0.05 \
	 --p-level-delimiter ';' \
	 --o-visualization da-barplot-sample-group-genus.qzv

	echo -e "$(date -R): Done"
	echo -e "$(date -R): Comparing alt-sample-group"

	qiime composition ancombc \
	--i-table table-genus.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula 'group2' \
	--o-differentials ancombc-alt-sample-group.qza 
	qiime composition da-barplot \
	 --i-data ancombc-alt-sample-group.qza \
	 --p-significance-threshold 0.05 \
	 --p-level-delimiter ';' \
	 --o-visualization da-barplot-alt-sample-group-genus.qzv

	echo -e "$(date -R): Done"
	echo -e "$(date -R): Comparing cohort"

	qiime composition ancombc \
	--i-table table-genus.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula 'cohort' \
	--o-differentials ancombc-cohort.qza 
	qiime composition da-barplot \
	 --i-data ancombc-cohort.qza \
	 --p-significance-threshold 0.05 \
	 --p-level-delimiter ';' \
	 --o-visualization da-barplot-cohort-genus.qzv

	 echo -e "$(date -R): Done"
}

qiime_ANCOM_BC_species () {
# Detects differential abundance using ANCOM-BC method

	echo -e "$(date -R): Running qiime differential abundance with ANCOM-BC at Species level"
	echo -e "$(date -R): Comparing sample-group"

	qiime taxa collapse \
	--i-table table.qza \
	--i-taxonomy taxonomy.qza \
	--p-level 7 \
	--o-collapsed-table table-species.qza

	qiime composition ancombc \
	--i-table table-species.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula 'group' \
	--o-differentials ancombc-sample-group.qza
	qiime composition da-barplot \
	 --i-data ancombc-sample-group.qza \
	 --p-significance-threshold 0.001 \
	 --p-level-delimiter ';' \
	 --o-visualization da-barplot-sample-group-species.qzv

	echo -e "$(date -R): Done"
	echo -e "$(date -R): Comparing alt-sample-group"

	qiime composition ancombc \
	--i-table table-species.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula 'group2' \
	--p-reference-levels group2::Alistipes-2 \
	--o-differentials ancombc-alt-sample-group.qza 
	qiime composition da-barplot \
	 --i-data ancombc-alt-sample-group.qza \
	 --p-significance-threshold 0.001 \
	 --p-level-delimiter ';' \
	 --o-visualization da-barplot-alt-sample-group-species.qzv

	echo -e "$(date -R): Done"
	# echo -e "$(date -R): Comparing cohort"

	# qiime composition ancombc \
	# --i-table table-species.qza \
	# --m-metadata-file sample-metadata.tsv \
	# --p-formula 'cohort' \
	# --o-differentials ancombc-cohort.qza 
	# qiime composition da-barplot \
	#  --i-data ancombc-cohort.qza \
	#  --p-significance-threshold 0.001 \
	#  --p-level-delimiter ';' \
	#  --o-visualization da-barplot-cohort-species.qzv

	#  echo -e "$(date -R): Done"
}

qiime_diversity_significance () {

	qiime diversity adonis \
	--i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
	--m-metadata-file sample-metadata.tsv \
	--p-formula "group2" \
	--o-visualization core-metrics-results/jaccard_alt-sample-group-by-cohort_adonis.qzv

}

qiime_phyloseq_convert () {

	# Export OTU table
	mkdir -p phyloseq
	qiime tools export \
	--input-path table-deblur.qza \
	--output-path phyloseq

	# Convert biom format to tsv format
	biom convert \
	-i phyloseq/feature-table.biom \
	-o phyloseq/otu_table.tsv \
	--to-tsv
	cd phyloseq || exit
	sed -i '1d' otu_table.tsv
	sed -i 's/#OTU ID//' otu_table.tsv
	cd ../

	# Export representative sequences
	qiime tools export \
	--input-path rep-seqs-deblur.qza \
	--output-path phyloseq

	# Export tree files
	qiime tools export \
	--input-path unrooted-tree.qza \
	--output-path phyloseq
	cd phyloseq || exit
	mv tree.nwk unrooted_tree.nwk
	cd ../

	qiime tools export \
	--input-path rooted-tree.qza \
	--output-path phyloseq
	cd phyloseq || exit
	mv tree.nwk rooted_tree.nwk

	echo -e "$(date -R): Done"

}

rdp_classify () {

	classifier_path="/home/ethan/work/rdp_classifier/rdp_classifier_2.13/dist/classifier.jar"
	sequences_path=$(find . -name "dna-sequences.fasta" -exec realpath {} \;)

	java -Xmx4g -jar "${classifier_path}" classify \
	-c 0.5 \
	-f db \
	-g 16srrna \
	-h heirarchy_table.tsv \
	-o db_classification_table.tsv \
	"${sequences_path}"

	echo -e "$(date -R): Done"

}

## Main

## Data import
# qiime_import
## qiime_import_visualize # Incorporated into qiime_import

## Data processing
# qiime_dada2
## qiime_dada2_visualize # Incorporated into qiime_dada2
## qiime_dada2_feature_table # Incorporated into qiime_dada2

## Data analyis
# qiime_phylo_tree
# qiime_table_filter
# qiime_core_metrics
# qiime_alpha_diversity
# qiime_alpha_rarefaction
# qiime_beta_diversity
# qiime_classifier
# qiime_ANCOM_BC_genus


qiime_ANCOM_BC_species


# qiime_diversity_significance


# rdp_classify

## Data export
# qiime_phyloseq_convert

######################################################################
##########                      Notes                       ##########
######################################################################

# java -Xmx4g -jar "${classifier_path}" classify
# usage:  [options] <samplefile>[,idmappingfile] ...
#  -b,--bootstrap_outfile <arg>   the output file containing the number of
#                                 matching assignments out of 100 bootstraps for
#                                 major ranks. Default is null
#  -c,--conf <arg>                assignment confidence cutoff used to determine
#                                 the assignment count for each taxon. Range
#                                 [0-1], Default is 0.8.
#  -d,--metadata <arg>            the tab delimited metadata file for the samples,
#                                 with first row containing attribute name and
#                                 first column containing the sample name
#  -f,--format <arg>              tab-delimited output format:
#                                 [allrank|fixrank|biom|filterbyconf|db]. Default
#                                 is allRank.
#                                 allrank: outputs the results for all ranks
#                                 applied for each sequence: seqname, orientation,
#                                 taxon name, rank, conf, ...
#                                 fixrank: only outputs the results for fixed
#                                 ranks in order: domain, phylum, class, order,
#                                 family, genus
#                                 biom: outputs rich dense biom format if OTU or
#                                 metadata provided
#                                 filterbyconf: only outputs the results for major
#                                 ranks as in fixrank, results below the
#                                 confidence cutoff were bin to a higher rank
#                                 unclassified_node
#                                 db: outputs the seqname, trainset_no, tax_id,
#                                 conf.
#  -g,--gene <arg>                16srrna, fungallsu, fungalits_warcup,
#                                 fungalits_unite. Default is 16srrna. This option
#                                 can be overwritten by -t option
#  -h,--hier_outfile <arg>        tab-delimited output file containing the
#                                 assignment count for each taxon in the
#                                 hierarchical format. Default is null.
#  -m,--biomFile <arg>            the input clluster biom file. The classification
#                                 result will replace the taxonomy of the
#                                 corresponding cluster id.
#  -o,--outputFile <arg>          tab-delimited text output file for
#                                 classification assignment.
#  -q,--queryFile                 legacy option, no longer needed
#  -s,--shortseq_outfile <arg>    the output file containing the sequence names
#                                 that are too short to be classified
#  -t,--train_propfile <arg>      property file containing the mapping of the
#                                 training files if not using the default. Note:
#                                 the training files and the property file should
#                                 be in the same directory.
#  -w,--minWords <arg>            minimum number of words for each bootstrap
#                                 trial. Default(maximum) is 1/8 of the words of
#                                 each sequence. Minimum is 5

