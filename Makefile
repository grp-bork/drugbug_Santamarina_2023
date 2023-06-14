

all: data/table_sensitization_protection.tsv figures/fig1.png figures/fig2.png figures/fig3.png figures/fig4.png figures/fig5.png \
	suppl_figures/QC_param_evaluation_pvalue_concentration_dependency.png suppl_figures/suppl_figure_conditions_colon_concentrations.pdf suppl_figures/QC_corrs.pdf \
	suppl_figures/nifurtimox_classification.png suppl_figures/nifurtimox1.pdf

suppl_figures/nifurtimox1.pdf: data/Nifur_GM.tsv data/AUCs_Supernatants_Nifur_GM.tsv
	./suppl_fig_nifurtimox.R

suppl_figures/nifurtimox_classification.png: data/mean_effects.tsv
	./nifurtimox_classification.R

suppl_figures/sfig_rare_species.pdf: data/16S_counts.tsv data/combined_monoculture_aucs.tsv
	./rare_species.R

suppl_figures/QC_corrs.pdf:
	./QC_16S.R

suppl_figures/QC_param_evaluation_pvalue_concentration_dependency.png: data/effects.tsv data/combined_monoculture_aucs.tsv
	./QC_thresholds.R

suppl_figures/suppl_figure_conditions_colon_concentrations.pdf: data/drug_suppl_data.tsv data/metabolomics_conditions.tsv data/16S_counts.tsv data/Plates_Layout.tsv data/combined_monoculture_aucs.tsv data/Maier_et_al_hits.tsv data/closest_to_gut_concentration.tsv data/table_sensitization_protection.tsv
	./suppl_table_drugs.R


figures/fig1.png: illustrations/fig1_drugbug_community_schema.png panels/fig1_rel_abundance_barchart.png panels/fig1_drug_vs_control_example.png panels/fig1_single_species_growth_curves.png illustrations/fig1_monoculture_schema.png
	./fig1_combine_panels.R

panels/fig1_rel_abundance_barchart.png: data/counts.tsv
	./fig1_example_abundance.R

panels/fig2_classification_explanation_color_bg.png: panels/fig1_drug_vs_control_example.png

panels/fig1_drug_vs_control_example.png: data/mean_effects.tsv
	./fig2_classification_explanation.R

panels/fig1_single_species_growth_curves.png: data/single_species_growth_curves.tsv 
	./fig1_growth_curves.R


figures/fig2.png: panels/fig2_classification_explanation_text.png panels/fig2_classification_explanation_color_bg.png panels/fig2_classification_closest_conc_barchart.png panels/fig2_concentration_boxplot_Protection.png
	./fig2_combine_panels.R

panels/fig2_classification_explanation_text.png:
	./fig2_classification_textual_explanation.R

data/table_sensitization_protection.tsv: panels/fig2_classification_closest_conc_barchart.png

panels/fig2_classification_closest_conc_barchart.png: data/effects.tsv data/effect_counts_hit.tsv data/species_effect_counts_hit.tsv data/closest_to_gut_concentration.tsv data/metabolomics_conditions.tsv
	./fig2_classification.R

panels/fig2_concentration_boxplot_Protection.png: data/effect_counts_hit.tsv
	./fig2_plot_protection_sensitisation_vs_concentration.R


figures/fig3.png: panels/fig3_metab_examples.png panels/fig3_metab_explanation.png panels/fig3_outcome_protected_2.5.png
	./fig3_combine_panels.R

panels/fig3_metab_examples.png: panels/fig3_outcome_protected_2.5.png

panels/fig3_metab_explanation.png: panels/fig3_outcome_protected_2.5.png

panels/fig3_outcome_protected_2.5.png: data/metabolomics_treatment_normalised.tsv data/metabolomics_MGAM_normalised.tsv data/metab_effect_correlation.tsv data/metab_effects.tsv
	./fig3_metab.R


figures/fig4.png: panels/fig4_c_comes_community.png panels/fig4_metabolomics_growth1.png panels/fig4_metabolomics_growth2.png panels/fig4_niclosamide_classification.png illustrations/fig4_schema.png illustrations/fig4_schema_ccomes.png panels/fig4_ecoli_metabolomics.png
	./fig4_combine_panels.R

panels/fig4_ecoli_metabolomics.png: data/ecoli_metabolomics.tsv
	./fig4_ecoli_metabolomics.R

panels/fig4_c_comes_community.png: data/c_comes_community.tsv
	./fig4_c_comes_community.R

panels/fig4_metabolomics_growth2.png: panels/fig4_metabolomics_growth1.png

panels/fig4_metabolomics_growth1.png: data/niclosamide_metabolomics/Metabolomics.tsv data/niclosamide_metabolomics/AUCs_Supernatants_Matching_Metabolomics.tsv
	./fig4_metabolomics_growth.R

panels/fig4_niclosamide_classification.png: data/mean_effects.tsv
	./fig4_classification.R


figures/fig5.png: panels/fig5_nitroreductase_overexpression.png panels/fig5_proteomics.png illustrations/fig5_schema.png panels/fig5_growth.png
	./fig5_combine_panels.R

panels/fig5_nitroreductase_overexpression.png: data/nitroreductase_overexpression.tsv
	./fig5_overexpression_mics.R 

panels/fig5_proteomics.png: data/raw\ data_Bvulgatus_R.xlsx data/Roseburia_raw_R_sorted\ by\ log\ 10\ untreated.xlsx 
	./fig5_proteomics.R

panels/fig5_growth.png: data/Bv_OE_supernatants.tsv
	./fig5_growth.R


.DELETE_ON_ERROR:

.SECONDARY:

