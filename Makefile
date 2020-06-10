new_species:
	@echo "Cloning species directory structure..."
	@mkdir new_species_directory
	@touch new_species_directory/Makefile
	@mkdir new_species_directory/zip
	@mkdir new_species_directory/bin
	@mkdir new_species_directory/sequencing
	@mkdir new_species_directory/sequencing/fprimers
	@mkdir new_species_directory/sequencing/rprimers
	@mkdir new_species_directory/sequencing/ref
	@mkdir new_species_directory/sequencing/fasta
	@mkdir new_species_directory/sequencing/fastq
	@mkdir new_species_directory/sequencing/genbank
	@mkdir new_species_directory/sequencing/protein
	@mkdir new_species_directory/sequencing/bed
	@mkdir new_species_directory/sequencing/sam
	@mkdir new_species_directory/sequencing/bam
	@mkdir new_species_directory/sequencing/vcf
	@mkdir new_species_directory/raw
	@mkdir new_species_directory/go
	@mkdir new_species_directory/src
	@mkdir new_species_directory/images
	@mkdir new_species_directory/misc
	@echo "Done."
