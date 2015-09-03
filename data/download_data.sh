# download GreenGenes reference sequences
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz

# download initial primer pairs from Klindworth et al., NAR 2013
wget http://www.arb-silva.de/fileadmin/silva_databases/primer_evaluation/supplementary_tables/supplementary_table_01.xlsx
wget http://www.arb-silva.de/fileadmin/silva_databases/primer_evaluation/supplementary_tables/supplementary_table_08.xlsx

# choose cluster similarity level(s) of the GreenGenes OTUs that will be used as representative set(s)
# possible similarity levels: 61 64 67 70 73 76 79 82 85 88 91 94 97 99
SIM="85"

gzip -d gg_13_5_updated_tax.txt.gz
for lev in ${SIM}; do
	# extract the representative set and taxonomy files from the archive
	echo "Extracting files for level "${lev}
	tar -xf gg_13_5_otus.tar.gz gg_13_5_otus/rep_set/${lev}_otus.fasta
	tar -xf gg_13_5_otus.tar.gz gg_13_5_otus/taxonomy/${lev}_otu_taxonomy.txt
	
	# correct taxonomy file based on our updated taxonomy
	echo "Updating taxonomy for level "${lev}
	sed 's/; /	/g' gg_13_5_otus/taxonomy/${lev}_otu_taxonomy.txt > gg_13_5_otus/taxonomy/${lev}_tmp_tax.txt
	sort -n gg_13_5_otus/taxonomy/${lev}_tmp_tax.txt > gg_13_5_otus/taxonomy/${lev}_sort_tax.txt
	awk 'BEGIN{FS="\t";i=0;j=0};NR==FNR{a[i++]=$1};(NR!=FNR && $1==a[j]){j++;print $0}' gg_13_5_otus/taxonomy/${lev}_sort_tax.txt gg_13_5_updated_tax.txt > gg_13_5_otus/taxonomy/${lev}_upd_tax.txt
	rm gg_13_5_otus/taxonomy/${lev}_sort_tax.txt gg_13_5_otus/taxonomy/${lev}_tmp_tax.txt
done
gzip gg_13_5_updated_tax.txt
