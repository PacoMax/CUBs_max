#!/usr/bin/env bash
if [ "$1" == "-h" ]; then
echo "Welcome to CUBs_max.sh part of the CUBs_go_max pipeline"
echo "To run this program it's necesary to install these programs:"
echo ""
echo "	ENCprime-master	link"
echo "	codonW	link"
echo "	infernal-1.1.2-linux-intel-gcc	link"
echo "	hmmer-3.2.1	link"
echo "	tRNAscan-SE-2.0	link"
echo ""
echo "Usage: `basename $0` [king] [code] [list1] [list2] [sizes]"
echo ""
	echo "	king:"
		echo  "		E	: search for eukaryotic tRNAs"
		echo  "		B	: search for bacterial tRNAs"
		echo  "		A	: search for archaeal tRNAs"
	echo "	code:"
		echo  "		Number between 1 and 25 indicates a Genbank code"
	echo "	list1:"
		echo  "		Folder that contain the genomes"
	echo "	list2:"
		echo  "		Folder that contain the cds"
        echo "  sizes:"
                echo  "         The size in CDS sequences to cut the files contaning the CDS sequences"
		echo  "		This option depends on the memory ram of your system"
	echo ""
	echo " If it was not successfully running, it's necessary to delete all created directories"
echo ""
exit 0
fi

king=${1?Error: No superkingdom specified. Please, ask for help (CUBs_max.sh -h)}
code=${2?Error: No genetic code specified. Please, ask for help (CUBs_max.sh -h)}
list1=${3?Error: No genome folder input. Please, ask for help (CUBs_max.sh -h)}
list2=${4?Error: No CDS folder input. Please, ask for help (CUBs_max.sh -h)}
sizes=${5?Error: No size for CDS input selected. Please, ask for help (CUBs_max.sh -h)}

mkdir trnas_tables
mkdir ENC_table
mkdir triplets
mkdir CU_cvs

if ! { cd $list1; }; then
    echo "$0: Cannot find the directory $list1"
    exit 1
fi
for genome in $(ls)
	do
	tRNAscan-SE -$king --codons -Q -o# $genome
	get_trnas_from_trnascan.pl -i $genome.out > ../trnas_tables/$genome.out.trna
	done

if ! { cd ../$list2; }; then
    echo "$0: Cannot find the directory $list2"
    exit 1
fi
for cds in $(ls)
	do
	perl -pe '/^>/ ? print "\n" : chomp' $cds | tail -n +2 > $cds.fasta
        get_table_protein.pl -f ${cds}.fasta > ${cds}.table.txt
        cut -d " " -f1 ${cds}.fasta | sed 's/>lcl|/>/g' > ${cds}.fasta_cut
        codonw ${cds}.fasta_cut -all_indices -nomenu -silent
        cut -f2-15 ${cds}.out | paste ${cds}.table.txt - | tr ' ' '_'| sed 's/\t\t/\t/g' | sed 's/\*\*\*\*\*/NA/g' > ../CU_cvs/${cds}.table_codon_usage.cvs
	awk -v size=$sizes -v pre="$cds.r" -v pad=7 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' ${cds}.fasta_cut
	touch $cds._results_r
	for cds2 in $(ls | grep $cds.r)
		do
		l=$(grep -c ">" $cds2)
		SeqCount -c $cds2 $l
		SeqCount -n $cds2 $l
		ENCprime $cds2.codcnt $cds2.acgtfreq $code $cds2.results 1 -q
		cat $cds2.results >> $cds._results_r
		rm $cds2.*codcnt
		rm $cds2.*acgtfreq
		rm $cds2.results
		done
	grep -v "Totals:" $cds._results_r | grep -v "Name Nc Ncp" | sed 's/://g' > ../ENC_table/ENC_${cds%_cds_from_genomic.fna.fasta_cut.results}_table
	codonM.pl $cds.fasta_cut ../triplets/${cds}.table_triplets.txt
	done
