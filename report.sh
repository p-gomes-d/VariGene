#!/bin/bash

awk 'NR>'$(grep -n "#CHROM" $(ls *_annotated.vcf) | cut -d ":" -f 1)' {print}' $(ls *_annotated.vcf) > no_header

#all variants
all_count=$(awk '{print}' no_header | wc -l)
#known variants
known_count=$(awk '{print}' no_header | grep "ALLELEID" | wc -l)
#pathogenic variants
patho_count=$(awk '{print}' no_header | grep -o "CLNSIG=[^;]*" | grep "Pathogenic" | wc -l)
#benign variants
benign_count=$(awk '{print}' no_header | grep -o "CLNSIG=[^;]*" | grep -E "Benign|Benign/Likely_benign|Likely_benign" | wc -l)
#vus variants
vus_count=$(awk '{print}' no_header | grep -o "CLNSIG=[^;]*" | grep "Uncertain_significance" | wc -l)
#other variants
other_count=$(awk '{print}' no_header | grep -o "CLNSIG=[^;]*" | grep "other" | wc -l)
#lines where pathogenic variants occur
lines=$(awk '{print}' no_header | grep -n "CLNSIG=Pathogenic" | cut -d ":" -f 1)

echo ""
echo "========================================================================================================================================================================================================="
echo ""
echo ""
echo "  [][] VariGene - Relatório Clínico [][]"
echo ""
echo ""
echo ""
echo "    [] Informação sobre a amostra []"
echo ""
echo ""
echo -e "        Amostra:""\t"$(ls *_annotated.vcf | grep -o "^[^_]*_[^_]*")
echo ""
echo "  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  "
echo ""
echo "    [] Sumário sobre as variantes []"
echo ""
echo ""
echo -e "        Variantes detectadas:""\t""\t""\t"$(echo $all_count)
echo -e "        Variantes identificadas:""\t""\t"$(echo $known_count)" ("$(echo $(echo "$(echo "scale=3;$known_count / $all_count" | bc) * 100" | bc | sed 's/..$//')"%")")"
echo "        |"
echo -e "        |_Variantes patogénicas:""\t""\t"$(echo $patho_count)
echo -e "        |_Variantes benignas:""\t""\t""\t"$(echo $benign_count)
echo -e "        |_Variantes de significado incerto:""\t"$(echo $vus_count)
echo -e "        |_Variantes (outras):""\t""\t""\t"$(echo $other_count)
echo ""
echo ""
echo "  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  "
echo ""
echo "    [] Variantes patogénicas []"
echo ""
echo ""
echo -e "        Cromossoma""\t""Gene""\t""\t""Variante""\t""\t""Tipo""\t""\t""\t""\t""\t""Doenças associadas"
for i in $lines;do CLNDN_i=$(awk 'NR=='$i' {print}' no_header | grep -o "CLNDN=[^;]*" | cut -d "=" -f 2);if echo $CLNDN_i | grep -qE "\|not_provided|\|not_specified";then if echo $CLNDN_i | grep -q "|not_specified";then CLNDN_i=$(echo $CLNDN_i | sed 's/|not_specified//');fi;if echo $CLNDN_i | grep -q "|not_provided";then CLNDN_i=$(echo $CLNDN_i | sed 's/|not_provided//');fi;fi;echo -e "        "$(awk 'NR=='$i' {print$1}' no_header | sed 's/^...//')"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "GENEINFO=[^;]*" | cut -d "=" -f 2 | cut -d ":" -f 1)"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "CLNHGVS=[^;]*" | cut -d "=" -f 2 | cut -d ":" -f 2)"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "CLNVC=[^;]*" | cut -d "=" -f 2)"\t""\t"$(echo $CLNDN_i);done
echo ""
echo ""
echo "========================================================================================================================================================================================================="
echo ""

rm no_header

#for i in $lines;do CLNDN_i=$(awk 'NR=='$i' {print}' no_header | grep -o "CLNDN=[^;]*" | cut -d "=" -f 2);if echo $CLNDN_i | grep -qE "\|not_provided|\|not_specified";then if echo $CLNDN_i | grep -q "|not_specified";then CLNDN_i=$(echo $CLNDN_i | sed 's/|not_specified//');fi;if echo $CLNDN_i | grep -q "|not_provided";then CLNDN_i=$(echo $CLNDN_i | sed 's/|not_provided//');fi;fi;echo -e "      "$(awk 'NR=='$i' {print$1}' no_header | sed 's/^...//')"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "GENEINFO=[^;]*" | cut -d "=" -f 2 | cut -d ":" -f 1)"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "CLNHGVS=[^;]*" | cut -d "=" -f 2 | cut -d ":" -f 2)"\t""\t"$(awk 'NR=='$i' {print}' no_header | grep -o "CLNVC=[^;]*" | cut -d "=" -f 2)"\t""\t"$(echo $CLNDN_i);done
