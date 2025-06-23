#!/bin/bash

dia_file="snps_dia_significativos.txt"
ad_file="snps_ad_significativos.txt"
output_file="interseccao_variantes_significativas_ad_dia.txt"

# Processamento dos arquivos
awk '
BEGIN {OFS="\t"}
NR==1 {next} 
{
    # Processamento especial para o arquivo de diabetes
    chr = $1
    bp = $3
    
    # Extrai o primeiro rsID quando existir
    split($2, parts, ";")
    rsid = ""
    for (i in parts) {
        if (parts[i] ~ /^rs[0-9]+/) {
            rsid = parts[i]
            break
        }
    }
    
    # Se não encontrou rsID, usa o formato chr:pos
    if (rsid == "") {
        split($2, parts, ":")
        if (parts[1] == chr && parts[2] == bp) {
            rsid = $2
        }
    }
    
    print chr, bp, rsid, $0
}' "$dia_file" > dia_processed.tmp

awk '
BEGIN {OFS="\t"}
NR==1 {next}
{
    # Processamento para o arquivo de Alzheimer
    split($2, parts, ":")
    chr = parts[1]
    bp = parts[2]
    rsid = ($2 ~ /^rs[0-9]/) ? $2 : ""
    
    print chr, bp, rsid, $0
}' "$ad_file" > ad_processed.tmp

# Encontrando interseções por CHR:BP
join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11 \
     <(awk -F'\t' '{print $1":"$2 "\t" $0}' dia_processed.tmp | sort -k1,1) \
     <(awk -F'\t' '{print $1":"$2 "\t" $0}' ad_processed.tmp | sort -k1,1) > joined.tmp

# Formatando a saída final
echo -e "CHR\tBP\tSNP_DIABETES\tA1_DIABETES\tF_A_DIABETES\tF_U_DIABETES\tA2_DIABETES\tP_DIABETES\tOR_DIABETES\tSNP_ALZHEIMER\tA1_ALZHEIMER\tF_A_ALZHEIMER\tF_U_ALZHEIMER\tA2_ALZHEIMER\tP_ALZHEIMER\tOR_ALZHEIMER" > "$output_file"

awk -F'\t' '{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
    $2, $3, $8, $10, $11, $12, $13, $15, $16, $17, $19, $20, $21, $22, $24, $25
}' joined.tmp >> "$output_file"

# Verificação final
if [ $(wc -l < "$output_file") -eq 1 ]; then
    echo "ATENÇÃO: Nenhuma variante em comum foi encontrada."
    echo "Possíveis causas:"
    echo "1. Os arquivos realmente não compartilham variantes significativas"
    echo "2. Os formatos dos SNPs são muito diferentes"
    echo ""
    echo "Vamos verificar manualmente 10 exemplos do arquivo Alzheimer no arquivo Diabetes:"
    
    head -n 10 ad_processed.tmp | while read -r line; do
        chr=$(echo "$line" | cut -f1)
        bp=$(echo "$line" | cut -f2)
        echo "Buscando CHR$chr:$bp no arquivo Diabetes..."
        grep -c "$chr[[:space:]]$bp" "$dia_file"
    done
fi

# Limpeza
rm dia_processed.tmp ad_processed.tmp joined.tmp 2>/dev/null
