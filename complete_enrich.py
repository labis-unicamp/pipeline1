#!/usr/bin/env python3
"""
Complete Biological Pathway Enrichment Analysis
This script performs enrichment analysis across multiple databases
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

try:
    import gseapy as gp
    from gseapy import enrichr
except ImportError:
    print("Installing gseapy...")
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "gseapy"])
    import gseapy as gp
    from gseapy import enrichr

def clean_gene_list(gene_string):
    """Cleans and processes the gene list"""
    genes = gene_string.replace(',', ' ').split()
    unique_genes = []
    seen = set()
    
    for gene in genes:
        if ',' in gene:
            for subgene in gene.split(','):
                subgene = subgene.strip()
                if subgene and subgene not in seen:
                    unique_genes.append(subgene)
                    seen.add(subgene)
        elif ';' in gene:
            subgene = gene.split(';')[0].strip()
            if subgene and subgene not in seen:
                unique_genes.append(subgene)
                seen.add(subgene)
        else:
            gene = gene.strip()
            if gene and gene not in seen:
                unique_genes.append(gene)
                seen.add(gene)
    
    return unique_genes

def perform_comprehensive_enrichment(gene_list):
    """
    Performs enrichment across multiple databases
    """
    results = {}
    databases = [
        'KEGG_2021_Human',
        'WikiPathway_2023_Human',
        'Reactome_2022',
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'BioPlanet_2019',
        'MSigDB_Hallmark_2020',
        'Panther_2016',
        'HumanCyc_2016'
    ]
    
    for db in databases:
        try:
            print(f"Analyzing database: {db}...")
            enr = enrichr(gene_list=gene_list,
                         gene_sets=db,
                         organism='Human',
                         cutoff=0.1)  # More permissive p-value
            
            if not enr.results.empty:
                results[db] = enr.results
                print(f"  - {len(enr.results)} pathways found")
            else:
                print(f"  - No pathways found")
                
        except Exception as e:
            print(f"  - Error with {db}: {e}")
            continue
    
    return results

def categorize_pathways(pathway_name):
    """
    Categorizes pathways based on their content
    """
    name_lower = pathway_name.lower()
    
    # Main categories
    categories = {
        'Metabolism': ['metabolic', 'glycolysis', 'gluconeogenesis', 'fatty acid', 'lipid', 
                       'oxidative phosphorylation', 'tca cycle', 'pentose phosphate'],
        'Cell Signaling': ['signaling', 'pathway', 'pi3k', 'akt', 'mapk', 'mtor', 'nf-kb', 
                               'wnt', 'notch', 'jak-stat', 'tgfb', 'camp', 'calcium'],
        'Diseases': ['cancer', 'diabetes', 'alzheimer', 'parkinson', 'cardiovascular', 'viral',
                   'infectious', 'autoimmune', 'inflammatory'],
        'Cellular Processes': ['cell cycle', 'apoptosis', 'autophagy', 'dna repair', 
                               'replication', 'transcription', 'translation'],
        'Immune System': ['immune', 'immunity', 't cell', 'b cell', 'interferon', 
                               'complement', 'antigen', 'cytokine'],
        'Neuroscience': ['neuron', 'synaptic', 'neurotransmitter', 'dopamine', 'serotonin',
                        'gaba', 'glutamate', 'long-term potentiation'],
        'Development': ['development', 'differentiation', 'morphogenesis', 'embryonic',
                           'organogenesis', 'stem cell'],
        'Transport': ['transport', 'ion channel', 'membrane', 'vesicle', 'endocytosis',
                      'exocytosis']
    }
    
    for category, keywords in categories.items():
        if any(keyword in name_lower for keyword in keywords):
            return category
    
    return 'Other'

def get_color_map(categories):
    """
    Generates a color map for categories
    """
    # More robust color palette
    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
        '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
    ]
    
    color_map = {}
    for i, category in enumerate(categories):
        color_map[category] = colors[i % len(colors)]
    
    return color_map

def create_pathway_visualization(all_results):
    """
    Creates visualization of the most relevant pathways found
    """
    # Combine all results
    all_pathways = []
    for db_name, results in all_results.items():
        for _, row in results.iterrows():
            all_pathways.append({
                'Term': row['Term'],
                'P-value': row['P-value'],
                'Database': db_name,
                'Genes': row['Overlap'],
                'Combined_Score': row.get('Combined Score', 0),
                'Category': categorize_pathways(row['Term'])
            })
    
    if not all_pathways:
        print("No pathways found for visualization")
        return
    
    # Convert to DataFrame and sort by p-value
    pathway_df = pd.DataFrame(all_pathways)
    pathway_df = pathway_df.sort_values('P-value').head(30)  # Top 30
    
    # Style configuration
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
    
    # 1. Bar plot of most significant pathways
    top_pathways = pathway_df.head(15)
    y_pos = np.arange(len(top_pathways))
    
    # Get all unique categories and create color map
    all_categories = pathway_df['Category'].unique()
    color_map = get_color_map(all_categories)
    colors = [color_map[cat] for cat in top_pathways['Category']]
    
    bars = ax1.barh(y_pos, -np.log10(top_pathways['P-value']), color=colors, alpha=0.7)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([f"{term[:45]}..." if len(term) > 45 else term 
                         for term in top_pathways['Term']], fontsize=10)
    ax1.set_xlabel('-log10(P-value)', fontsize=12)
    ax1.set_title('Top 15 Most Significant Pathways', fontsize=14)
    ax1.grid(axis='x', alpha=0.3)
    
    # Add values to bars
    for i, (bar, pval) in enumerate(zip(bars, top_pathways['P-value'])):
        width = bar.get_width()
        ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                f'{pval:.1e}', ha='left', va='center', fontsize=8)
    
    # Add category legend
    legend_elements = [plt.Rectangle((0,0),1,1, color=color_map[cat], label=cat) 
                      for cat in all_categories]
    ax1.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    # 2. Distribution by database
    db_counts = pathway_df['Database'].value_counts()
    ax2.pie(db_counts.values, labels=db_counts.index, autopct='%1.1f%%', startangle=90)
    ax2.set_title('Pathway Distribution by top5 Database', fontsize=14)
                 
    # 3. Scatter plot: P-value vs Combined Score
    scatter_colors = [color_map[cat] for cat in pathway_df['Category']]
    
    scatter = ax3.scatter(pathway_df['Combined_Score'], -np.log10(pathway_df['P-value']), 
               c=scatter_colors, alpha=0.7, s=60)
    ax3.set_xlabel('Combined Score', fontsize=12)
    ax3.set_ylabel('-log10(P-value)', fontsize=12)
    ax3.set_title('Significance vs Combined Score', fontsize=14)
    ax3.grid(True, alpha=0.3)
    
    # Add p=0.05 cutoff line
    ax3.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
    ax3.legend()
    
    # 4. Distribution by category
    category_counts = pathway_df['Category'].value_counts()
    colors_cat = [color_map[cat] for cat in category_counts.index]
    
    bars_cat = ax4.bar(range(len(category_counts)), category_counts.values, color=colors_cat, alpha=0.7)
    ax4.set_xticks(range(len(category_counts)))
    ax4.set_xticklabels(category_counts.index, rotation=45, ha='right')
    ax4.set_ylabel('Number of Pathways', fontsize=12)
    ax4.set_title('Distribution by Biological Category', fontsize=14)
    ax4.grid(axis='y', alpha=0.3)
    
    # Add values to bars
    for i, v in enumerate(category_counts.values):
        ax4.text(i, v + 0.1, str(v), ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('comprehensive_pathway_analysis.png', dpi=300, bbox_inches='tight')
    print("Pathway visualization saved at: comprehensive_pathway_analysis.png")

def create_detailed_pathway_report(all_results):
    """
    Creates a detailed pathway report in image format
    """
    # Combine and organize all pathways
    all_pathways = []
    for db_name, results in all_results.items():
        for _, row in results.iterrows():
            # Extract overlapping genes
            overlap_info = row['Overlap'].split('/')
            gene_count = int(overlap_info[0])
            total_genes = int(overlap_info[1])
            
            all_pathways.append({
                'Pathway': row['Term'][:60] + '...' if len(row['Term']) > 60 else row['Term'],
                'P-value': row['P-value'],
                'Gene_Count': gene_count,
                'Total_Genes': total_genes,
                'Gene_Ratio': gene_count / total_genes,
                'Database': db_name.replace('_', ' '),
                'Combined_Score': row.get('Combined Score', 0),
                'Category': categorize_pathways(row['Term'])
            })
    
    if not all_pathways:
        return
    
    # Sort by p-value and get top 20
    pathway_df = pd.DataFrame(all_pathways).sort_values('P-value').head(20)
    
    # Get color map for all categories
    all_categories = pathway_df['Category'].unique()
    color_map = get_color_map(all_categories)
    
    # Create figure for detailed report
    fig, ax = plt.subplots(figsize=(18, 14))
    ax.axis('tight')
    ax.axis('off')
    
    # Create detailed table
    table_data = []
    headers = ['Rank', 'Pathway', 'P-value', 'Genes\nOverlap', 'Gene\nRatio', 'Database', 'Category']
    
    for i, (_, row) in enumerate(pathway_df.iterrows()):
        table_data.append([
            f"{i+1}",
            row['Pathway'],
            f"{row['P-value']:.2e}",
            f"{row['Gene_Count']}/{row['Total_Genes']}",
            f"{row['Gene_Ratio']:.3f}",
            row['Database'],
            row['Category']
        ])
    
    # Create table
    table = ax.table(cellText=table_data, colLabels=headers, 
                    cellLoc='left', loc='center', bbox=[0, 0, 1, 1])
    
    # Format table
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    
    # Color header
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#2E86AB')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Color alternating rows and by category
    for i in range(1, len(table_data) + 1):
        for j in range(len(headers)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f8f8f8')
            
            # Color by category in last column
            if j == 6:  # category column
                category = table_data[i-1][6]
                table[(i, j)].set_facecolor(color_map.get(category, '#F0F0F0'))
            
            # Highlight significant p-values
            if j == 2:  # p-value column
                try:
                    pval_str = table_data[i-1][2].replace('e-', 'e-')
                    pval = float(pval_str)
                    if pval < 0.001:
                        table[(i, j)].set_facecolor('#ffcccc')
                    elif pval < 0.01:
                        table[(i, j)].set_facecolor('#ffe6cc')
                    elif pval < 0.05:
                        table[(i, j)].set_facecolor('#ffffcc')
                except:
                    pass
    
    plt.title('Detailed Report of Most Significant Pathways', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Add legend
    legend_text = ("Color legend:\n"
                  "🔴 Red: p < 0.001 (very significant)\n"
                  "🟠 Orange: p < 0.01 (significant)\n" 
                  "🟡 Yellow: p < 0.05 (moderately significant)\n\n"
                  "Colors by biological category")
    
    plt.figtext(0.02, 0.02, legend_text, fontsize=9, 
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
    
    plt.savefig('detailed_pathway_report.png', dpi=300, bbox_inches='tight')
    print("Detailed report saved at: detailed_pathway_report.png")

def create_database_comparison(all_results):
    """
    Creates comparative visualization between databases
    """
    if not all_results:
        return
    
    # Prepare data for comparison
    comparison_data = []
    for db_name, results in all_results.items():
        if not results.empty:
            # Get top 5 pathways per database
            top_pathways = results.head(5)
            for _, row in top_pathways.iterrows():
                comparison_data.append({
                    'Database': db_name,
                    'Pathway': row['Term'],
                    'P-value': row['P-value'],
                    'Combined_Score': row.get('Combined Score', 0)
                })
    
    if not comparison_data:
        return
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 1. Number of significant pathways per database
    db_counts = {db: len(results) for db, results in all_results.items() if not results.empty}
    ax1.bar(range(len(db_counts)), list(db_counts.values()), color='skyblue', alpha=0.7)
    ax1.set_xticks(range(len(db_counts)))
    ax1.set_xticklabels(list(db_counts.keys()), rotation=45, ha='right')
    ax1.set_ylabel('Number of Pathways Found')
    ax1.set_title('Significant Pathways by Database')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add values to bars
    for i, v in enumerate(db_counts.values()):
        ax1.text(i, v + 0.5, str(v), ha='center', va='bottom')
    
    # 2. Best p-value per database
    best_pvalues = {}
    for db, results in all_results.items():
        if not results.empty:
            best_pvalues[db] = results['P-value'].min()
    
    if best_pvalues:
        dbs = list(best_pvalues.keys())
        pvals = [-np.log10(best_pvalues[db]) for db in dbs]
        
        bars = ax2.bar(range(len(dbs)), pvals, color='lightcoral', alpha=0.7)
        ax2.set_xticks(range(len(dbs)))
        ax2.set_xticklabels(dbs, rotation=45, ha='right')
        ax2.set_ylabel('-log10(P-value)')
        ax2.set_title('Best P-value by Database')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add values to bars
        for i, (bar, pval) in enumerate(zip(bars, list(best_pvalues.values()))):
            ax2.text(i, bar.get_height() + 0.1, f'{pval:.1e}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('database_comparison.png', dpi=300, bbox_inches='tight')
    print("Database comparison saved at: database_comparison.png")

def main():
    """Main function for complete analysis"""
    print("=== COMPLETE PATHWAY ENRICHMENT ANALYSIS ===\n")
    
    # Provided gene list
    gene_string = """MTCL1
AGBL4
C1orf74
GIGYF1
PHF21A
OR10D4P
STK31
NOC2L
SSU72
PLCH2
PRDM16
LRRC47
RNF207
RERE
PEX14
VPS13D
PRAMEF7
FBLIM1
UQCRHL
UBR4
UBR4
RNU4-28P
MICOS10
EIF4G3
ELOA
FUCA1
LINC02800
NIPAL3
RHD
CRYBG2
RPS6KA1
GPATCH3
THEMIS2
TMEM39B
THRAP3
SH3D21
INPP5B
UTP11
LOC105378654
LOC105378669
SCMH1-DT
HIVEP3
MAST2
CYP4Z2P
AGBL4
GLIS1
GLIS1
TCEANC2
USP24
LINC01755
LINC01755
DOCK7
DOCK7-DT
ROR1
ROR1
PDE4B
PDE4B
GNG12
ST6GALNAC3
ZZZ3
GIPC2
LPAR3
SELENOF
SELENOF
TGFBR3
TGFBR3
RNPC3-DT
RNPC3
CELSR2
SORT1
ATXN7L2
LINC01525
HYDIN2
FALEC
CTSS
GABPB2
GABPB2
PIP5K1A
CCDST
CCDST
LCE1F
SPRR4
TDRD10
UBE2Q1
HDGF
INSRR
OR10AA1P
IFI16
IFI16
DUSP23
Y_RNA
ATF6
RGS5
RXRG
POU2F1
DCAF6
F5
F5
TNN
STX6
CACNA1E
CFHR4
ASCL5
PIK3C2B
LOC105371692
NFASC
IL24
LINC02942
LINC02942
CD46
HHAT
HHAT
LPGAT1-AS1
LOC105372916
LOC105372916
RAB3GAP2
BROX
RNU6-1319P
MRPL55
MRPL55
H2AC25
RNF187
RNF187
RNA5S3
C1orf198
DISC1
ARID4B
EDARADD
ACTN2
SMYD3
OR6R1P
OR2L6P
OR2L3
OR14I1
CALML3-AS1
IL2RA
CELF2
SEPHS1
SEPHS1
APBB1IP
CCDC7
PARD3
ZNF37A
MARCHF8
MAPK8
JMJD1C
LINC02671
PALD1
ADAMTS14
PCBD1
MICU1
MICU1
KAT6B
KAT6B
LRMDA
LRMDA
NRG3
LIPK
HECTD2
HECTD2
CYP2C19
SORBS1
SORBS1
SORBS1
PIK3AP1
RRP12
R3HCC1L
HPSE2
HPSE2
CNNM1
CNNM1
CHUK
ARMH3
PCGF6
LINC01435
CCDC186
ATRNL1
CPXM2
CPXM2
PHRF1
IRF7
EPS8L2
MUC2
MUC5AC
LSP1
TRPM5
ZNF195
MMP26
MMP26
MMP26
OR52S1P
HBE1
HBE1
OR52N5
OR52N2
OR52N2
TRIM5
OR56A4
SCUBE2
SBF2-AS1
AMPD3
LYVE1
ZBED5
ABCC8
USH1C
LDHC
PTPN5
PTPN5
MRGPRX11P
ANO3
SLC5A12
WT1-AS
WT1-AS
TRAF6
ALX4
AMBRA1
MADD
CELF1
OR4A11P
OOSP1
MS4A2
SCGB1D4
ASRGL1
EEF1G
MARK2
MACROD1
MACROD1
CAPN1
FRMD8
PCNX3
PCNX3
PACS1
CCS
KMT5B
PPP6R3
SHANK2
ZNF705EP
CLPB
ARRB1
PAK1
GAB2
PIWIL4
CEP57
NCAM1
RNF214
TMPRSS4
TTC36
DDX6
JHY
PKNOX2
RPUSD4
PPP1R10P1
PPP1R10P1
PPP1R10P1
PPP1R10P1
KDM5A
GALNT8
TAPBPL
ACRBP
MLF2
KLRG1
KLRF1
LINC02446
TAS2R14
TAS2R42
EPS8
EPS8
PIK3C2G
SLCO1B7
SLCO1A2
SLCO1A2
LDHB
TMTC1
RNU6-494P
PRICKLE1
PRICKLE1
OR5BJ1P
SPMIP11
TUBA1B-AS1
TMBIM6
SLC4A8
KRT79
CSAD
CSAD
ZNF740
PCBP2
ATP5MC2
ATP5MC2
ARHGEF25
SRGAP1
MSRB3-AS1
OTOGL
MGAT4C
TMCC3
NDUFA12
ANO4
C12orf75
SVOP
RILPL1
ZNF664
TMEM132B
TMEM132B
GALNT9
POLE
ZDHHC20
RNF17
RPL21
PAN3
NBEA
LHFPL6
LHFPL6
LINC00428
RB1
LINC00434
OBI1-AS1
LINC01080
ABCC4
GGACT
MCF2L
TMEM255B
TMEM255B
OSGEP
PPP1R3E
IL25
STXBP6
HECTD1
NUBPL
PRORP
SLC25A21
TRIM9
DDHD1-DT
Y_RNA
ZBTB25
RPPH1-2P
RPL36AP2
GALNT16
NUMB
MIDEAS
FBLN5
IFI27L1
VRK1
MEG8
LOC105370673
LOC105370673
PACS2
TEX22
TEX22
TEX22
TEX22
SNRPN
LINC02346
GABRG3
OCA2
HERC2
GPR176
RPUSD2
KNL1
PLA2G4E
PLA2G4F
STARD9
UBR1
STRCP1
PDIA3
FRMD5
GATM
AFG2B
SQOR
SQOR
FBN1
FGF7
SLC27A2
UNC13C
FAM81A
FAM81A
IGDCC3
CCDC33
PTPN9
NRG4
UBE2Q2P1
AKAP13
AKAP13
AGBL1
AGBL1
AGBL1
MIR9-3HG
MIR9-3HG
MIR9-3HG
AP3S2
SEMA4B
GDPGP1
GABARAPL3
BLM
CRAT37
LOC105371026
LINC02348
TARS3
OR4F14P
AXIN1
FBXL16
FBXL16
ANTKMT
TELO2
TELO2
TMEM204
IFT140
MEIOB
NHERF2
ZSCAN10
NAA60
ROGDI
PPL
EMP2
TEKT5
NOMO3
NOMO3
TMC7
ACSM1
ACSM1
HS3ST2
NSMCE1
SULT1A1
SULT1A1
ATP2A1
NKD1
RBL2
MT1F
MT1H
MT1H
MT1H
MIR138-2
CDH11
SLC12A4
PRMT7
CDH3
SF3B3
PHLPP2
PHLPP2
PMFBP1
WDR59
BCAR1
TERF2IP
ADAMTS18
PKD1L2
PKD1L2
PKD1L2
NECAB2
SLC38A8
KCNG4
USP10
USP10
LINC02139
JPH3
CBFA2T3
CDH15
FAM157C
DOC2B
VPS53
SMG6
SMG6
OR1D2
OR1E3
OR1A1
ASPA
SHPK
CTNS
CTNS
ZZEF1
UBE2G1
SPNS3
ALOX15
MINK1
ZFP3
RABEP1
DERL2
PITPNM3
ZBTB4
VAMP2
CCDC42
MYH2
MYH2
SHISA6
COX10-DT
KRT17P1
RAI1
DRC3
GID4
LLGL1
SPECC1
LYRM9
RSKR
MYO18A
TAOK1
MYO1D
ASIC2
CCL16
LHX1-DT
AATF
ACACA
PLXDC1
PLXDC1
PLXDC1
ARL5C
WIPF2
KRTAP9-2
KRTAP9-4
NKIRAS2
STAT5B
RUNDC1
LINC01976
G6PC3
G6PC3
HDAC5
SLC4A1
NSF
CBX1
SKAP1
ANKFN1
ANKFN1
ANKFN1
LOC101927510
CLTC
VMP1
VMP1
MED13
KCNH6
TEX2
PRKCA
HELZ
SLC39A11
SLC39A11
SLC39A11
SDK2
LINC02074
DNAI2
SRP68
RHBDF2
ST6GALNAC2
LINC01993
SOCS3
TIMP2
RBFOX3
TBC1D16
RPTOR
RPTOR
GREB1L
IMPACT
ASXL3
ASXL3
ZBTB7C
TCF4
CCBE1
CCBE1
CDH20
CDH19
DSEL-AS1
GZMM
ABCA7
SF3A2
TMPRSS9
ZNF77
TLE2
SH3GL1
FUT6
RANBP3
CRB3
C3
FBN3
CERS4
NFILZ
MBD3L1
MUC16
MUC16
MUC16
MUC16
MUC16
OR1M4P
OR7D1P
OLFM2
PPAN
PDE4A
SLC44A2
SLC44A2
SLC44A2
ILF3
DOCK6
ECSIT
ZNF700
ZNF799
MTATP6P27
ZNF490
HOOK2
CC2D1A
PODNL1
OR7A2P
NOTCH3
BRD4
OR10H2
CYP4F24P
UCA1-AS1
UCA1-AS1
CYP4F11
MYO9B
ABHD8
UNC13A
FCHO1
ELL
FKBP8
GATAD2A
GMIP
ZNF506
ZNF682
ZNF100
ZNF726
PEPD
LINC02965
PROSER3
ZNF529-AS1
ZNF875
SIPA1L3
SIPA1L3
SIPA1L3
DPF1
FBXO17
FBXO27
RPS29P25
PSMC4
SPTBN4
TGFB1
PSG3
IGSF23
PVR
CBLC
APOC1P1
FBXO46
FOXA3
DHDH
SNRNP70
TRPM4
DKKL1
TBC1D17
SIGLEC16
VRK3
ZNF473CR
MYH14
MYH14
KLK3
CD33
SIGLEC24P
LIM2
ZNF615
ZNF615
ZNF578
ZNF845
MIR512-2
MBOAT7
LILRB1
KIR2DL3
PPP1R12C
TNNT1
NLRP5
ZNF587B
ZNF587B
ZNF587
ZNF587
C19orf18
ZNF329
ZNF329
LPIN1
SDC1
DTNB
OTOF
OTOF
NRBP1
IFT172
RBKS
CLIP4
ALK
ALK
ALK
CDC42EP3-AS1
LHCGR
NRXN1-DT
EML6
EHBP1
WDPCP
ARHGAP25
ZNF638
ZNF638
ALMS1
STAMBP
SEMA4F
KCMF1
TCF7L1
IGKV2-26
VWA3B
MRPL30
SH3RF3
SEPTIN10
LINC02966
ACTR3-AS1
LINC01956
LOC107985817
EPB41L5
TSN
TSN
FMNL2
LOC105373691
UPP2
TANC1
PLA2R1
LOC105373718
KLHL23
UBR3
MYO3B
SLC25A12
FKBP7
ZNF385B
ZNF385B
ZNF385B
PDE1A
ITGAV
CAVIN2-AS1
LINC01923
FLACC1
BMPR2
PTH2R
MARCHF4
STK36
STK11IP
NYAP2
DAW1
SPHKAP
LINC01807
ARMC9
DIS3L2
DIS3L2
DGKD
LOC105373936
SH3BP4
PER2
PER2
LINC01937
HDAC4
UICLM
FARP2
FARP2
ING5
DEFB126
DEFB132
SNPH
SIRPD
SIRPG
PTPRA
DNAAF9
DNAAF9
ATRN
MAVS
PLCB1
ANKEF1
SLX4IP
LINC00261
APMAP
PYGB
SAMHD1
SRC
PTPRT
HNF4A
SPINT4
LOC105376989
CDH4
CDH4
MTG2
LAMA5
LINC02970
TCFL5
C1QBPP1
ADAMTS5
LINC01695
LOC105372772
HUNK
RUNX1
B3GALT5
DSCAM
DSCAM
WDR4
WDR4
CRYAA
AIRE
SLC19A1
PCBP3
PSMA6P3
COL6A2
LSS
PPP1R26P5
TOP3B
BCR
CABIN1
LRRC75B
GRK3
GRK3
MIATNB
MTMR3
OSBP2
OSBP2
RNA5SP496
RFPL2
SYN3
MFNG
XRCC6
CCDC134
NAGA
EFCAB6
EFCAB6
NUP50-DT
KIAA0930
UPK3A
ATXN10
GRAMD4
TBC1D22A
TBC1D22A
TBC1D22A
PPP6R2
TYMP
SHANK3
SNRPA1P2
CNTN6
IL5RA
IL5RA
ITPR1
GRM7
RAD18
CIDEC
NUP210
ANKRD28
DPH3
DPH3
KCNH8
SLC4A7
STAC
GOLGA4
DLEC1
SCN10A
SCN11A
ULK4
ULK4
ULK4
ULK4
ANO10
ZDHHC3
RTP3
KLHL18
MAP4
PLXNB1
QRICH1
MIR135A1
TASOR
PRICKLE2
EPHA6
CD47
ABHD10
CCDC80
CFAP44
ZDHHC23
IGSF11
B4GALT4
ARHGAP31
STXBP5L
STXBP5L
STXBP5L
HCLS1
CASR
KPNA1
PARP15
MUC13
ALDH1L1-AS2
MGLL
EEFSEC
HMCES
LINC02006
LINC02006
LINC02006
RSRC1
RSRC1
IL12A
LINC01322
FNDC3B
MFN1
LINC01994
MCF2L2
TP63
PPP1R2
SDHAP1
NRROS
FYTTD1
IDUA
SPON2
SPON2
FAM53A
TACC3
RGS12
DOK7
OTOP1
OTOP1
OTOP1
EVC
EVC
TBC1D14
SH3TC1
GPR78
ZNF518B
HS3ST1
CC2D2A
FBXL5
FGFBP2
CLRN2
LAP3
LCORL
LGI2
SLC34A2
RBPJ
PCDH7
TLR6
LIMCH1
SLC30A9
KCTD8
GABRA2
UGT2B7
CSN2
AREG
PCAT4
PRKG2
PRKG2
LINC02994
AFF1
IBSP
MANBA
LINC02503
EGF
PITX2
LINC01438
AFG2A
SCLT1
ARHGAP10
DCHS2
NPY2R
TMEM144
TMEM192
NEK1
SORBS2
SDHA
AHRR
AHRR
SLC9A3
SLC9A3
SLC9A3
MTRR
DNAH5
MYO10
AMACR
C7
MROH2B
PLCXD3
LINC02105
ANKRD55
LINC02122
SV2C
IQGAP2
MIR3660
RHOBTB3
STARD4-AS1
REEP5
CCDC112
MARCHF3
CCDC192
P4HA2
PFDN1
SLC35A4
PCDHGA1
PRELID2
STK32A
DPYSL3
SH3TC2
ARHGEF37
ARHGEF37
PPARGC1B
SGCD
SGCD
IL12B-AS1
TENM2
TENM2
TENM2
CPLX2
THOC3
UNC5A
RMND5B
ADAMTS2
IRF4
PXDC1
SIRT5
SLC17A1
H1-1
H2AC4
H2AC5P
H2AC6
H4C8
BTN2A1
H2BC15
GABBR1
HCG4
HCG4
HCG4
DDX39BP1
RN7SKP186
DDR1
RNU6-1133P
C6orf15
HCP5
TNXB
TNXB
TNXB
TSBP1
TAP2
BRD2
RPS18
ANKS1A
MAPK14
CPNE5
TREML2
CNPY3
PPP2R5D
MRPL14
SUPT3H
CYP39A1
MEP1A
CRISP2
DST
DST
ASCC3
CRYBG1
SNX3
WASF1
TRAF3IP2
HDAC2-AS2
DSE
KPNA5
SAMD3
LINC01013
UTRN
LATS1
ESR1
SYNE1
FNDC1
FNDC1
SLC22A3
MAP3K4
RPS6KA2
SMOC2
WDR27
LOC105375123
IQCE
TTYH3
RADIL
RADIL
TNRC18
TNRC18
CYTH3
RAC1
KDELR2
AGR2
ABCB5
HOTTIP
MINDY4
MINDY4
BBS9
GLI3
TNS3
ABCA13
STYXL1
GNAI1
PCLO
ELAPOR2
CASD1
ZSCAN21
STAG3
PILRA
PILRA
GIGYF1
POLR2J2
RBM28
METTL2B
STRIP2
CALD1
KDM7A
TCAF2
CUL1
ZNF775
LOC107986860
LOC107986860
LOC107986860
LOC124901774
GIMAP2
PRKAG2
CNPY1
MIR595
RPL21P76
DLGAP2
DLGAP2
CLN8
ARHGEF10
ARHGEF10
CSMD1
TNKS
TUSC3
MTUS1
MTUS1-DT
XPO7
PDLIM2
RHOBTB2
ADAM28
DOCK5
XKR4
LINC01592
CHMP4C
DPYS
ZFPM2
MRPL13
TRIB1AL
IQANK1
ARHGAP39
ZNF517
ZNF16
DOCK8
DMRT2
IL11RA
VCP
NPR2
GNAQ
ISCA1
AUH
PHF2
ERP44
LOC105376187
PALM2AKAP2
PRPF4
PAPPA
PAPPA
PAPPA
GAPVD1
ASS1
PRRC2B
AK8
VAV2
RXRA
UBAC1
NACC2
QSOX2
CCDC187
CCDC187
ABCA2
EHMT1
"""    

    gene_list = clean_gene_list(gene_string)
    print(f"Total unique genes: {len(gene_list)}\n")
    
    # Enrichment across multiple databases
    print(f"=== ENRICHMENT ACROSS MULTIPLE DATABASES ===")
    all_results = perform_comprehensive_enrichment(gene_list)
    
    # Create visualizations of found pathways
    if all_results:
        print(f"\n=== GENERATING VISUALIZATIONS ===")
        create_pathway_visualization(all_results)
        create_detailed_pathway_report(all_results)
        create_database_comparison(all_results)
    
    # Save complete results
    if all_results:
        print(f"\n=== SAVING COMPLETE RESULTS ===")
        with pd.ExcelWriter('comprehensive_pathway_results.xlsx', engine='openpyxl') as writer:
            # Pathways by database
            for db_name, results in all_results.items():
                if not results.empty:
                    sheet_name = db_name.replace('_', '')[:31]  # Excel character limit
                    results.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print("Complete results saved at: comprehensive_pathway_results.xlsx")
    
    # Analysis summary
    print(f"\n=== ANALYSIS SUMMARY ===")
    print(f"• Total genes analyzed: {len(gene_list)}")
    
    total_pathways = 0
    for db_name, results in all_results.items():
        if not results.empty:
            count = len(results)
            total_pathways += count
            print(f"• {db_name}: {count} pathways found")
    
    print(f"• Total enriched pathways: {total_pathways}")
    print(f"• Databases with results: {len([r for r in all_results.values() if not r.empty])}")

if __name__ == "__main__":
    main()

