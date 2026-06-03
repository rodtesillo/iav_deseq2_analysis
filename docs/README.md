# IAV DESeq2 Analysis

Programa en Python para identificar genes diferencialmente expresados a partir de resultados de DESeq2 durante infección por Influenza A Virus.

```
iav_deseq2_analysis/
├── data/
├── results/
├── docs/
├── analyze_iav.py
└── README.md
```

# Como correr el programa

```python
python analyze_iav.py data/iav_deseq2_results.tsv results/iav_significant_genes.tsv
```

# Usos de threeholds opcionales

python analyze_iav.py data/iav_deseq2_results.tsv results/iav_significant_genes.tsv --lfc_threshold 2.0 --padj_threshold 0.01 (En caso de no asiganar valores el programa usara 1 y 0.05 por defecto)

# Ejemplos

1. Usar valores por defecto (LFC ≥ 1.0, padj < 0.05)
bashpython analyze_iav.py data/iav_deseq2_results.tsv results/resultados.tsv
2. Usar LFC más estricto
bashpython analyze_iav.py data/iav_deseq2_results.tsv results/resultados.tsv --lfc_threshold 1.5
3. Personalizar ambos umbrales
bashpython analyze_iav.py data/iav_deseq2_results.tsv results/resultados.tsv --lfc_threshold 2.0 --padj_threshold 0.01
4. Ver ayuda
bashpython analyze_iav.py --help

# Ejemplo de entrada:

```
gene_id	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
MX1	1250.5	4.2	0.15	28.0	1e-10	0.0001
GENE1	890.3	-3.0	0.12	-25.0	1e-9	0.001
MYC	450.2	0.8	0.18	4.4	1e-5	0.15
```

# Ejemplo de salida:

```
gene	log2FoldChange	padj	status
MX1	4.2	0.0001	upregulated
GENE1	-3.0	0.001	downregulated
TP53	-1.8	0.0001	downregulated
```

# Interpretación de log2FoldChange:

log2FC > 0: Gen sobreexpresado (upregulated) en la condición de tratamiento
log2FC < 0: Gen subexpresado (downregulated) en la condición de tratamiento
log2FC = 2.0: La expresión se multiplica por 4x (2^2)
log2FC = 1.0: La expresión se multiplica por 2x (2^1)
log2FC = -1.0: La expresión se divide por 2x (2^-1)

