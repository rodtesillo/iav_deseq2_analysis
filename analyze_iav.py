# Importar librerias
import os
import argparse


# Crear funcion para filtrar significancia
def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    """Evalúa si un gen cumple los criterios de significancia."""

    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True

    return False


# Crear función para clasificar genes en upregulated y downregulated
def classify_gene(log2_fold_change):
    """Evalúa si un gen aumenta o disminuye su expresión."""

    if log2_fold_change > 0:
        return "upregulated"
    elif log2_fold_change < 0:
        return "downregulated"
    else:
        return "Sin cambios"
