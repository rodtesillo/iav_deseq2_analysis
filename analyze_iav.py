# Importar librerias
import os
import argparse


# Crear funcion para filtrar significancia
def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    """Evalúa si un gen cumple los criterios de significancia."""

    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True

    return False
