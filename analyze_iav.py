# Importar librerias
import sys


# Crear funcion para filtrar significancia
def is_significant(log2_fold_change, padj, lfc_threshold=1.0, padj_threshold=0.05):
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


# Crear funcción para abrir el archivo, leerlo y devolver una lista de genes validos
def load_deseq2_results(file_path):
    valid_genes = []
    invalid_lines = []

    try:
        # Abrir archivo
        with open(file_path, "r") as file:
            # Leer línea por línea
            for line_number, line in enumerate(file, start=1):

                # Ignorar líneas vacías
                line = line.strip()
                if not line:
                    continue

                # Ignorar encabezado
                if line_number == 1:
                    continue

                # Separar columnas
                columns = line.split("\t")

                # Validar las 3 columnas de gene, log2fc y padj
                if len(columns) < 3:
                    invalid_lines.append((line_number, "Columnas insuficientes"))
                    continue

                try:
                    # Convertir valores numéricos
                    gene_id = columns[0].strip()
                    log2_fold_change = float(columns[2])
                    padj = float(columns[6])

                    # Ignorar lineas invalidas
                    if not (
                        isinstance(log2_fold_change, float) and isinstance(padj, float)
                    ):
                        raise ValueError("Valores no numéricos")

                    if padj < 0 or padj > 1:
                        raise ValueError("p-valor fuera de rango [0, 1]")

                    # Si todo es válido, agregar a la lista
                    valid_genes.append(
                        {
                            "gene_id": gene_id,
                            "log2_fold_change": log2_fold_change,
                            "padj": padj,
                        }
                    )

                # Ignorar líneas inválidas
                except (ValueError, IndexError) as e:
                    invalid_lines.append((line_number, str(e)))
                    continue

    except FileNotFoundError:
        print(f"Error: No se encontró el archivo '{file_path}'")
        return []

    except Exception as e:
        print(f"Error al leer el archivo: {e}")
        return []

    # Mostrar resumen
    print(f"Genes válidos cargados: {len(valid_genes)}")
    if invalid_lines:
        print(f"Líneas ignoradas: {len(invalid_lines)}")
        for line_num, reason in invalid_lines[:5]:  # Mostrar primeras 5
            print(f"   Línea {line_num}: {reason}")
        if len(invalid_lines) > 5:
            print(f"   ... y {len(invalid_lines) - 5} más")

    return valid_genes


# Crear función para crear una lista con los genes significativos
def filter_genes(genes, lfc_threshold=1.0, padj_threshold=0.05):
    """Mete los genes filtrados en una lista."""

    # Crear una lista vacía para los resultados filtrados
    filtered_results = []

    # Recorrer cada gen
    for gene in genes:
        gene_id = gene["gene_id"]
        log2_fold_change = gene["log2_fold_change"]
        padj = gene["padj"]

        # Aplicar is_significant()
        if is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):

            # Si es significativo, aplicar classify_gene()
            classification = classify_gene(log2_fold_change)

            # Guardar una tupla con el estado final
            filtered_results.append((gene_id, log2_fold_change, padj, classification))

    return filtered_results


# Función para escribir los genes filtrados
def write_results(filtered_genes, output_file):
    """Escribe los genes filtrados en un archivo de salida."""

    try:
        with open(output_file, "w") as file:
            # Escribir encabezado
            file.write("gene\tlog2FoldChange\tpadj\n")

            # Escribir cada gen filtrado
            for gene_id, log2_fold_change, padj in filtered_genes:
                file.write(f"{gene_id}\t{log2_fold_change}\t{padj}\n")

        print(f"Archivo guardado: {output_file}")
        print(f"Total de genes escritos: {len(filtered_genes)}")

    except Exception as e:
        print(f"Error al escribir el archivo: {e}")


# Función para crear resumen
def print_summary(filtered_genes):
    """Imprime un resumen de los genes filtrados."""

    # Contar genes por clasificación
    upregulated_count = sum(1 for gene in filtered_genes if gene[3] == "upregulated")
    downregulated_count = sum(
        1 for gene in filtered_genes if gene[3] == "downregulated"
    )
    total_significant = len(filtered_genes)

    # Imprimir resumen
    print("RESUMEN DE ANÁLISIS DE EXPRESIÓN GÉNICA")
    print(f"Genes significativos: {total_significant}")
    print(f"Genes sobreexpresados: {upregulated_count}")
    print(f"Genes subexpresados: {downregulated_count}")


# Función principal que cordinara las demas funciones
def main():
    """Función principal para cordinar las demas funciones."""

    # Leer argumentos del usuario
    if len(sys.argv) < 3:
        print("Error: Argumentos insuficientes")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # 2. Definir umbrales
    try:
        lfc_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
        padj_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    except ValueError:
        print("Error: Los umbrales deben ser números")
        sys.exit(1)

    try:
        # Llamar a load_deseq2_results()
        genes = load_deseq2_results(input_file)

        if not genes:
            print("No se cargaron genes.")
            sys.exit(1)

        # Llamar a filter_genes()
        filtered_genes = filter_genes(genes, lfc_threshold, padj_threshold)

        if not filtered_genes:
            print(
                "No se encontraron genes significativos con los umbrales especificados."
            )
            sys.exit(1)

        # Llamar a write_results()
        write_results(filtered_genes, output_file)

        # Llamar a print_summary()
        print_summary(filtered_genes)

    except Exception as e:
        print(f"Error general: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
