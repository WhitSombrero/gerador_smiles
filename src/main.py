"""
Ponto de entrada do pipeline de geração de SMILES.

Uso:
    python main.py
"""

from pipelines.pipeline import executar_pipeline

if __name__ == "__main__":
    executar_pipeline(
        caminho_saida="data/banco_smiles_aromatizados.parquet",
        repeticoes_por_anel=200,
        n_combinacoes=10_000,
    )
