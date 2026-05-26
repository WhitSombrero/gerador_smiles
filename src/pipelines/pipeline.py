"""Orquestração do pipeline completo de geração de SMILES."""
import os
import pandas as pd

from fragmentos.aneis import ANEIS_6_MEMBROS
from fragmentos.auxocromos import  AUXOCROMOS
from fragmentos.pontes import PONTES

from geracao.combinacoes import  combinar_fragmentos_via_pontes
from geracao.fragmentos import gerar_fragmentos_com_auxocromos

from utils.aromatizacao import aplicar_aromatizacao_e_filtrar
from utils.validacao import     filtrar_smiles_validos, remover_duplicatas_moleculares


def executar_pipeline(
    caminho_saida: str = "data/banco_smiles_aromatizados.parquet",
    repeticoes_por_anel: int = 200,
    n_combinacoes: int = 10_000,
) -> pd.DataFrame:
    """
    Executa o pipeline completo de geração de SMILES.

    Etapas
    ------
    1. Gera fragmentos (anéis + auxocromos aleatórios)
    2. Filtra SMILES inválidos
    3. Remove duplicatas moleculares
    4. Combina pares de fragmentos via pontes
    5. Aromatiza anéis de 6 membros e filtra inválidos
    6. Exporta para Parquet

    Parâmetros
    ----------
    caminho_saida : str
        Caminho do arquivo Parquet de saída.
    repeticoes_por_anel : int
        Variantes geradas por anel na Etapa 1.
    n_combinacoes : int
        Moléculas combinadas alvo na Etapa 4.

    Retorna
    -------
    pd.DataFrame com coluna 'smiles'.
    """
    # Etapa 1 — Geração de fragmentos
    print("Etapa 1: Gerando fragmentos com auxocromos...")
    smiles_gerados = gerar_fragmentos_com_auxocromos(
        aneis=ANEIS_6_MEMBROS,
        auxocromos=AUXOCROMOS,
        repeticoes=repeticoes_por_anel,
    )

    # Etapa 2 — Filtragem de SMILES inválidos
    df_fragmentos = filtrar_smiles_validos(smiles_gerados)
    print(f"  → {len(df_fragmentos)} fragmentos válidos gerados.")

    # Etapa 3 — Deduplicação
    print("Etapa 3: Removendo duplicatas...")
    df_fragmentos = remover_duplicatas_moleculares(df_fragmentos)
    print(f"  → {len(df_fragmentos)} fragmentos únicos após deduplicação.")

    # Etapa 4 — Combinação via pontes
    print(f"Etapa 4: Combinando fragmentos via pontes (alvo: {n_combinacoes})...")
    df_combinados = combinar_fragmentos_via_pontes(
        fragmentos_df=df_fragmentos,
        pontes=PONTES,
        n_combinacoes=n_combinacoes,
    )
    print(f"  → {len(df_combinados)} combinações geradas.")

    # Etapa 5 — Aromatização e filtragem final
    print("Etapa 5: Aromatizando anéis de 6 membros e filtrando...")
    df_final = aplicar_aromatizacao_e_filtrar(df_combinados)

    # Etapa 6 — Exportação
    os.makedirs(os.path.dirname(caminho_saida), exist_ok=True)
    df_final.to_parquet(caminho_saida, index=False)
    print(f"Etapa 6: Exportado → '{caminho_saida}' ({len(df_final)} moléculas).")

    return df_final
