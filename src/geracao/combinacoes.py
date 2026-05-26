"""Combinação vetorizada de pares de fragmentos conectados por pontes."""

import numpy as np
import pandas as pd
from rdkit import Chem


def _substituir_numeracao_anel(smiles: str, de: str = "1", para: str = "2") -> str:
    """Substitui o índice de fechamento de anel para evitar conflitos de numeração."""
    return smiles.replace(de, para)


def combinar_fragmentos_via_pontes(
    fragmentos_df: pd.DataFrame,
    pontes: dict[int, str],
    n_combinacoes: int = 10_000,
    tamanho_lote: int = 200_000,
    coluna_smiles: str = "smiles",
) -> pd.DataFrame:
    """
    Combina aleatoriamente pares de fragmentos com pontes intermediárias.

    Usa geração vetorizada (NumPy) e valida cada candidato com RDKit.
    Continua gerando lotes até atingir `n_combinacoes` moléculas únicas.

    Parâmetros
    ----------
    fragmentos_df : pd.DataFrame
        DataFrame com coluna de SMILES dos fragmentos.
    pontes : dict[int, str]
        Dicionário de grupos de ligação.
    n_combinacoes : int
        Número alvo de moléculas únicas válidas.
    tamanho_lote : int
        Candidatos gerados por iteração (para eficiência vetorial).
    coluna_smiles : str
        Nome da coluna de SMILES no DataFrame.

    Retorna
    -------
    pd.DataFrame com coluna 'smiles'.
    """
    smiles_grupo_1 = fragmentos_df[coluna_smiles].to_numpy(dtype=object)
    smiles_grupo_2 = np.vectorize(_substituir_numeracao_anel)(smiles_grupo_1)
    array_pontes = np.array(list(pontes.values()), dtype=object)

    n1 = len(smiles_grupo_1)
    n2 = len(smiles_grupo_2)
    n_pontes = len(array_pontes)

    vistos: set[str] = set()
    resultados: list[str] = []

    while len(vistos) < n_combinacoes:
        idx1 = np.random.randint(0, n1, tamanho_lote)
        idx2 = np.random.randint(0, n2, tamanho_lote)
        idx_p = np.random.randint(0, n_pontes, tamanho_lote)

        candidatos = smiles_grupo_1[idx1] + array_pontes[idx_p] + smiles_grupo_2[idx2]

        for smiles in candidatos:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            smiles_canonico = Chem.MolToSmiles(mol, canonical=False)
            if smiles_canonico not in vistos:
                vistos.add(smiles_canonico)
                resultados.append(smiles_canonico)

            if len(vistos) >= n_combinacoes:
                break

    return pd.DataFrame({"smiles": resultados})
