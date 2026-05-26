"""Validação de SMILES e remoção de duplicatas moleculares via RDKit."""

import pandas as pd
from rdkit import Chem


def filtrar_smiles_validos(smiles_dict: dict[int, str]) -> pd.DataFrame:
    """
    Converte dicionário {id: smiles} em DataFrame mantendo apenas SMILES válidos.

    Parâmetros
    ----------
    smiles_dict : dict[int, str]
        Dicionário gerado pelas funções de geração de fragmentos.

    Retorna
    -------
    pd.DataFrame com colunas 'id' e 'smiles'.
    """
    registros = [
        {"id": idx, "smiles": smiles}
        for idx, smiles in smiles_dict.items()
        if Chem.MolFromSmiles(smiles) is not None
    ]
    return pd.DataFrame(registros)


def remover_duplicatas_moleculares(
    df: pd.DataFrame,
    coluna_smiles: str = "smiles",
) -> pd.DataFrame:
    """
    Remove moléculas duplicadas por SMILES canônico (RDKit).

    Preserva a primeira ocorrência de cada molécula única e descarta
    SMILES inválidos silenciosamente.

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com coluna de SMILES.
    coluna_smiles : str
        Nome da coluna que contém os SMILES.

    Retorna
    -------
    pd.DataFrame sem duplicatas, com índice resetado.
    """
    vistos: set[str] = set()
    indices_unicos: list[int] = []

    for idx, smiles in df[coluna_smiles].items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        smiles_canonico = Chem.MolToSmiles(mol, canonical=True)
        if smiles_canonico not in vistos:
            vistos.add(smiles_canonico)
            indices_unicos.append(idx)

    return df.loc[indices_unicos].reset_index(drop=True)
