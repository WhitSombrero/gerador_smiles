"""Aromatização de anéis de 6 membros em moléculas geradas."""

import pandas as pd
from rdkit import Chem


def _aromatizar_aneis_6_membros(smiles: str) -> str | None:
    """
    Tenta aromatizar todos os anéis de 6 membros de uma molécula.

    Para cada anel de 6 membros encontrado, marca as ligações como
    aromáticas e chama SanitizeMol para validar a aromaticidade.

    Parâmetros
    ----------
    smiles : str
        SMILES da molécula a ser aromatizada.

    Retorna
    -------
    str com SMILES aromatizado, ou None se a aromatização falhar.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    for anel in mol.GetRingInfo().AtomRings():
        if len(anel) != 6:
            continue

        for i in range(len(anel)):
            atomo1 = anel[i]
            atomo2 = anel[(i + 1) % len(anel)]
            ligacao = mol.GetBondBetweenAtoms(atomo1, atomo2)
            if ligacao is None:
                return None
            ligacao.SetBondType(Chem.rdchem.BondType.AROMATIC)

        for idx_atomo in anel:
            mol.GetAtomWithIdx(idx_atomo).SetIsAromatic(True)

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None

    return Chem.MolToSmiles(mol)


def aplicar_aromatizacao_e_filtrar(
    df: pd.DataFrame,
    coluna_smiles: str = "smiles",
) -> pd.DataFrame:
    """
    Aplica aromatização e remove moléculas inválidas.

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com coluna de SMILES.
    coluna_smiles : str
        Nome da coluna de entrada.

    Retorna
    -------
    pd.DataFrame com coluna 'smiles' contendo apenas moléculas aromatizáveis.
    """
    df = df.copy()
    df["smiles_aromatico"] = df[coluna_smiles].apply(_aromatizar_aneis_6_membros)

    df_valido = df.dropna(subset=["smiles_aromatico"]).reset_index(drop=True)
    print(
        f"[aromatização] {len(df_valido)} / {len(df)} "
        "moléculas válidas após aromatização."
    )

    return df_valido[["smiles_aromatico"]].rename(columns={"smiles_aromatico": "smiles"})
